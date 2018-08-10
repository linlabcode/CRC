"""SET OF GENERAL UTILITY FUNCTIONS FOR SEQ DATA."""

# Functions require samtools to be callable with a 'samtools' command and bamliquidator to be
# callable with a 'bamliquidator' command

import gzip
import os
import subprocess
from collections import defaultdict

# ==================================================================
# ==========================I/O FUNCTIONS===========================
# ==================================================================


def open_plus(file_name, mode='r'):
    """Open function that can handle gzip files."""
    if file_name.split('.')[-1] == 'gz':
        return gzip.open(file_name, mode + 'b')
    return open(file_name, mode)


def parse_table(file_name, sep, header=False, excel=False):
    """Opens standard delimited files and outputs a nested list."""
    with open_plus(file_name) as infile:
        line = infile.readline()

        # First figure out if this table came from excel and has \r for line breaks
        if line.count('\r') > 0:
            excel = True

    with open_plus(file_name) as infile:
        if header:
            header = infile.readline()  # disposes of the header

        table = []
        if excel:
            lines_raw = infile.readlines()
            for line in lines_raw:
                line = line.replace('\r', '')
                line = line.rstrip()
                table.append(line.split(sep))
        else:
            for line in infile:
                line = line.rstrip().split(sep)
                table.append(line)

        return table


def unparse_table(table, output, sep):
    """Transforms a 'parse_table' output into a text file."""
    with open_plus(output, 'w') as fh_out:
        if not sep:
            for i in table:
                fh_out.write(str(i) + '\n')
        else:
            for line in table:
                fh_out.write(sep.join([str(x) for x in line]) + '\n')


def format_folder(folder_name, create=False):
    """Checks if a folder exists and returns a bool.

    If 'create=True' it creates a folder if it doesn't exist already.

    """
    if folder_name[-1] != '/':
        folder_name += '/'

    try:
        os.listdir(folder_name)
        return folder_name
    except OSError:
        print('Folder {} does not exist'.format(folder_name))
        if create:
            os.mkdir(folder_name)
            return folder_name
    return False


# ==================================================================
# ===================ANNOTATION FUNCTIONS===========================
# ==================================================================

def import_bound_region(bound_region_file, name):
    """Imports bound regions in either bed format or in error model format."""
    bound = parse_table(bound_region_file, '\t')
    loci_list = []
    ticker = 1
    bed = bool(bound_region_file.split('.')[-1] == 'bed')
    if bed:
        for line in bound:
            if len(line) < 3:
                continue
            if ticker % 1000 == 0:
                print(ticker)
            loci_list.append(
                Locus(
                    line[0],
                    int(line[1]),
                    int(line[2]),
                    '.',
                    ID=name + '_' + str(ticker),
                )
            )
            ticker = ticker + 1
    else:
        for line in bound:
            if ticker % 1000 == 0:
                print(ticker)
            loci_list.append(
                Locus(
                    'chr' + line[0],
                    int(line[1]),
                    int(line[2]),
                    '.',
                    ID=name + '_' + str(ticker),
                )
            )
            ticker = ticker + 1

    return LocusCollection(loci_list, 500)


# ==================================================================
# ========================LOCUS INSTANCE============================
# ==================================================================

# Locus and LocusCollection instances courtesy of Graham Ruby

class Locus(object):
    """Standard locus class for tracking genomic loci."""
    __chr_dict = dict()
    __sense_dict = {'+': '+', '-': '-', '.': '.'}

    def __init__(self, chr, start, end, sense, ID='', score=0):
        """Initialize attributes."""
        coords = sorted([int(start), int(end)])
        # This method for assigning chromosome should help avoid storage of redundant strings.
        if chr not in self.__chr_dict:
            self.__chr_dict[chr] = chr
        self.chr = self.__chr_dict[chr]
        self.sense = self.__sense_dict[sense]
        self.start = int(coords[0])
        self.end = int(coords[1])
        self.id = ID
        self.score = score

    def len(self):
        """get locus length."""
        return self.end - self.start + 1

    def get_antisense_locus(self):
        """Get antisense locus."""
        if self.sense == '.':
            return self
        switch = {'+': '-', '-': '+'}
        return Locus(self.chr, self.start, self.end, switch[self.sense])

    def coords(self):
        """Returns a sorted list of the coordinates."""
        return [self.start, self.end]

    def overlaps(self, other_locus):
        """Find if loci overlap.

        True if two loci share any coordinates in common.

        """
        if self.chr != other_locus.chr:
            return False
        if not(self.sense == '.' or
               other_locus.sense == '.' or
               self.sense == other_locus.sense):
            return False
        if self.start > other_locus.end or other_locus.start > self.end:
            return False
        return True

    def contains(self, other_locus):
        """Find if locus contains other locus.

        True if all the nucleotides of the given locus overlap with the self locus.

        """
        if self.chr != other_locus.chr:
            return False
        if not(self.sense == '.' or
               other_locus.sense == '.' or
               self.sense == other_locus.sense):
            return False
        if self.start > other_locus.start or other_locus.end > self.end:
            return False
        return True

    def overlaps_antisense(self, other_locus):
        """Same as overlaps, but considers the opposite strand."""
        return self.get_antisense_locus().overlaps(other_locus)

    def contains_antisense(self, other_locus):
        """Same as contains, but considers the opposite strand."""
        return self.get_antisense_locus().contains(other_locus)

    def __hash__(self):
        """Define hash."""
        return self.start + self.end

    def __eq__(self, other):
        """Define Locus objects equality."""
        if self.__class__ != other.__class__:
            return False
        if self.chr != other.chr:
            return False
        if self.start != other.start:
            return False
        if self.end != other.end:
            return False
        if self.sense != other.sense:
            return False
        return True

    def __ne__(self, other):
        """Opposite to equality."""
        return not self.__eq__(other)

    def __str__(self):
        """Locus string representation."""
        return self.chr + '(' + self.sense + '):' + '-'.join(map(str, self.coords()))


class LocusCollection(object):
    """A collection of locus objects used for querying large sets of loci."""
    def __init__(self, loci, window_size=50):
        """Initialize attributes."""
        # Top-level keys are chr, then strand, no space
        self.__chr_to_coord_to_loci = dict()
        self.__loci = dict()
        self.__win_size = window_size
        for lcs in loci:
            self.__add_locus(lcs)

    def __add_locus(self, lcs):
        """Add locus to the collection."""
        if lcs not in self.__loci:
            self.__loci[lcs] = None
            if lcs.sense == '.':
                chr_key_list = [lcs.chr + '+', lcs.chr + '-']
            else:
                chr_key_list = [lcs.chr + lcs.sense]
            for chr_key in chr_key_list:
                if chr_key not in self.__chr_to_coord_to_loci:
                    self.__chr_to_coord_to_loci[chr_key] = dict()
                for num in self.__get_key_range(lcs):
                    if num not in self.__chr_to_coord_to_loci[chr_key]:
                        self.__chr_to_coord_to_loci[chr_key][num] = []
                    self.__chr_to_coord_to_loci[chr_key][num].append(lcs)

    def __get_key_range(self, locus):
        """Get key range."""
        start = locus.start / self.__win_size
        end = locus.end / self.__win_size + 1  # add 1 because of the range
        return range(int(start), int(end))

    def __len__(self):
        """Define length."""
        return len(self.__loci)

    def append(self, new):
        """Append locus."""
        self.__add_locus(new)

    def has_locus(self, locus):
        """Check if locus in the loci dict"""
        return locus in self.__loci

    def remove(self, old):
        """Remove locus from the collection."""
        if old not in self.__loci:
            raise ValueError("Requested locus isn't in collection.")
        del self.__loci[old]
        if old.sense == '.':
            sense_list = ['+', '-']
        else:
            sense_list = [old.sense]
        for k in self.__get_key_range(old):
            for sense in sense_list:
                self.__chr_to_coord_to_loci[old.chr + sense][k].remove(old)

    def get_loci(self):
        """Get loci."""
        return self.__loci.keys()

    def __subset_helper(self, locus, sense):
        """"""
        sense = sense.lower()
        if ['sense', 'antisense', 'both'].count(sense) != 1:
            raise ValueError("Sense command invalid: '{}'.".format(sense))
        matches = dict()
        senses = ['+', '-']
        if locus.sense == '.' or sense == 'both':
            def lamb(s): return True
        elif sense == 'sense':
            def lamb(s): return bool(s == locus.sense)
        elif sense == 'antisense':
            def lamb(s): return bool(s != locus.sense)
        else:
            raise ValueError("Sense value was inappropriate: '{}'.".format(sense))

        for s in filter(lamb, senses):
            chr_key = locus.chr + s
            if chr_key in self.__chr_to_coord_to_loci:
                for n in self.__get_key_range(locus):
                    if n in self.__chr_to_coord_to_loci[chr_key]:
                        for lcs in self.__chr_to_coord_to_loci[chr_key][n]:
                            matches[lcs] = None

        return matches.keys()

    def get_overlap(self, locus, sense='sense'):
        """Returns all members of the collection that overlap the locus.

        :param sense: can be 'sense' (default), 'antisense', or 'both'.

        """
        matches = self.__subset_helper(locus, sense)
        # Now, get rid of the ones that don't really overlap
        real_matches = dict()
        if sense == 'sense' or sense == 'both':
            for i in filter(lambda lcs: lcs.overlaps(locus), matches):
                real_matches[i] = None
        if sense == 'antisense' or sense == 'both':
            for i in filter(lambda lcs: lcs.overlaps_antisense(locus), matches):
                real_matches[i] = None
        return real_matches.keys()

    def stitch_collection(self, stitch_window=1, sense='both'):
        """Reduces the collection by stitching together overlapping loci.

        Returns a new collection. Initializing stichWindow to 1 helps collect directly adjacent
        loci.

        :param sense: can be 'sense' (default), 'antisense', or 'both'.

        """
        locus_list = self.get_loci()
        old_collection = LocusCollection(locus_list, 500)

        stitched_collection = LocusCollection([], 500)

        for locus in locus_list:
            if old_collection.has_locus(locus):
                old_collection.remove(locus)
                overlapping_loci = old_collection.get_overlap(
                    Locus(
                        locus.chr,
                        locus.start - stitch_window,
                        locus.end + stitch_window,
                        locus.sense,
                        locus.id,
                    ),
                    sense,
                )

                stitch_ticker = 1
                while overlapping_loci:
                    stitch_ticker += len(overlapping_loci)
                    overlap_coords = locus.coords()

                    for overlapping_locus in overlapping_loci:
                        overlap_coords += overlapping_locus.coords()
                        old_collection.remove(overlapping_locus)
                    overlap_coords = [int(x) for x in overlap_coords]
                    if sense == 'both':
                        locus = Locus(
                            locus.chr,
                            min(overlap_coords),
                            max(overlap_coords),
                            '.',
                            locus.id,
                        )
                    else:
                        locus = Locus(
                            locus.chr,
                            min(overlap_coords),
                            max(overlap_coords),
                            locus.sense,
                            locus.id,
                        )
                    overlapping_loci = old_collection.get_overlap(
                        Locus(
                            locus.chr,
                            locus.start - stitch_window,
                            locus.end + stitch_window,
                            locus.sense,
                        ),
                        sense
                    )
                locus.id = '{}_{}_lociStitched'.format(stitch_ticker, locus.id)

                stitched_collection.append(locus)

            else:
                continue

        return stitched_collection


# ==================================================================
# ========================LOCUS FUNCTIONS===========================
# ==================================================================

def make_search_locus(locus, up_search, down_search):
    """Takes a locus and expands it by a fixed upstream/downstream amount.

    Spits out the new larger locus.

    """
    if locus.sense == '-':
        search_locus = Locus(
            locus.chr,
            locus.start - down_search,
            locus.end + up_search,
            locus.sense,
            locus.id,
        )
    else:
        search_locus = Locus(
            locus.chr,
            locus.start - up_search,
            locus.end + down_search,
            locus.sense,
            locus.id,
        )
    return search_locus


# ==================================================================
# ==========================BAM CLASS===============================
# ==================================================================

def convert_bitwise_flag(flag):
    """Get strand from the flag."""
    if flag & 16:
        return "-"
    return "+"


class Bam(object):
    """A class for a sorted and indexed bam file that allows easy analysis of reads."""
    def __init__(self, bam_file):
        """Initialize attributes."""
        self._bam = bam_file

        # Get the number of reads use idxstats instead of flagstat
        command = '{} idxstats {}'.format('samtools', self._bam)

        stats = subprocess.Popen(
            command,
            stdin=subprocess.PIPE,
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
            shell=True,
        )
        stat_lines = stats.stdout.readlines()
        stats.stdout.close()

        self._mapped_reads = sum(
            [int(line.decode('utf-8').rstrip().split('\t')[-2]) for line in stat_lines]
        )
        self._total_reads = (
            self._mapped_reads +
            int(stat_lines[-1].decode('utf-8').rstrip().split('\t')[-1])
        )

        # Now get the readlength #check the first 1000 reads
        view_command = '{} view {} chr1:90000000-100000000'.format('samtools', self._bam)
        read_stats = subprocess.Popen(
            view_command,
            stdin=subprocess.PIPE,
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
            shell=True,
        )
        read_stat_lines = read_stats.stdout.readlines()
        read_stats.stdout.close()

        if not read_stat_lines:  # try chrom names w/o the chr
            view_command = '{} view {} 1:90000000-100000000'.format('samtools', self._bam)
            read_stats = subprocess.Popen(
                view_command,
                stdin=subprocess.PIPE,
                stderr=subprocess.PIPE,
                stdout=subprocess.PIPE,
                shell=True,
            )
            read_stat_lines = read_stats.stdout.readlines()
            read_stats.stdout.close()

        self._read_lengths = uniquify([
            len(line.decode('utf-8').split('\t')[9]) for line in read_stat_lines if len(line) > 0
        ])

    def path(self):
        """Return bam file path."""
        return self._bam

    def get_read_lengths(self):
        """Get read lengths."""
        return self._read_lengths

    def get_raw_reads(self, locus, sense, unique=False, include_jxn_reads=False,
                      print_command=False):
        """Gets raw reads from the bam using samtools view.

        Can enforce uniqueness and strandedness.

        """
        locus_line = locus.chr + ':' + str(locus.start) + '-' + str(locus.end)

        command = '{} view {} {}'.format('samtools', self._bam, locus_line)
        if print_command:
            print(command)
        get_reads = subprocess.Popen(
            command,
            stdin=subprocess.PIPE,
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
            shell=True,
        )
        reads = get_reads.communicate()
        reads = reads[0].decode('utf-8').split('\n')[:-1]
        reads = [read.split('\t') for read in reads]
        if not include_jxn_reads:
            reads = filter(lambda x: x[5].count('N') < 1, reads)

        kept_reads = []
        seq_dict = defaultdict(int)
        if sense == '-':
            strand = ['+', '-']
            strand.remove(locus.sense)
            strand = strand[0]
        else:
            strand = locus.sense
        for read in reads:
            read_strand = convert_bitwise_flag(read[1])

            if sense == 'both' or sense == '.' or read_strand == strand:
                if unique and seq_dict[read[9]] == 0:
                    kept_reads.append(read)
                elif not unique:
                    kept_reads.append(read)
            seq_dict[read[9]] += 1

        return kept_reads

    def liquidate_locus(self, locus, n=1, sense='.', extension=0, mmr=False):
        """Use bamliquidator on a locus."""
        bamliquidator_cmd = (
            'bamliquidator {} {} {} {} {} {} {}'.format(
                self._bam,
                locus.chr,
                str(locus.start),
                str(locus.end),
                sense,
                n,
                extension,
            )
        )
        bamliquidator_out = subprocess.Popen(
            bamliquidator_cmd,
            stdout=subprocess.PIPE,
            shell=True,
        )
        score = [int(x.decode('utf-8').rstrip()) for x in bamliquidator_out.stdout.readlines()]
        if mmr:
            mmr_count = float(self._mapped_reads) / 1000000
            score = [round(x/mmr_count, 4) for x in score]

        return score


# ==================================================================
# ========================MISC FUNCTIONS============================
# ==================================================================

# uniquify function
# by Peter Bengtsson
# Used under a creative commons license
# sourced from  here: http://www.peterbe.com/plog/uniqifiers-benchmark

def uniquify(seq, idfun=None):
    """Makes a list unique.

    Order preserving.

    """
    seen = {}
    result = []
    for item in seq:
        if idfun is None:
            marker = item
        else:
            marker = idfun(item)
        if marker in seen:
            continue
        seen[marker] = 1
        result.append(item)
    return result


# ==================================================================
# =======================SEQUENCE FUNCTIONS=========================
# ==================================================================

def fetch_seq(directory, chrom, start, end, UCSC=False, line_breaks=True, header=True):
    """Fetches a sequence from a genome directory.

    The directory must contain individual chrom fasta files.

    """
    file_name = directory + chrom + '.fa'
    header_offset = 0
    n_start = 0
    n_end = 0
    with open_plus(file_name, 'r') as infile:
        if header:
            infile.seek(0)
            header_offset = len(infile.readline())
        if line_breaks:
            n_start = (start - 1) // 50
            n_end = (end - 1) // 50
        if UCSC:
            infile.seek((start + n_start + header_offset))
        else:
            infile.seek((start - 1 + n_start + header_offset))
        span = ((end + n_end - 1) - (start + n_start - 1))
        read = infile.read(span)
        if line_breaks:
            read = read.replace('\n', '')
        return read
