"""CRC analysis utilities."""

# Functions require fasta-get-markov to be callable with a 'fasta-get-markov' command
# and fimo with 'fimo' command. Both are part of the MEME suite.

import pickle
import subprocess
import sys
from collections import defaultdict

import networkx as nx
import numpy
from crc.definitions import ROOT_DIR
from crc.resources import utils
from networkx.algorithms.clique import find_cliques_recursive

# ================================================================================
# ===================================CLASSES======================================
# ================================================================================


class Genome(object):
    """Genome has attributes of a build name, a fasta directory and an annotation file."""
    __chrDict = dict()
    __featureDict = dict()

    def __init__(self, name, genome_directory):
        """Initialize attributes."""
        self._name = name
        self._directory = genome_directory

    def name(self):
        """Name of the genome."""
        return self._name

    def directory(self):
        """Genome directory."""
        return self._directory

    def add_feature(self, feature, path):
        """Add a genome feature to the feature dict."""
        if feature in self.__featureDict:
            print('WARNING OVERRIDING {} PATH WITH {}'.format(feature, path))
        self.__featureDict[feature] = path

    def return_feature(self, feature):
        """Tries to load the selected feature from the feature dictionary."""
        if feature not in self.__featureDict:
            print('ERROR: GENOME {} DOES NOT HAVE FEATURE {}'.format(self.name(), feature))
            sys.exit()
        else:
            return self.__featureDict[feature]

    def has_feature(self, feature):
        """Check if the feature exists."""
        return bool(feature in self.__featureDict)


# ================================================================================
# =================================FUNCTIONS======================================
# ================================================================================

def load_genome(genome_build, chrom_path, mask_file=None, config_file=''):
    """Loads annotation for a genome into a genome object."""

    # This nested dictionary has all of the useful information and likely will have to be
    # edited so it can be configured any time
    genome_dict = {
        'HG19': {
            'tf_file': '{}/annotation/TFlist_NMid_hg19.txt'.format(ROOT_DIR),
            'mask': '{}/annotation/hg19_encode_blacklist.bed'.format(ROOT_DIR),
            'motif_convert': '{}/annotation/MotifDictionary.txt'.format(ROOT_DIR),
            'motif_database': '{}/annotation/VertebratePWMs.txt'.format(ROOT_DIR),
        },
        'RN6': {
            'tf_file': '{}/annotation/TFlist_NMid_rn6.txt'.format(ROOT_DIR),
            'motif_convert': '{}/annotation/MotifDictionary.txt'.format(ROOT_DIR),
            'motif_database': '{}/annotation/VertebratePWMs.txt'.format(ROOT_DIR),
        },
        'MM10': {
            'tf_file': '{}/annotation/TFlist_NMid_mm10.txt'.format(ROOT_DIR),
            'motif_convert': '{}/annotation/MotifDictionary.txt'.format(ROOT_DIR),
            'motif_database': '{}/annotation/VertebratePWMs.txt'.format(ROOT_DIR),
        },
    }

    genome_build = genome_build.upper()

    # Allow an optional config file to overwrite default paths
    if config_file:
        config_table = utils.parse_table(config_file, '\t')
        for line in config_table[1:]:
            (build, field, feature_path) = line[0].split(':')
            genome_dict[build.upper()][field.lower()] = feature_path

    if genome_build not in genome_dict:
        print('ERROR: UNSUPPORTED GENOME BUILD {}. EXITING NOW'.format(genome_build))
        sys.exit()
    else:
        print('USING BUILD {} WITH FOLLOWING FIELDS:'.format(genome_build))
        print(genome_dict[genome_build])

    # Now attempt to load the genome
    genome = Genome(genome_build, chrom_path)

    # Adding additional optional features
    genome.add_feature('tf_file', genome_dict[genome_build]['tf_file'])
    genome.add_feature('motif_convert', genome_dict[genome_build]['motif_convert'])
    genome.add_feature('motif_database', genome_dict[genome_build]['motif_database'])
    if mask_file:
        genome.add_feature('mask', mask_file)
    else:
        genome.add_feature('mask', genome_dict[genome_build]['mask'])

    return genome


def gene_to_enhancer(genome, enhancer_file, activity_path):
    """Assign each Super-Enhancer to the closest active TSS to its center.

    Return a dictionary keyed by TF that points to a list of loci.

    """
    print('Identifying enhancers and target genes from {}'.format(enhancer_file))
    # Should this do gene assignment????
    # For now assume gene assignment has been done
    # Can later toggle to do gene assignment

    # First load the TF lists
    tf_table = utils.parse_table(genome.return_feature('tf_file'), '\t')

    motif_table = utils.parse_table(genome.return_feature('motif_convert'), '\t')

    # This gives all tfs that have a motif
    motif_tfs = utils.uniquify([line[1] for line in motif_table])

    # Intersect w/ the activity table
    if activity_path:
        activity_table = utils.parse_table(activity_path, '\t')

        # Figure out the right column for actual gene names
        # (basically not NM or NR and not a numeral)
        for i in range(len(activity_table[0])):
            # Assumes refseq
            if (activity_table[0][i][0:2] != 'NM' and activity_table[0][i][0:2] != 'NR'
                    and not activity_table[0][i].isdigit()):
                gene_col = i
                break
        print(
            'using column {} of {} gene activity table for common names'
            ''.format(gene_col + 1, activity_path)
        )

        active_gene_list = [line[gene_col].upper() for line in activity_table]

        tf_list_name = utils.uniquify(
            [line[1] for line in tf_table
             if active_gene_list.count(line[1]) > 0 and motif_tfs.count(line[1]) > 0]
        )
    else:
        tf_list_name = [line[1] for line in tf_table if motif_tfs.count(line[1]) > 0]

    print(
        'Identified {} TFs from {} that have motifs'
        ''.format(len(tf_list_name), genome.return_feature('tf_file'))
    )

    # Keyed by gene with loci objects in the list
    gene_to_enhancer_dict = defaultdict(list)
    enhancer_to_gene_dict = defaultdict(list)

    # Assuming id,chrom,start,stop w/ gene names in the last 3 columns per standard ROSE output
    enhancer_table = utils.parse_table(enhancer_file, '\t')
    print('Analyzing {} cis-regulatory regions'.format(len(enhancer_table)))

    # Now let's make the enhancer table by region and then by gene
    enhancer_region_table = [['ENHANCER_ID', 'CHROM', 'START', 'STOP', 'GENE_LIST']]
    enhancer_tf_region_table = [['ENHANCER_ID', 'CHROM', 'START', 'STOP', 'GENE_LIST']]
    gene_region_table = [['GENE', 'TF', 'CHROM', 'START', 'STOP', 'ENHANCER_ID']]
    gene_tf_region_table = [['GENE', 'CHROM', 'START', 'STOP', 'ENHANCER_ID']]
    gene_summary_table = [['GENE', 'TF', 'ENHANCER_LIST']]

    # Will need to track which ones are TFs
    candidate_tf_list = []
    # Find the columns for gene assignment
    header = enhancer_table[0]
    header_length = len(enhancer_table[0])
    closest_index = header.index('CLOSEST_GENE')
    proximal_index = header.index('PROXIMAL_GENES')
    overlap_index = header.index('OVERLAP_GENES')
    for line in enhancer_table[1:]:
        # Don't bother trying to figure out lines w/o target genes
        if len(line) != header_length:
            continue
        enhancer_locus = utils.Locus(line[1], line[2], line[3], '.', line[0])
        closest_gene_list = line[closest_index].split(',') if line[closest_index] else []
        proximal_gene_list = line[proximal_index].split(',') if line[proximal_index] else []
        overlap_gene_list = line[overlap_index].split(',') if line[overlap_index] else []
        all_gene_list = closest_gene_list + proximal_gene_list + overlap_gene_list
        all_gene_list = [gene.upper() for gene in all_gene_list]

        # Gets a unique list of all tfs
        if activity_path:
            all_gene_list = utils.uniquify(
                [gene for gene in all_gene_list if active_gene_list.count(gene) > 0]
            )
        else:
            all_gene_list = utils.uniquify(all_gene_list)

        candidate_gene_list = utils.uniquify(
            [gene for gene in all_gene_list if tf_list_name.count(gene) > 0]
        )
        if all_gene_list:
            for gene in all_gene_list:
                gene_to_enhancer_dict[gene].append(enhancer_locus)
                enhancer_to_gene_dict[enhancer_locus].append(gene)
            newline = line[0:4] + [','.join(all_gene_list)]
        else:
            newline = line[0:4] + ['']
        enhancer_region_table.append(newline)

        if candidate_gene_list:
            tf_line = line[0:4] + [','.join(candidate_gene_list)]
            enhancer_tf_region_table.append(tf_line)

    # Now iterate through each gene and list the enhancers
    gene_list = list(gene_to_enhancer_dict.keys())
    print(gene_list)
    gene_list.sort()
    for gene in gene_list:
        if tf_list_name.count(gene) > 0:
            tf_status = 1
            candidate_tf_list.append(gene)
        else:
            tf_status = 0
        enhancer_loci = gene_to_enhancer_dict[gene]
        enhancer_string = ','.join([enhancer.id for enhancer in enhancer_loci])
        gene_summary_table.append([gene, tf_status, enhancer_string])
        for enhancer in enhancer_loci:
            newline = [
                gene,
                tf_status,
                enhancer.chr,
                enhancer.start,
                enhancer.end,
                enhancer.id,
            ]
            gene_region_table.append(newline)
            if tf_status == 1:
                newline = [gene, enhancer.chr, enhancer.start, enhancer.end, enhancer.id]
                gene_tf_region_table.append(newline)

    return (
        gene_region_table,
        gene_tf_region_table,
        enhancer_region_table,
        enhancer_tf_region_table,
        gene_summary_table,
        candidate_tf_list,
        gene_to_enhancer_dict,
    )


def gaussian_smooth(read_list, degree=5):
    """Smoothing function for raw bamliquidator output."""
    window = degree * 2 - 1
    weight = numpy.array([1.0] * window)
    weight_gauss = []

    for i in range(window):
        i = i - degree + 1
        frac = i / float(window)
        gauss = 1 / (numpy.exp((4 * (frac)) ** 2))
        weight_gauss.append(gauss)

    weight = numpy.array(weight_gauss) * weight
    smoothed = [0.0] * (len(read_list) - window)
    for i in range(len(smoothed)):
        smoothed[i] = sum(numpy.array(read_list[i:i+window]) * weight) / sum(weight)
    smoothed = [0, 0, 0, 0, 0] + smoothed + [0, 0, 0, 0]  # return an array of the same length

    return smoothed


def score_valley(locus, bam_list, max_read_length):
    """Calculate valley scores for a locus.

    Based on this refernce:
    http://bioinformatics.oxfordjournals.org/content/26/17/2071.full

    Takes in a bamDict where for each dataset you hold the path, the mmr, and the readlength so
    you can match extension pull in w subprocess.

    """
    # Average density for all in the group
    n_bins = locus.len() // 10

    # Want to make an average density gram
    density_matrix = []
    for bam in bam_list:
        # Calculate the extension
        extension = max_read_length - bam.get_read_lengths()[0]
        # This gives the normalized signal vector
        signal_vector = bam.liquidate_locus(locus, n_bins, '.', extension, mmr=True)
        density_matrix.append(signal_vector)

    # Now get the average
    if len(density_matrix) > 1:
        density = [
            numpy.average([line[i] for line in density_matrix])
            for i in range(len(density_matrix[0]))
        ]
    else:
        density = density_matrix[0]

    smooth_density = gaussian_smooth(density, 5)

    score_array = []
    region_max = max(smooth_density)

    # Now take the smooth reads and calculate a valley score
    for i, smooth_density_value in enumerate(smooth_density):
        score = 0
        try:
            leftmax = max(smooth_density[i-25:i-10])
        except ValueError:
            leftmax = 'edge'
        try:
            rightmax = max(smooth_density[i+10:i+25])
        except ValueError:
            rightmax = 'edge'

        if rightmax == 'edge' and leftmax == 'edge':
            shoulder_height_max = 0
        elif leftmax == 'edge':
            shoulder_height_max = rightmax
        elif rightmax == 'edge':
            shoulder_height_max = leftmax
        else:
            shoulder_height_max = max(leftmax, rightmax)

        ratio = (shoulder_height_max - float(smooth_density_value)) / region_max
        if ratio > 0.3:
            score = 1
        else:
            score = 0

        score_array.append(score)

    return score_array


def stitch_valleys(valley_list):
    """Returns a stitched list of valleys to extract seq from."""
    valley_collection = utils.LocusCollection(valley_list, 1)
    stitched_valley_collection = valley_collection.stitch_collection()
    loci = []
    regions = []
    for valley in stitched_valley_collection.get_loci():
        if [valley.chr, valley.start, valley.end] not in regions:
            loci.append(valley)
            regions.append([valley.chr, valley.start, valley.end])
    return loci


def find_valleys(gene_to_enhancer_dict, bam_file_list, project_name, project_folder, cutoff=0.2):
    """Returns a dictionary of refseqs with all valley loci that are associated.

    Returns 2 kinds of bed files. 1 = all

    """
    # First make the bamDict
    all_valley_bed = []
    valley_dict = {}

    # Start w/ a bam_file_list and make a list of bam type objects
    bam_list = [utils.Bam(bam_path) for bam_path in bam_file_list]
    max_read_length = max([bam.get_read_lengths()[0] for bam in bam_list])

    gene_list = list(gene_to_enhancer_dict.keys())
    gene_list.sort()
    ticker = 0
    print('number of regions processed:')
    for gene in gene_list:

        valley_dict[gene] = []

        for region in gene_to_enhancer_dict[gene]:
            if ticker % 100 == 0:
                print(ticker)
            ticker += 1
            score_array = score_valley(
                region,
                bam_list,
                max_read_length,
            )
            for index, score in enumerate(score_array):
                if score > cutoff:
                    valley = utils.Locus(region.chr, region.start + index * 10,
                                         region.start + (index + 1) * 10, '.')
                    valley_dict[gene].append(valley)

        stitched_valleys = stitch_valleys(valley_dict[gene])
        for valley in stitched_valleys:
            all_valley_bed.append([valley.chr, valley.start, valley.end])
            valley_dict[gene] = stitched_valleys

    all_bed_path = project_folder + project_name + '_all_valleys.bed'
    utils.unparse_table(all_valley_bed, all_bed_path, '\t')

    return all_bed_path


def filter_subpeaks(subpeak_file, analysis_name, output_folder):
    """Takes the initial subpeaks in, stitches them."""
    # Stitch the subpeaks
    print(subpeak_file)
    subpeak_collection = utils.import_bound_region(subpeak_file, '%s_subpeak' % (analysis_name))

    subpeak_collection = subpeak_collection.stitch_collection()

    subpeak_loci = subpeak_collection.get_loci()

    all_sub_bed = []
    for locus in subpeak_loci:
        bed_line = [locus.chr, locus.start, locus.end, '.', locus.id]
        all_sub_bed.append(bed_line)

    all_bed_path = output_folder + analysis_name + '_all_subpeak.bed'
    utils.unparse_table(all_sub_bed, all_bed_path, '\t')

    return all_bed_path


def generate_subpeak_fasta(gene_to_enhancer_dict, subpeaks, genome, project_name, const_extension):
    """Generate a subpeak FASTA.

    From a BED file of constituents generate a FASTA for the consituients contained within the
    canidate supers.

    """
    genome_directory = genome.directory()
    subpeak_dict = {}
    subpeak_bed = [['track name=' + project_name + ' color=204,0,204']]
    subpeak_table = utils.parse_table(subpeaks, '\t')

    subpeak_loci = [utils.Locus(l[0], int(l[1]), int(l[2]), '.') for l in subpeak_table]
    subpeak_collection = utils.LocusCollection(subpeak_loci, 50)

    for gene in gene_to_enhancer_dict.keys():
        subpeak_dict[gene] = []
        for region in gene_to_enhancer_dict[gene]:
            overlaps = subpeak_collection.get_overlap(region)
            extended_overlaps = [
                utils.make_search_locus(x, const_extension, const_extension) for x in overlaps
            ]

            overlap_collection_temp = utils.LocusCollection(extended_overlaps, 50)
            overlap_collection = overlap_collection_temp.stitch_collection()
            for overlap in overlap_collection.get_loci():
                subpeak_bed.append([overlap.chr, overlap.start, overlap.end])
                subpeak_dict[gene].append(overlap)

    fasta = []
    for gene in subpeak_dict:
        for subpeak in subpeak_dict[gene]:
            fasta_title = '|'.join([gene, subpeak.chr, str(subpeak.start), str(subpeak.end)])
            fasta_line = utils.fetch_seq(genome_directory, subpeak.chr, int(subpeak.start + 1),
                                         int(subpeak.end + 1))

            fasta.append('>' + fasta_title)
            fasta.append(fasta_line.upper())

    return subpeak_bed, fasta


def make_motif_background(subpeak_fasta, project_folder, project_name):
    """Makes a 1st order markov background file for fimo."""
    bg_cmd = ''.join([
        'fasta-get-markov -m 1 < ',
        subpeak_fasta,
        '  > ',
        project_folder,
        project_name,
        '_bg.meme',
    ])
    bg_path = '{}{}_bg.meme'.format(project_folder, project_name)
    subprocess.call(bg_cmd, shell=True)

    return bg_path


def find_motifs(subpeak_fasta, bg_path, candidate_tf_list, project_folder, analysis_name,
                motif_convert_file, motif_database_file):
    """Find motifs.

    Takes the refseq to subpeak seq dict and returns the networkx object with all connections.

    """
    fimo_folder = utils.format_folder(project_folder + 'FIMO/', True)
    subpeak_name = subpeak_fasta.split('/')[-1].split('.')[0]
    output = '{}{}_fimo.txt'.format(fimo_folder, subpeak_name)

    # Create a dictionary to call motif names keyed on gene names
    motif_database = utils.parse_table(motif_convert_file, '\t')
    motif_database_dict = {}  # create a dict keyed by TF with multiple motifs

    for line in motif_database:
        motif_database_dict[line[1]] = []
    for line in motif_database:
        motif_database_dict[line[1]].append(line[0])

    candidate_tf_list.sort()

    print(candidate_tf_list)

    # Now make a list of all motifs
    motif_list = []
    for tf in candidate_tf_list:
        motif_list += motif_database_dict[tf]

    motif_list = utils.uniquify(motif_list)

    fimo_bash_path = '{}{}_fimo.sh'.format(fimo_folder, analysis_name)
    fimo_bash = open(fimo_bash_path, 'w')
    fimo_bash.write('#!/usr/bin/bash\n\n')

    fimo_cmd = 'fimo'
    for motif in motif_list:
        fimo_cmd += ' --motif ' + "'{}'".format(str(motif))

    # fimo_cmd += ' --thresh 1e-5'  # if you want to increase stringency
    fimo_cmd += ' -verbosity 1'
    fimo_cmd += ' -text'
    fimo_cmd += ' -oc ' + project_folder + 'FIMO'
    fimo_cmd += ' --bgfile {}'.format(bg_path)
    fimo_cmd += ' ' + motif_database_file + ' '
    fimo_cmd += subpeak_fasta
    fimo_cmd += ' > ' + output
    print(fimo_cmd)
    fimo_bash.write(fimo_cmd)
    fimo_bash.close()

    subprocess.call(fimo_cmd, shell=True)  # will wait that fimo is done to go on

    return output

#  After fimo, we should collapse the motif edges into beds before going into networks


def collapse_fimo(fimo_output, candidate_tf_list, output_folder, analysis_name,
                  motif_convert_file):
    """Collapses motifs from fimo.

    For each source node (TF) and each target node (gene enhancer regions), collapse motif
    instances then spit out a ginormous set of beds and a single crazy collapsed bed.

    """
    # First build up the motif name conversion database
    motif_database = utils.parse_table(motif_convert_file, '\t')
    motif_database_dict = defaultdict(list)

    # The reverse of the other dict, from motif name to gene name
    # A motif can go to multiple genes
    for line in motif_database:
        motif_database_dict[line[0]].append(line[1])

    # Make the folder to store motif beds
    utils.format_folder('{}motif_beds/'.format(output_folder), True)

    edge_dict = {}

    # First layer are source nodes
    for tf in candidate_tf_list:
        edge_dict[tf] = defaultdict(list)
    # Next layer are target nodes which are derived from the fimo output

    fimo_table = utils.parse_table(fimo_output, '\t')
    print(fimo_output)

    # fimo sometimes puts the region in either the first or second column
    fimo_line = fimo_table[1]
    if fimo_line[1].count('|') > 0:
        region_index = 1
    else:
        region_index = 2
    print('USING COLUMN {} OF FIMO OUTPUT FOR REGION'.format(region_index))

    for line in fimo_table[1:]:
        source_tfs = motif_database_dict[line[0]]  # motifId
        for source in source_tfs:
            if candidate_tf_list.count(source) == 0:
                continue
            region = line[region_index].split('|')

            target = region[0]
            if region_index == 2:
                target_locus = utils.Locus(
                    region[1],
                    int(region[2]) + int(line[3]),
                    int(region[2]) + int(line[4]),
                    '.'
                )
            else:
                target_locus = utils.Locus(
                    region[1],
                    int(region[2]) + int(line[2]),
                    int(region[2]) + int(line[3]),
                    '.'
                )

            # What's missing here is the enhancer id of the target locus
            try:
                edge_dict[source][target].append(target_locus)
            except KeyError:
                print('This motif is not in the network')
                print(line)
                sys.exit()

    # Now we actually want to collapse this down in a meaningful way
    # Overlapping motifs count as a single binding site. This way a TF with tons of motifs
    # that finds the same site over and over again doesn't get over counted
    all_bed = []
    all_bed_path = '{}{}_all_motifs.bed'.format(output_folder, analysis_name)
    for tf in candidate_tf_list:
        print(tf)
        target_nodes = edge_dict[tf].keys()
        bed_header = [
            'track name = "{}" description="{} motifs in {}"'.format(tf, tf, analysis_name)
        ]
        all_bed.append(bed_header)
        target_bed = [bed_header]
        target_bed_path = '{}motif_beds/{}_motifs.bed'.format(output_folder, tf)
        for target in target_nodes:
            edge_collection = utils.LocusCollection(edge_dict[tf][target], 50)
            edge_collection = edge_collection.stitch_collection()
            edge_loci = edge_collection.get_loci()
            edge_dict[tf][target] = edge_loci
            for locus in edge_loci:
                bed_line = [locus.chr, locus.start, locus.end, target, '', '+']
                target_bed.append(bed_line)
                all_bed.append(bed_line)

        utils.unparse_table(target_bed, target_bed_path, '\t')
    # Now the loci are all stitched up
    utils.unparse_table(all_bed, all_bed_path, '\t')
    return edge_dict


def build_graph(edge_dict, gene_to_enhancer_dict, output_folder, analysis_name, cutoff=1):
    """Build a target graph from the collapsed edge dictionary.

    Require at least n motifs to constitute an edge where n is set by cutoff.
    Default is 1.

    """
    node_list = list(edge_dict.keys())
    node_list.sort()

    # This is only edges between TFs
    graph = nx.DiGraph(name=analysis_name)
    graph.add_nodes_from(node_list)

    # This stores ALL edges identified by motifs
    edge_table = [['SOURCE', 'TARGET', 'CHROM', 'START', 'STOP', 'REGION_ID', 'TF_INTERACTION']]
    edge_output = '{}{}_EDGE_TABLE.txt'.format(output_folder, analysis_name)

    for source in node_list:
        print(source)
        target_list = list(edge_dict[source].keys())
        target_list.sort()
        for target in target_list:

            # Now we need to see which target regions this guy overlaps
            target_regions = gene_to_enhancer_dict[target]
            target_collection = utils.LocusCollection(target_regions, 50)

            # Get the edges hitting that target
            edge_loci = edge_dict[source][target]
            if node_list.count(target) > 0:
                tf_interaction = 1
            else:
                tf_interaction = 0
            # Only add to the graph if this is a TF/TF interaction
            if len(edge_loci) >= cutoff and node_list.count(target) > 0:
                graph.add_edge(source, target)

            # Now for each edge, add to the table
            for edge_locus in edge_loci:
                region_string = ','.join(
                    [locus.id for locus in target_collection.get_overlap(edge_locus)]
                )
                edge_line = [
                    source,
                    target,
                    edge_locus.chr,
                    edge_locus.start,
                    edge_locus.end,
                    region_string,
                    tf_interaction,
                ]
                edge_table.append(edge_line)

    utils.unparse_table(edge_table, edge_output, '\t')
    return graph


def get_clique_ranking(clique_list, out_degree_dict):
        """Clique score generator."""
        for clique in clique_list:
            score = 0
            for gene in clique:
                score += out_degree_dict[gene]
            score = float(score) / len(clique)
            if score > 0 and len(clique) > 2:
                yield (clique, score)


def pairs(self_loops, graph):
    """Recover bidirectional edges."""
    for n in self_loops:
        for m in self_loops:
            if n != m:
                if graph.has_edge(n, m) and graph.has_edge(m, n):
                    yield [n, m]


def format_network_output(graph, output_folder, analysis_name):
    """Takes the networkx graph and returns all figures, tables, etc."""

    # Output the network as a .ntx dictionary of lists
    network_filename = output_folder + analysis_name + '.ntx'
    with open(network_filename, 'wb') as network_file:
        network_dict_of_lists = nx.to_dict_of_lists(graph)
        pickle.dump(network_dict_of_lists, network_file)

    # Output the adjacency list and nodelist
    node_file = output_folder + analysis_name + '_NODELIST.txt'
    if nx.__version__[0] == '1':
        node_list = [[n] for n in graph.nodes_iter()]
    elif nx.__version__[0] == '2':
        node_list = [[n] for n in graph.nodes()]
    else:
        print('ERROR: UNSUPPORTED VERSION OF NETWORKX MODULE')
        sys.exit()
    utils.unparse_table(node_list, node_file, '\t')

    adj_file = output_folder + analysis_name + '_ADJ_LIST.txt'

    if nx.__version__[0] == '1':
        adj_list = graph.adjacency_list()
    elif nx.__version__[0] == '2':
        adj_list = [list(n[1].keys()) for n in graph.adjacency()]
    else:
        print('ERROR: UNSUPPORTED VERSION OF NETWORKX MODULE')
        sys.exit()

    utils.unparse_table(adj_list, adj_file, '\t')

    edges_table = [['From', 'To']]
    for i, gene in enumerate(node_list):
        for j in adj_list[i]:
            newline = [gene[0], j]
            edges_table.append(newline)

    edge_file = output_folder + analysis_name + '_EDGE_LIST.txt'
    utils.unparse_table(edges_table, edge_file, '\t')

    # Make the degree table
    deg_table = [['Tf', 'In_Degree', 'Out_Degree', 'Total_Connections']]
    deg_file = output_folder + analysis_name + '_DEGREE_TABLE.txt'

    # Shouldn't we output the table for the TFs that have motifs only?
    # For canidateMotifs in graph.nodes()....
    for node in graph.nodes():
        newline = [node, graph.in_degree()[node], graph.out_degree()[node], graph.degree()[node]]
        deg_table.append(newline)

    utils.unparse_table(deg_table, deg_file, '\t')

    print('DEFINING THE CORE REGULATORY CIRCUIT')

    autoreg = graph.selfloop_edges()
    self_loops = [x for x, y in autoreg]
    self_loop_file = output_folder + analysis_name + '_SELF_LOOPS.txt'
    utils.unparse_table(self_loops, self_loop_file, '')

    un_dir_graph = nx.from_edgelist(pairs(self_loops, graph))
    clique_gen = find_cliques_recursive(un_dir_graph)
    out_degree_dict = graph.out_degree()

    clique_ranking = get_clique_ranking(clique_gen, out_degree_dict)

    factor_enrichment_dict = {}
    for factor in self_loops:
        factor_enrichment_dict[factor] = 0

    clique_len = 0
    top_cliques = []
    min_clique = ()
    for clique, score in clique_ranking:
        clique_len += 1
        for factor in clique:
            factor_enrichment_dict[factor] += 1

        # Get top 100 cliques
        if clique_len <= 100:
            top_cliques.append((clique, score))
            continue

        if not min_clique:
            min_clique = min(top_cliques, key=lambda x: x[1])

        if score > min_clique[1]:
            top_cliques.remove(min_clique)
            top_cliques.append((clique, score))
            min_clique = min(top_cliques, key=lambda x: x[1])

    top_cliques.sort(reverse=True, key=lambda x: x[1])
    clique_file = output_folder + analysis_name + '_CLIQUE_SCORES_DEGREE.txt'
    utils.unparse_table(top_cliques, clique_file, '\t')

    factor_ranking_table = []
    for factor in self_loops:
        newline = [factor, factor_enrichment_dict[factor] / float(clique_len)]
        factor_ranking_table.append(newline)

    factor_ranking_file = output_folder + analysis_name + '_ENRICHED_CLIQUE_FACTORS.txt'
    utils.unparse_table(factor_ranking_table, factor_ranking_file, '\t')
