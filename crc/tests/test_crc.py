"""Test the `crc` main function."""

from crc.bin.crc3 import crc

import os
import pytest  # noqa: F401
import sys

TEST_FILES_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), 'files'))
TEST_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'tests'))


def test_crc():
    """Test crc functionality."""
    crc(
        '{}/test_enhancers.bed'.format(TEST_FILES_DIR),
        'HG19',
        '{}/Chromosomes/'.format(TEST_FILES_DIR),
        TEST_DIR,
        'test',
        subpeak_file='{}/mock_regions.bed'.format(TEST_FILES_DIR),
    )

    scores = []
    with open(os.path.join(TEST_DIR, 'test_CLIQUE_SCORES_DEGREE.txt')) as infile:
        for line in infile:
            scores.append(float(line.split('\t')[1].strip('\n')))

    if (sys.version_info > (3, 0)):
        test_scores = [8.25, 8.0, 7.75, 7.333333333333333]
    else:
        test_scores = [8.25, 8.0, 7.75, 7.33333333333]

    assert scores == test_scores, 'Clique scores do not match!'
