import nose
import numpy as np
import h5py
from Bio import SeqIO
import unittest
import scripts.sequence_aggregation as seq


class SequenceAggregationTest(unittest.TestCase):

    def setUp(self):
        filename = '../../shorter.fna'
        self.original_fasta = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
        seq.sequence_to_array(filename)
        self.result = h5py.File('../scripts/aggregated_dict.hdf5', 'r')

    def test_length(self):
        """ The shorter.fna file has 3 keys, so the result h5py file should have 3 keys. """
        assert len(self.result.keys()) == len(self.original_fasta.keys()) == 5

    # order: [A, T, G, C, N, Other]

    def test_that_first_letter_is_g(self):
        """ The first letter in the first dataset, which is G, should be [0, 0, 1, 0, 0, 0] """
        assert np.array_equal(self.result.get('AE014298.5')[0], [0, 0, 1, 0, 0, 0])

    def test_that_last_letter_is_t(self):
        """ The last letter in the first dataset, which is t, should be [0, 1, 0, 0, 0, 0] """
        assert np.array_equal(self.result.get('AE014298.5')[len(self.result.keys()) - 1], [0, 1, 0, 0, 0, 0])
