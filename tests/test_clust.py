import unittest

from clover import main

class TestClust(unittest.TestCase):

    def setUp(self) -> None:
        self.example_process = main.MyProcess("test",[],[])
        self.example_process.now_clust_threshold = 0
        self.example_process.read_len = 6
        self.example_process.dna_tree_nums = 4
        self.example_process.fuzz_list = [1,1,4]
        self.example_process.loc_nums = [0]
        self.example_process.fuzz_tree_nums = 1
        self.example_process.config_dict['read_len_min'] = 4


    def test_global_align(self):
        self.example_process.cluster("1 AAAAAA")
        self.example_process.cluster("2 AAAAAA")
        cluster_result = self.example_process.ref_dict[1]
        self.assertEqual(cluster_result,["1","2"])

