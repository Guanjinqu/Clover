import unittest

from clover import align

class TestGlobalAlign(unittest.TestCase):

    def setUp(self) -> None:
        self.example_read_1 = "AAA"
        self.example_read_2 = "ATA"


    def test_global_align(self):
        self.assertEqual([(1,"T")],align.global_align(self.example_read_1,self.example_read_2))


