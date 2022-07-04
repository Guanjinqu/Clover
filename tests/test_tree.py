import unittest

from clover import tree

class TestInsert(unittest.TestCase):
    
    def setUp(self) -> None:
        self.tree = tree.Trie()
        self.word = "AA"
        self.label = 1

    def test_insert(self):
        self.tree.insert(self.word,self.label)
        self.assertEqual(self.tree.children[0].children[0].isEnd,1)

class TestSearchPrefix(unittest.TestCase):
    
    def setUp(self) -> None:
        self.tree = tree.Trie()
        self.word = "AA"
        self.label = 1
        self.tree.insert(self.word,self.label)
        self.word_2 = "AA"

    def test_searchPrefix(self):
        search_result = self.tree.searchPrefix(self.word_2)
        self.assertEqual(search_result,1)

class TestDelete(unittest.TestCase):
    
    def setUp(self) -> None:
        self.tree = tree.Trie()
        self.word = "AA"
        self.label = 1
        self.tree.insert(self.word,self.label)
        self.word_2 = "AA"

    def test_delete(self):
        self.tree.delete(self.word)
        search_result = self.tree.searchPrefix(self.word_2)
        self.assertEqual(search_result,False)

class TestFuzzAlign(unittest.TestCase):
    
    def setUp(self) -> None:
        self.tree = tree.Trie()
        self.word = "AA"
        self.label = 1
        self.tree.insert(self.word,self.label)
        self.word_2 = "AT"

    def test_fuzz_align(self):
        search_result = self.tree.fuzz_align(self.word_2)
        self.assertEqual(search_result,[1,[0]])

class TestFuzzFin(unittest.TestCase):
    
    def setUp(self) -> None:
        self.tree = tree.Trie()
        self.word = "AAAA"
        self.label = 2
        self.tree.insert(self.word,self.label)
        self.word_2 = "ATAA"

    def test_fuzz_fin(self):
        search_result = self.tree.fuzz_fin(self.word_2,10)
        self.assertEqual(search_result,[2,1])