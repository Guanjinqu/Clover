"""Tree Structure Module

This module provides retrieval algorithms based on tree structures, 
including fuzzy search algorithms that allow for horizontal drift.


"""

class Trie:
    """Tree Structure Class

    This class is used to initialize a tree structure without the 
    need to additionally specify the depth of the tree.

    Attributes:
        dna_dict: dict,The elements contained in the sequence.
        node_nums: int,The number of elements contained in the sequence.
        children: int,Number of tree branches.
        isEnd: int,The value to determine if the tree is terminated.
    
    """

    def __init__(self):
        self.dna_dict = {"A":0,"T":1,"G":2,"C":3}
        self.node_nums = len(self.dna_dict)
        self.children = [None] * self.node_nums
        self.isEnd = False

    #Node retrieval function without drift.
    def searchPrefix(self, prefix: str) -> "Trie":
        dict=self.dna_dict
        node = self
        for ch in prefix:
            ch = dict[ch]
            if not node.children[ch]:
                return None
            node = node.children[ch]
        return node.isEnd

    def insert(self, word: str,label:str) -> None:
        """Add a branch to the tree.
        
        Args:
            word: str,The sequence added to the tree.
            label: str,Sequence of labels.

        """
        node = self
        dict=self.dna_dict
        for ch in word:
            ch = dict[ch]
            if not node.children[ch]:
                node.children[ch] = Trie()
            node = node.children[ch]
        node.isEnd = label

    
    def delete(self,word:str):
        """Deletes a branch from the tree.
        
        Args:
            word: str,Sequences deleted from the tree.

        """
        node = self
        dict=self.dna_dict
        for ch in word:
            ch = dict[ch]
            if not node.children[ch]:
                node.children[ch] = Trie()
            node = node.children[ch]
        node.isEnd = False

    def fuzz_align(self,word):
        """Horizontal drift function

        Args:
            word: Sequence of fuzzy retrieval.
        
        return:
            Returns a list with the positions that need to be drifted 
            laterally and the nodes that can be drifted laterally.
        
        """
        node = self 
        dict = self.dna_dict
        num=0
        list=[]
        fin_list=[]
        len_=len(word)
        for ch in word :

            ch = dict[ch]
            if not node.children[ch]:
                
                for i in range(self.node_nums):
                    if not node.children[i]:
                        pass
                    else:
                        list.append(i)
                #print(word,ch,list,num)
                if num +2 < len_ :
                    for k in list :
                        if not node.children[k].children[dict[word[num+1]]]:
                            pass
                        elif not node.children[k].children[dict[word[num+1]]].children[dict[word[num+2]]]:
                            pass
                        else:
                            fin_list.append(k)
                    return [num,fin_list]
                if num +2 == len_:
                    for k in list :
                        if not node.children[k].children[dict[word[num+1]]]:
                            pass
                        else:
                            fin_list.append(k)
                    return [num,fin_list]
                if num +1 == len_:
                    for k in list :
                        fin_list.append(k)
                        return [num,fin_list]
            else:
                node = node.children[ch]
            num = num + 1

        return node.isEnd

    def fuzz_fin(self,word,max_value):
        """Fuzzy search with horizontal drift.

        Args:
            word: str,Sequence of search
            max_value: int,The maximum number of horizontal drifts.

        return:
            Returns a list, the first element of which is the index of the final matched 
            core sequence, and the second element is the number of horizontal drifts.
        """

        list=[[word,0]]
        
        fin_list=["",1000]
        dict2  = {}
        for key in self.dna_dict:
            dict2[self.dna_dict[key]] = key
        error_list=[]
        while True :
            if list == [] or fin_list[1] == 0 :

                break
            dna = list[0]
            
            if dna[1] > max_value :
                break
            del list[0]
            a = self.fuzz_align(dna[0])
            error_list.append([dna,a])
            if type(a) == int :
                if dna[1] < fin_list[1] :
                    fin_list=[a,dna[1]]
            elif self.fuzz_align(word)[1] == [] :
                break
            elif a[0] == len(dna[0])-1:
                for i in range(len(a[1])):
                    i=a[1][i]
                    k = dna[0][:a[0]]+dict2[i]
                    list.append([k,dna[1]+1])                
            else:
                for i in range(len(a[1])):
                    i=a[1][i]
                    k = dna[0][:a[0]]+dict2[i]+dna[0][a[0]-len(dna[0])+1:]
                    list.append([k,dna[1]+1])
        fin_list.append(error_list)
        return fin_list[:2]
    

