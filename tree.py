from re import search
import numpy as np
import pandas as pd



class Trie:
    def __init__(self):
        self.children = [None] * 4
        self.isEnd = False

    
    def searchPrefix(self, prefix: str) -> "Trie":
        dict={"A":0,"T":1,"G":2,"C":3}
        node = self
        for ch in prefix:
            ch = dict[ch]
            if not node.children[ch]:
                return None
            node = node.children[ch]
        return node

    def fuzz_search(self,word,nums):
        dict={"A":0,"T":1,"G":2,"C":3}    
        node = self
        for ch in word:
            ch = dict[ch]
            if not node.children[ch]:
                if nums==0:
                    return None
                else:
                    nums=nums-1
                    for i in range(3):
                        if node.children[i] is not None :
                            ch = i
            if not node.children[ch]:
                return None
            node = node.children[ch]
        return node is not None and node.isEnd



    def insert(self, word: str,label:str) -> None:
        node = self
        dict={"A":0,"T":1,"G":2,"C":3}
        for ch in word:
            ch = dict[ch]
            if not node.children[ch]:
                node.children[ch] = Trie()
            node = node.children[ch]
        node.isEnd = label

    def search(self, word: str) -> bool:
        node = self.searchPrefix(word)
        return node is not None and node.isEnd

    def startsWith(self, prefix: str) -> bool:
        return self.searchPrefix(prefix) is not None

    def delete(self,word:str):
        node = self
        dict={"A":0,"T":1,"G":2,"C":3}
        for ch in word:
            ch = dict[ch]
            if not node.children[ch]:
                node.children[ch] = Trie()
            node = node.children[ch]
        node.isEnd = False

    def fuzz_align(self,word):
        node = self 
        dict = {"A":0,"T":1,"G":2,"C":3}
        num=0
        list=[]
        fin_list=[]
        len_=len(word)
        for ch in word :

            ch = dict[ch]
            if not node.children[ch]:
                
                for i in range(4):
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
        nums=0
        list=[[word,0]]
        
        fin_list=["",1000]
        dict2 = {0:"A",1:"T",2:"G",3:"C"}
        error_list=[]
        while True :
            #print(list)
            if list == [] or fin_list[1] == 0 :

                break
            dna = list[0]
            if dna[1] > max_value :
                break
            #print(dna)
            del list[0]
            a = self.fuzz_align(dna[0])
            #print(a)
            #print(a)
            error_list.append([dna,a])
            if type(a) == int :
                if dna[1] < fin_list[1] :
                    fin_list=[a,dna[1]]
            elif self.fuzz_align(word)[1] == [] :
                #print("yes")
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
        return fin_list
    

