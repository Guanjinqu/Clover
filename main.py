import os
from multiprocessing import Process, Queue
import time
from tqdm import tqdm
import align as ag
import load_config as lc
import tree as tr
from collections import  Counter
class MyProcess(Process):
    def __init__(self, name, data, q_output):
        Process.__init__(self)
        self.config_dict=lc.out_put_config()
        self.name = name
        self.data = data
        self.q_output = q_output
        self.q_output_temp = []
        self.a_tree = tr.Trie() 
        self.b_tree = tr.Trie()  
        self.c_tree = tr.Trie()
        if self.config_dict['other_tree_nums'] == 2 :
            self.d_tree = tr.Trie()

        self.Cluster_size_threshold = self.config_dict['Cluster_size_threshold']
        self.ref_list={}    #ref的集合
        self.ref_dict={}     #ref的字典
        self.ref_error_dict={}
        self.num_dict={}
        self.tag_dict={}
        self.index_list=[]   #该列表输出聚类的结果，即()
        self.now_clust_threshold = self.config_dict['now_clust_threshold']
        self.read_len= self.config_dict['read_len']   #序列的长度
        self.dna_tree_nums = self.config_dict['end_tree_len'] #构造前后缀的长度
        self.fuzz_list = [self.config_dict['thd_tree_loc'],self.config_dict['four_tree_loc'],self.config_dict['other_tree_len']] #构造模糊查看的集合
        self.loc_nums = self.config_dict['Vertical_drift'] #允许前后位移的程度
        self.tag_nums = self.config_dict['tag_nums']      #训练集中标签数量
        self.align_swicth= self.config_dict['align_fuc'] #是否开启比对功能
        self.fuzz_tree_nums= self.config_dict['Horizontal_drift']   #运行树模糊搜索的程度
        self.h_index=self.config_dict['h_index_nums']
        self.e_index=self.config_dict['e_index_nums']
        self.test_num = 0



    def cluster(self,read):

        self.test_num=self.test_num+1
        line_=read.split()

        dna_num=self.test_num
        if self.config_dict['Virtual_mode'] == False or self.config_dict['mmr_mode'] == True:
            dna_index=line_[0]
            dna_str = line_[1]
            dna_tag=line_[0]
        else:
            dna_str = line_[1]
            dna_tag = line_[0]
        if self.h_index == 0 :
            dna_a_str = dna_str[:self.dna_tree_nums]
        else:
            dna_a_str = dna_str[self.h_index:self.h_index+self.dna_tree_nums]
        if self.e_index == 0 : 
            dna_b_str = dna_str[-self.dna_tree_nums:]
        else:
            dna_b_str = dna_str[-self.e_index-self.dna_tree_nums:-self.e_index]
        dna_str_num=len(dna_str)
        if dna_str_num < self.config_dict['read_len_min'] or "N" in dna_str:
            pass 
        else:
            a_align = self.a_tree.fuzz_fin(dna_a_str,self.config_dict['tree_threshold']) 
            
            if a_align[1]<self.fuzz_tree_nums:
                if self.config_dict['Virtual_mode'] == False:
                    
                    self.index_list.append((dna_index,self.ref_dict[a_align[0]][0]))
                else:
                    self.ref_dict[a_align[0]].append(dna_tag)
                
                if self.align_swicth is True:
                    error_list = []
                    align_list=a_align
                    if self.config_dict["now_align_alg"] == True:
                        error_list=ag.align(self.ref_list[align_list[0]],dna_str)
                    else:
                        if dna_str_num==self.read_len :
                            if self.ref_list[align_list[0]] == dna_str:
                                pass
                            else:
                                error_list=ag.align(self.ref_list[align_list[0]],dna_str)
                        
                    for line in error_list:

                        if align_list[0] in self.ref_error_dict : 
                            self.ref_error_dict[align_list[0]]["nums"]=self.ref_error_dict[align_list[0]]["nums"] + 1
                            if line[0] in  self.ref_error_dict[align_list[0]]:
                                self.ref_error_dict[align_list[0]][line[0]]=self.ref_error_dict[align_list[0]][line[0]]+1
                            else:
                                self.ref_error_dict[align_list[0]][line[0]]=1
                            if self.ref_error_dict[align_list[0]][line[0]]/self.ref_error_dict[align_list[0]]["nums"] >0.5 and self.ref_error_dict[a_align[0]][line[0]]>5 :
                                now_read=self.ref_list[align_list[0]][:line[0]]+line[1]+self.ref_list[align_list[0]][line[0]+1:]

                                self.a_tree.insert(now_read[:self.dna_tree_nums],align_list[0])
                                self.ref_list[align_list[0]]=now_read

                                self.ref_error_dict[align_list[0]]={}
                                self.ref_error_dict[align_list[0]]["nums"]=1
                                self.ref_error_dict[align_list[0]][line[0]]=1                                   
                        else:
                            self.ref_error_dict[align_list[0]]={}
                            self.ref_error_dict[align_list[0]]["nums"]=1
                            self.ref_error_dict[align_list[0]][line[0]]=1
                        
            elif self.b_tree.fuzz_fin(dna_b_str,self.config_dict['tree_threshold'])[1] < self.fuzz_tree_nums:
                b_align=self.b_tree.fuzz_fin(dna_b_str,self.config_dict['tree_threshold'])
                
                if self.config_dict['Virtual_mode'] == False:
                    self.index_list.append((dna_index,self.ref_dict[b_align[0]][0]))
                else:
                    self.ref_dict[b_align[0]].append(dna_tag)

                if self.align_swicth is True:
                    error_list = []
                    align_list=b_align
                    if self.config_dict["now_align_alg"] == True:
                        error_list=ag.align(self.ref_list[align_list[0]],dna_str)
                    else:
                        if dna_str_num==self.read_len :
                            if self.ref_list[align_list[0]] == dna_str:
                                pass
                            else:
                                error_list=ag.align(self.ref_list[align_list[0]],dna_str)
                        
                    for line in error_list:
                        if align_list[0] in self.ref_error_dict : 
                            self.ref_error_dict[align_list[0]]["nums"]=self.ref_error_dict[align_list[0]]["nums"] + 1
                            if line[0] in  self.ref_error_dict[align_list[0]]:
                                self.ref_error_dict[align_list[0]][line[0]]=self.ref_error_dict[align_list[0]][line[0]]+1
                            else:
                                self.ref_error_dict[align_list[0]][line[0]]=1
                            if self.ref_error_dict[align_list[0]][line[0]]/self.ref_error_dict[align_list[0]]["nums"] >0.5 and self.ref_error_dict[a_align[0]][line[0]]>5 :
                                now_read=self.ref_list[align_list[0]][:line[0]]+line[1]+self.ref_list[align_list[0]][line[0]+1:]

                                self.a_tree.insert(now_read[:self.dna_tree_nums],align_list[0])
                                self.ref_list[align_list[0]]=now_read

                                self.ref_error_dict[align_list[0]]={}
                                self.ref_error_dict[align_list[0]]["nums"]=1
                                self.ref_error_dict[align_list[0]][line[0]]=1                                   
                        else:
                            self.ref_error_dict[align_list[0]]={}
                            self.ref_error_dict[align_list[0]]["nums"]=1
                            self.ref_error_dict[align_list[0]][line[0]]=1
                
            else:
                if dna_str_num >= self.config_dict['read_len_min'] :
                    fin_align=["",1000]
                    for i in self.loc_nums:
                        if self.h_index == 0 :
                            dna_c_str=dna_str[self.fuzz_list[0]-i:self.fuzz_list[0]+self.fuzz_list[2]-i]
                        else:
                            dna_c_str=dna_str[self.h_index+self.fuzz_list[0]-i:self.h_index+self.fuzz_list[0]+self.fuzz_list[2]-i]
                        c_align = self.c_tree.fuzz_fin(dna_c_str,self.config_dict['tree_threshold']) 
                        if c_align[1]<fin_align[1] :
                                fin_align=c_align
                        if self.config_dict['other_tree_nums'] == 2 :
                            if self.e_index == 0 :
                                dna_d_str=dna_str[self.read_len-2-self.fuzz_list[1]-i:self.read_len-2-self.fuzz_list[1]+self.fuzz_list[2]-i]
                            else:
                                dna_d_str=dna_str[self.read_len-2-self.fuzz_list[1]-i-self.e_index:self.read_len-2-self.fuzz_list[1]+self.fuzz_list[2]-i-self.e_index]
                            d_align = self.d_tree.fuzz_fin(dna_d_str,self.config_dict['tree_threshold'])
                            if d_align[1]<fin_align[1] :
                                fin_align=d_align

                    if fin_align[1] < self.fuzz_tree_nums :
                        if self.config_dict['Virtual_mode'] == False:
                            self.index_list.append((dna_index,self.ref_dict[fin_align[0]][0]))
                        else:
                            self.ref_dict[fin_align[0]].append(dna_tag)
                        
                        if self.align_swicth is True:
                            error_list = []
                            align_list=fin_align
                            if self.config_dict["now_align_alg"] == True:
                                error_list=ag.align(self.ref_list[align_list[0]],dna_str)
                            else:
                                if dna_str_num==self.read_len :
                                    if self.ref_list[align_list[0]] == dna_str:
                                        pass
                                    else:
                                        error_list=ag.align(self.ref_list[align_list[0]],dna_str)
                                
                            for line in error_list:
                                if align_list[0] in self.ref_error_dict : 
                                    self.ref_error_dict[align_list[0]]["nums"]=self.ref_error_dict[align_list[0]]["nums"] + 1
                                    if line[0] in  self.ref_error_dict[align_list[0]]:
                                        self.ref_error_dict[align_list[0]][line[0]]=self.ref_error_dict[align_list[0]][line[0]]+1
                                    else:
                                        self.ref_error_dict[align_list[0]][line[0]]=1
                                    if self.ref_error_dict[align_list[0]][line[0]]/self.ref_error_dict[align_list[0]]["nums"] >0.5 and self.ref_error_dict[a_align[0]][line[0]]>5 :
                                        now_read=self.ref_list[align_list[0]][:line[0]]+line[1]+self.ref_list[align_list[0]][line[0]+1:]

                                        self.a_tree.insert(now_read[:self.dna_tree_nums],align_list[0])
                                        self.ref_list[align_list[0]]=now_read

                                        self.ref_error_dict[align_list[0]]={}
                                        self.ref_error_dict[align_list[0]]["nums"]=1
                                        self.ref_error_dict[align_list[0]][line[0]]=1                                   
                                else:
                                    self.ref_error_dict[align_list[0]]={}
                                    self.ref_error_dict[align_list[0]]["nums"]=1
                                    self.ref_error_dict[align_list[0]][line[0]]=1



                    """#if fin_align[1] >= 2 and fin_align[1] <=5 :
                        if Levenshtein.distance(dna_str,self.ref_list[fin_align[0]]) <30 :
                            self.ref_dict[fin_align[0]].append(dna_tag)
                        else:
                            #self.ref_list[dna_num]=dna_str
                            self.ref_dict[dna_num]=[dna_tag]
                            self.a_tree.insert(dna_a_str,dna_num)
                            self.b_tree.insert(dna_b_str,dna_num)
                            self.c_tree.insert(dna_str[self.fuzz_list[0]-i:self.fuzz_list[0]+self.fuzz_list[2]-i],dna_num)
                            self.d_tree.insert(dna_str[self.read_len-2-self.fuzz_list[1]-i:self.read_len-2-self.fuzz_list[1]+self.fuzz_list[2]-i],dna_num)"""                            
                    if fin_align[1] >= self.now_clust_threshold :
                        if self.config_dict['align_fuc'] == True:
                            self.ref_list[dna_num]=dna_str
                        self.ref_dict[dna_num]=[dna_tag]
                        self.a_tree.insert(dna_a_str,dna_num)
                        self.b_tree.insert(dna_b_str,dna_num)
                        self.c_tree.insert(dna_str[self.fuzz_list[0]-i:self.fuzz_list[0]+self.fuzz_list[2]-i],dna_num)
                        self.d_tree.insert(dna_str[self.read_len-2-self.fuzz_list[1]-i:self.read_len-2-self.fuzz_list[1]+self.fuzz_list[2]-i],dna_num)


    def run(self):

        self.num_dict[self.name+"sum_read_num"]=0
        self.num_dict[self.name+"error_num"]=0
        self.num_dict[self.name+"sum_cluster_num"]=0
        self.num_dict[self.name+"sum_tag"]=[]
        len_=len(self.name) 
        
        st1=time.time()
        if self.config_dict['fast_mode'] == True :
            
            for line in tqdm(self.data):
                self.cluster(line)
            if 'output_file' in self.config_dict :

                self.num_dict[self.name+"index_list"]=self.index_list
            if self.config_dict['Virtual_mode'] == True :
                tag_sum=0
                tag_error=0
                nums_=0
                nums_sum=0
                if self.config_dict['Virtual_mode'] == True:
                    for key in self.ref_dict:
                        tag_len=len(self.ref_dict[key])
                        
                        if tag_len > self.Cluster_size_threshold :
                            nums_sum=nums_sum+1 
                            tag_c=Counter(self.ref_dict[key]).most_common(1)[0][1]
                            tag_error=tag_error+tag_len-tag_c
                            tag_sum=tag_sum+tag_len   
                            tag=self.ref_dict[key][0]
                            if tag in self.tag_dict:
                                pass
                            else:
                                nums_=nums_+1
                                self.tag_dict[tag]=1
                        self.ref_dict[key]=[]
                    
                    self.num_dict[self.name+"sum_read_num"]+=tag_sum
                    self.num_dict[self.name+"error_num"]+=tag_error
                    self.num_dict[self.name+"sum_cluster_num"]=nums_sum
                    for key in self.tag_dict:
                        self.num_dict[self.name+"sum_tag"].append(key)
        elif self.config_dict['mmr_mode'] == True:

            output_nums = 0
            f=open(self.config_dict['input_path'],"r")
            
            if 'output_file' in self.config_dict :
                output_file = open(self.name+"_"+self.config_dict['output_file'],'a')
                output_file.write("[")
            while True :
                output_nums += 0
                line = f.readline()
                if line == [] or line == "" :
                    break
                if self.config_dict['processes_nums'] == 0:
                    self.cluster(line)
                elif line.split()[1][:len_] == self.name :
                    self.cluster(line)
                    if 'output_file' in self.config_dict and output_nums%100000 == 0:
                        for i in self.index_list:
                            output_file.write(str(i))
                            output_file.write(',')
                            self.index_list=[]
            if 'output_file' in self.config_dict :
                for i in self.index_list:
                    output_file.write(str(i))
                    output_file.write(',')
                output_file.close()
                with open(self.name+"_"+self.config_dict['output_file'], 'rb+') as output_file:
                    output_file.seek(-1, os.SEEK_END)
                    output_file.truncate()
                output_file = open(self.name+"_"+self.config_dict['output_file'],'a')
                output_file.write("]")    
                output_file.close()
        else :
            if self.config_dict['Virtual_mode'] == True :
                f=open(self.config_dict['input_path'],"r")
                while True :
                    line = f.readline()
                    if line == [] or line == "" :
                        break
                    if self.config_dict['processes_nums'] == 0:
                        self.cluster(line)
                    elif line.split()[1][:len_] == self.name :
                        self.cluster(line)

                if 'output_file' in self.config_dict :
                    self.num_dict[self.name+"index_list"]=self.index_list
                tag_sum=0
                tag_error=0
                nums_=0
                nums_sum=0
                for key in self.ref_dict:
                    tag_len=len(self.ref_dict[key])
                    
                    if tag_len > self.Cluster_size_threshold :
                        nums_sum=nums_sum+1 
                        tag_c=Counter(self.ref_dict[key]).most_common(1)[0][1]
                        tag_error=tag_error+tag_len-tag_c
                        tag_sum=tag_sum+tag_len   
                        tag=self.ref_dict[key][0]
                        if tag in self.tag_dict:
                            pass
                        else:
                            nums_=nums_+1
                            self.tag_dict[tag]=1
                    self.ref_dict[key]=[]
                
                self.num_dict[self.name+"sum_read_num"]+=tag_sum
                self.num_dict[self.name+"error_num"]+=tag_error
                self.num_dict[self.name+"sum_cluster_num"]=nums_sum
                for key in self.tag_dict:
                    self.num_dict[self.name+"sum_tag"].append(key)
            else:
                f=open(self.config_dict['input_path'],"r")
                while True :
                    line = f.readline()
                    if line == [] or line == "" :
                        break
                    if self.config_dict['processes_nums'] == 0:
                        self.cluster(line)
                    elif line.split()[1][:len_] == self.name :
                        self.cluster(line)
                    if 'output_file' in self.config_dict :
                        self.num_dict[self.name+"index_list"]=self.index_list
        self.q_output.put(self.num_dict)
def all_permutations(items, length):
    res = []
    def track_back(tmp_permutation):
        if len(tmp_permutation) == length:
            res.append(tmp_permutation[:])
        else:
            for i in items:
                tmp_permutation.append(i)
                track_back(tmp_permutation)
                tmp_permutation.pop()
    track_back([])
    return ["".join(i) for i in res]


if __name__ == '__main__':
    # 按照序列前n位创建进程
    #******************************************************************************************
    config_dict=lc.out_put_config()

    tag_nums=config_dict['tag_nums']
    PROCESS_INDEX= int(config_dict['processes_nums'])
    
    #******************************************************************************************
    
    print(config_dict['tag'])
    if PROCESS_INDEX == 0 and config_dict['fast_mode'] == True :
        N_PROCESS=1
        process_names=['all']
        print("Generating data")
        data_dict = {}
        data_dict['all']=[]
        f= open(config_dict['input_path'],'r').readlines()
        for i in tqdm(f):
            if "N" in i :
                pass
            else:
                data_dict['all'].append(i)
    elif PROCESS_INDEX == 0 and config_dict['fast_mode'] == False :
        N_PROCESS=1
        process_names=['all']
        data_dict = {}
        data_dict['all']=[]
    elif PROCESS_INDEX > 0 :
        N_PROCESS = 4**PROCESS_INDEX
        process_names = all_permutations(["A","T","G","C"],PROCESS_INDEX)


        
        st_pre = time.time()
        data_dict = {}.fromkeys(process_names)
        for i in data_dict:
            data_dict[i] = []


        if config_dict['fast_mode'] == True:
            print("Pre-process the data")
            f= open(config_dict['input_path'],'r').readlines()
            for i in tqdm(f):
                if "N" in i or "*" in i :
                    pass
                else:
                    data_dict[i.split()[1][:PROCESS_INDEX]].append(i)
        else:
            data_dict[i] = []




    q_output = Queue(N_PROCESS*2)

    st = time.time()
    process_dict = {}.fromkeys(process_names)

    for i in process_dict:
        process_dict[i] = MyProcess(i, data_dict[i], q_output)

    for i in process_dict:
        process_dict[i].start()

    st3=time.time()
    count_dict={}
    i=0
    while True :
        if i == N_PROCESS :
            break
        count_dict.update(q_output.get())
        i=i+1
    print("Time:",time.time()-st3)

    for i in process_dict:
        process_dict[i].join()    
    
    new_count_dict={}
    new_count_dict["sum_read_num"]=0
    new_count_dict["error_num"]=0
    new_count_dict["sum_cluster_num"]=0
    new_count_dict["sum_tag"]={}
    new_count_dict["index_list"]=[]
    for key in process_dict :
        new_count_dict["sum_read_num"]+=count_dict[key+"sum_read_num"]
        new_count_dict["error_num"]+=count_dict[key+"error_num"]
        new_count_dict["sum_cluster_num"]+=count_dict[key+"sum_cluster_num"]
        tag_list=count_dict[key+"sum_tag"]
        
        for i in tag_list:
            if i in new_count_dict["sum_tag"]:
                pass
            else:
                new_count_dict["sum_tag"][i]=1
    if 'output_file' in config_dict and config_dict["mmr_mode"] is not True:
        for key in process_dict :
            new_count_dict["index_list"]+=count_dict[key+"index_list"]
        output_file = open(config_dict['output_file'],'w')
        output_file.write(str(new_count_dict["index_list"]))
        output_file.close()
    if config_dict['Virtual_mode'] == True :
        #print("Number of reads processed:",new_count_dict["sum_read_num"],new_count_dict["error_num"],new_count_dict["sum_cluster_num"])
        print("Number of reads processed:",new_count_dict["sum_read_num"])
        print("Accuracy：",(new_count_dict["sum_read_num"]-new_count_dict["error_num"])/new_count_dict["sum_read_num"])
        if config_dict['tag_mode'] == True :
            print("Coverage：",len(new_count_dict["sum_tag"])/int(tag_nums))
            print("Redundancy Rate：",(new_count_dict["sum_cluster_num"]-int(tag_nums))/int(tag_nums))
    elif config_dict['Statistical_model'] == True :
        pass

