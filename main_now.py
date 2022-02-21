import align as ag
import time
import tree as tr
import tqdm
from collections import  Counter
import os
import psutil

test_num=0
ref_list={}    #ref的集合
ref_dict={}     #ref的字典
ref_error_dict={}
ref_minhash_dict={}
minhash_list=["*"]
a_tree = tr.Trie() #构造前向树
b_tree = tr.Trie()  #构造后向树
c_tree = tr.Trie()
d_tree = tr.Trie()


a_hash_list = [10,20,30,40,50,60,70,71,72,73,74,75,76] #构造前项哈希取值
b_hash_list = [71,72,73,74,75,76] #构造后项哈希取值

loc_nums = [-1,-2,0,1,2] #允许前后位移的程度
hash_error_nums = 3 #允许哈希特征冗余的值

align_swicth= False #是否开启比对功能
fuzz_tree_nums= 3   #运行树模糊搜索的程度
memory_num=0
test_num2=0
time2=0

st=time.time()

def get_current_memory_gb() -> int:
# 获取当前进程内存占用。
    pid = os.getpid()
    p = psutil.Process(pid)
    info = p.memory_full_info()
    return info.uss / 1024. / 1024. / 1024.

memory_num=0

memory_list=[]
st_p=time.time()
#**********************************************************************

f = open('sort_DNA_1000w.txt','r').readlines()  #文件名称
len_dna=152  #117 152  #dna长度
tag_nums = 49119      #训练集中标签数量 51、 514 、4948、49119
#g = open("memory_goldman6.txt","w")
del_loc=[0,0]  #前后端删除的碱基数量#注意！我把86行删掉了 使用的时候别忘加上！！！
dna_tree_nums =  15#构造前后缀的长度
fuzz_list = [40,40,15] #构造模糊查看的集合
#**********************************************************************
start_time = time.time()

print("预处理时间：",time.time()-st_p)
for i in tqdm.tqdm(f):
    #line = f.readline()
    #if line == "" :
    #    break

    test_num=test_num+1
    #if test_num == 50000000:
    #    break
    line_=i.split()
    memory_num+=1
    if  line_== [] :
        continue
    
    dna_num=test_num
    dna_str = line_[1]
    dna_tag = line_[0]
    dna_str_num=len(dna_str)
    if 'N' in dna_str or dna_str_num<len_dna-2  or dna_tag == "*":
        pass
    else:

        dna_a_str = dna_str[del_loc[0]:dna_tree_nums+del_loc[0]]
        dna_b_str = dna_str[-dna_tree_nums:]
        
        a_align = a_tree.fuzz_fin(dna_a_str) 
        
        if a_align[1]<fuzz_tree_nums:

            ref_dict[a_align[0]].append(dna_tag)
            
            if align_swicth is True:
                error_list = []
                if dna_str_num==len_dna :
                    if ref_list[a_align[0]] == dna_str:
                        pass
                    else:
                        error_list=ag.align(ref_list[a_align[0]],dna_str)
                    
                    for line in error_list:
                        #print(line)
                        if a_align[0] in ref_error_dict : 
                            ref_error_dict[a_align[0]]["nums"]=ref_error_dict[a_align[0]]["nums"] + 1
                            if line[0] in  ref_error_dict[a_align[0]]:
                                ref_error_dict[a_align[0]][line[0]]=ref_error_dict[a_align[0]][line[0]]+1
                            else:
                                ref_error_dict[a_align[0]][line[0]]=1
                            if ref_error_dict[a_align[0]][line[0]]/ref_error_dict[a_align[0]]["nums"] >0.5 and ref_error_dict[a_align[0]][line[0]]>5 :
                                now_read=ref_list[a_align[0]][:line[0]]+line[1]+ref_list[a_align[0]][line[0]+1:]
                                #print(ref_error_dict[a_align[0]][line[0]])
                                #print(ref_list[a_align[0]],dna_a_str)
                                #a_tree.delete(ref_list[a_align[0]][:dna_tree_nums])
                                a_tree.insert(now_read[:dna_tree_nums],a_align[0])
                                ref_list[a_align[0]]=now_read
                                #print(a_align[0],dna_tag)
                                ref_error_dict[a_align[0]]={}
                                ref_error_dict[a_align[0]]["nums"]=1
                                ref_error_dict[a_align[0]][line[0]]=1                                   
                        else:
                            ref_error_dict[a_align[0]]={}
                            ref_error_dict[a_align[0]]["nums"]=1
                            ref_error_dict[a_align[0]][line[0]]=1
                    
        elif b_tree.fuzz_fin(dna_b_str)[1] < fuzz_tree_nums :
            b_align=b_tree.fuzz_fin(dna_b_str)
            
            ref_dict[b_align[0]].append(dna_tag)

            if align_swicth is True:
                error_list = []
                if dna_str_num==len_dna :
                    if ref_list[b_align[0]] == dna_str:
                        pass
                    else:
                        error_list=ag.align(ref_list[b_align[0]],dna_str)
                    
                    for line in error_list:
                        #print(line)
                        if b_align[0] in ref_error_dict : 
                            ref_error_dict[b_align[0]]["nums"]=ref_error_dict[b_align[0]]["nums"] + 1
                            if line[0] in  ref_error_dict[b_align[0]]:
                                ref_error_dict[b_align[0]][line[0]]=ref_error_dict[b_align[0]][line[0]]+1
                            else:
                                ref_error_dict[b_align[0]][line[0]]=1
                            if ref_error_dict[b_align[0]][line[0]]/ref_error_dict[b_align[0]]["nums"] >0.5 and ref_error_dict[b_align[0]][line[0]]>5 :
                                now_read=ref_list[b_align[0]][:line[0]]+line[1]+ref_list[b_align[0]][line[0]+1:]
                                #print(ref_error_dict[a_align[0]][line[0]])
                                #print(ref_list[a_align[0]],dna_a_str)
                                #a_tree.delete(ref_list[a_align[0]][:dna_tree_nums])
                                b_tree.insert(now_read[-dna_tree_nums:],b_align[0])
                                ref_list[b_align[0]]=now_read
                                ref_error_dict[b_align[0]]={}
                                ref_error_dict[b_align[0]]["nums"]=1
                                ref_error_dict[b_align[0]][line[0]]=1                                
                                #print(a_align[0],dna_tag)
                                
                        else:
                            ref_error_dict[b_align[0]]={}
                            ref_error_dict[b_align[0]]["nums"]=1
                            ref_error_dict[b_align[0]][line[0]]=1

            
        else:
            if dna_str_num == len_dna :
                fin_align=["",1000]
                for i in loc_nums:
                    dna_c_str=dna_str[fuzz_list[0]-i+del_loc[0]:fuzz_list[0]+fuzz_list[2]-i+del_loc[0]]
                    c_align = c_tree.fuzz_fin(dna_c_str) 
                    if c_align[1]< 2 :
                        ref_dict[c_align[0]].append(dna_tag)
                        fin_align=c_align
                    else:    
                        if c_align[1]<fin_align[1] :
                            fin_align=c_align
                    
                        dna_d_str=dna_str[len_dna-2-fuzz_list[1]-i-del_loc[1]:len_dna-2-fuzz_list[1]+fuzz_list[2]-i-del_loc[1]]
                        d_align = d_tree.fuzz_fin(dna_d_str)
                        if d_align[1]< 2 :
                            ref_dict[d_align[0]].append(dna_tag)
                            fin_align=d_align
                        else:
                            if d_align[1]<fin_align[1] :
                                fin_align=d_align
                #if fin_align[1] >= 2 and fin_align[1] <=8 :
                #    if Levenshtein.distance(dna_str,ref_list[fin_align[0]]) <30 :
                #        ref_dict[fin_align[0]].append(dna_tag)
                if fin_align[1] >= 8:

                    ref_list[dna_num]=dna_str
                    ref_dict[dna_num]=[dna_tag]
                    a_tree.insert(dna_a_str,dna_num)
                    b_tree.insert(dna_b_str,dna_num)
                    c_tree.insert(dna_str[fuzz_list[0]-i+del_loc[0]:fuzz_list[0]+fuzz_list[2]-i+del_loc[0]],dna_num)
                    d_tree.insert(dna_str[len_dna-2-fuzz_list[1]-i-del_loc[1]:len_dna-2-fuzz_list[1]+fuzz_list[2]-i-del_loc[1]],dna_num)


#g.write(str(memory_list))        
#g.close()   
#print(get_current_memory_gb())

                            
et=time.time()
tag_sum=0
tag_error=0
nums_=0
nums_sum=0
test_tag=""

#print(ref_dict)
num_tag_dict={}

for key in ref_dict:
    if len(ref_dict[key]) >10 :
        tag_len=len(ref_dict[key])
        tag_c=Counter(ref_dict[key]).most_common(1)[0][1]
        tag_error=tag_error+tag_len-tag_c
        tag_sum=tag_sum+tag_len
    #print(Counter(ref_dict[key]))
        nums_sum=nums_sum+1
        num_tag=ref_dict[key][0]
        if num_tag in num_tag_dict:
            pass
        else:
            nums_=nums_+1
            num_tag_dict[num_tag]=1




#print(time2)
print(nums_sum)
print("耗时：",et-start_time)
print("准确率：",(tag_sum-tag_error)/tag_sum)
print("覆盖率：",nums_/tag_nums)
print("冗余率：",(nums_sum-tag_nums)/tag_nums,len(ref_dict))



