import  difflib
import time 
"""def align(a,b):
    list=[]
    num_=0
    diff_content=difflib.Differ().compare(a,b)
    a=''.join(diff_content)
    #print(a)
    num_2=0
    while True:
        b=a.find("+",num_)
        c=a.find("-",num_)
        maxbc=max(b,c)
        minbc=min(b,c)
        if b==-1 and c == -1 :
            #print(list)
            break
        elif minbc+3 ==maxbc:
            error_lco=minbc/3-num_2
            error_str = a[maxbc+2]
            error_type = "*"
            list.append([error_lco,error_type,error_str])
            num_2=num_2+1
            num_=maxbc+1
        elif minbc==-1:
            error_lco=maxbc/3-num_2
            error_str = a[maxbc+2]
            error_type = a[maxbc]
            if  error_type == "+":
                num_2=num_2+1
            list.append([error_lco,error_type,error_str])
            num_=maxbc+1
        else:
            error_lco=minbc/3-num_2
            error_str = a[minbc+2]
            error_type = a[minbc]
            if  error_type == "+":
                num_2=num_2+1
            list.append([error_lco,error_type,error_str])
            num_=+1
    return list"""


def align(a,b):
    error_list=[]
    for i in range(len(a)):
        if a[i] == b[i]:
            pass
        else :
            error_list.append((i,b[i]))
    return error_list


