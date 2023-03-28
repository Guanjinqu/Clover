"""Module for input parameters.

This module is responsible for importing the input parameters to Clover.

"""

import json
import getopt
import sys

config_dict={
    "read_len" : 152,
    "end_tree_len" : 15, 
    "other_tree_len" : 15,
    "other_tree_nums" : 2,
    "thd_tree_loc" : 40,
    "four_tree_loc" : 40,
    "Vertical_drift" : 2,
    "Horizontal_drift" : 3,
    "tree_threshold" : 10,
    "now_clust_threshold" : 8,
    "tag_nums" : 1,
    "processes_nums" : 0,
    "Cluster_size_threshold" : 1,
    "h_index_nums" : 0,
    "e_index_nums" : 0,
    "read_len_min" : 0,

    "align_fuc" : False ,
    "mmr_mode" : False ,
    "Virtual_mode" : True,
    "fast_mode" : True ,
    "tag_mode" : False,
    "Statistical_model" : False,
    "same_tree_len" : True,
    "now_align_alg" : False
}

opt,args = getopt.getopt(sys.argv[1:],'-I:-L:-D:-V:-H:-T:-P:-O:-h',['help','low','no-fast','no-tag','stat'])


#Read input info
def load_json(path):
    lines = []
    with open(path) as f :
        for row in f.readlines():
            if row.strip().startswith('//'):
                continue
            lines.append(row)
    return json.loads("\n".join(lines))

#Generate a sequence of vertical drifts
def generate_vertical_drifts_list(x):
    list=[]
    i=0
    k=-x-1
    while True:
        k=k+1
        list.append(k)
        i+=1
        if k == x :
            break
    return list

#Write the input to config.json
def out_put_config():
    opt,args = getopt.getopt(sys.argv[1:],'-I:-L:-D:-V:-H:-T:-P:-O:-h',['help','low','no-fast','no-tag','stat'])

    for opt_name,opt_value in opt :
        if '-h' in opt_name or '--help' in opt_name:
            print("Please see readme.md")
        if '-I' in opt_name :
            config_dict['input_path'] = opt_value
        if '-L' in opt_name :
            config_dict['read_len'] =int(opt_value)
        if '-D' in opt_name :
            config_dict['end_tree_len'] = int(opt_value)
            if config_dict['same_tree_len'] == True :
                config_dict['other_tree_len'] = int(opt_value)
        if '-V' in opt_name :
            config_dict['Vertical_drift'] = generate_vertical_drifts_list(int(opt_value))
        if '-H' in opt_name :
            config_dict['Horizontal_drift'] = int(opt_value)
        if  '-T' in opt_name :
            config_dict['tag_nums'] = int(opt_value)
            config_dict['tag_mode'] = True
        if '-P' in opt_name :
            config_dict['processes_nums'] = int(opt_value)
        if '-O' in opt_name :
            config_dict['output_file'] = opt_value+'.txt'
        if '--align' in opt_name:
            config_dict['align_fuc'] = True 
        if '--no-fast' in opt_name:
            config_dict['fast_mode'] = False
        if '--no-tag' in opt_name :
            config_dict['Virtual_mode'] = False
        if '--stat' in opt_name :
            config_dict['Statistical_model'] = True
        if '--low' in opt_name:
            config_dict['mmr_mode'] = True
            config_dict['fast_mode'] = False
            config_dict['align_fuc'] = False
            config_dict['Statistical_model'] = False
            config_dict['Virtual_mode'] = False
    if config_dict['read_len_min'] == 0 :
        config_dict['read_len_min'] = config_dict['read_len'] - 5
    if type(config_dict['Vertical_drift']) == int :
        config_dict['Vertical_drift'] = generate_vertical_drifts_list(config_dict['Horizontal_drift']) 
    config_dict['tag']="""
      ___           ___       ___           ___           ___           ___     
     /\  \         /\__\     /\  \         /\__\         /\  \         /\  \    
    /::\  \       /:/  /    /::\  \       /:/  /        /::\  \       /::\  \   
   /:/\:\  \     /:/  /    /:/\:\  \     /:/  /        /:/\:\  \     /:/\:\  \  
  /:/  \:\  \   /:/  /    /:/  \:\  \   /:/__/  ___   /::\~\:\  \   /::\~\:\  \ 
 /:/__/ \:\__\ /:/__/    /:/__/ \:\__\  |:|  | /\__\ /:/\:\ \:\__\ /:/\:\ \:\__\\
 \:\  \  \/__/ \:\  \    \:\  \ /:/  /  |:|  |/:/  / \:\~\:\ \/__/ \/_|::\/:/  /
  \:\  \        \:\  \    \:\  /:/  /   |:|__/:/  /   \:\ \:\__\      |:|::/  / 
   \:\  \        \:\  \    \:\/:/  /     \::::/__/     \:\ \/__/      |:|\/__/  
    \:\__\        \:\__\    \::/  /       ~~~~          \:\__\        |:|  |    
     \/__/         \/__/     \/__/                       \/__/         \|__|    """
    return config_dict


