[![PythonVersion](https://img.shields.io/badge/python-3.7-blue)](https://img.shields.io/badge/python-3.7-blue)

# Clover
Clover:  Tree structure-based efficient DNA clustering for DNA-based data storage


## What is Clover?
Clover is an efficient DNA sequence clustering algorithm, which applies to a large number of disordered DNA sequences generated after DNA sequencing in the DNA storage field. Clover has many advantages such as high efficiency, easy to use, easy to expand and low memory consumption. In our experiment, Clover could cluster 10 million sequences in about 10 seconds. In addition, we managed to cluster 10 billion DNA sequences using a home computer. Clover is written in Python, making it easy to extend(Although python is slower than C, Clover's algorithms are very fast, and even faster if you use pypy!). We need to emphasize that Clover is a clustering algorithm optimized specifically for data in the DNA storage domain, and clustering datasets from other domains is difficult to produce good results.Welcome to communicate with us if you have any suggestions or questions:D

## Requirements
- tqdm

## Quick Start
### Installation

Install it using pip:
> pip install dna-clover

Or you can also install it from source.

> git clone https://github.com/Guanjinqu/Clover.git

> cd Clover

> pip install -r requirements.txt

> python setup.py install develop --user

### Clustering labeled data:
> python -m clover.main -I example_index_data.txt -O output_file -L 150 -P 0 --no-tag
- **-I [input_file]** Input file, as the paper is currently under review, we currently only allow input in txt,fasta,fastq, and each line is 
> [index], [read]
- **-O [output_file_name]** Output file name,Output file name.The output file will be located in the Clover folder
- **-L [int]** Length of read
- **-P [int]** Select the process mode, if p=0 means single process operation. p=1, 2 means 4 processes, 16 processes, and so on,respectively. We do not recommend that p is greater than 2.
- **--no-tag** This selection needs to be added when the input file is in untagged format.
### Clustering labeled data:
> python -m clover.main -I example_tag_data.txt -L 150 -T 3
- **-I [input_file]** Input files, as the paper is currently under review, we currently only allow input in txt format, and each line is 
> [tag], [read]
- **-T [int]** The total number of clusters in the file, this selection can also be left out, although it will result in no coverage in the final statistics.
In this mode, we will output Clover's clustering statistics, such as elapsed time, accuracy, coverage (if the total number of tags is entered), redundancy, etc.
## Customization
### Startup Argvs

- **-I [input_file]** Input file
- **-O [output_file_name]** Output file name
- **-L [int]** Length of read
- **-P [int]** Select the process mode
- **-T [int]** The total number of clusters in the file
- **-D [int]** The depth of the tree, the higher the value the higher the accuracy of the clustering, but it will greatly increase the memory.
- **-V [int]** Vertical drift value, the higher the value the higher the accuracy, but it will increase the elapsed time.
- **-H [int]** Horizontal drift values, too large or small values will reduce the accuracy
- **--no-tag** Be sure to add this option if you use Clover for sequence clustering, which means that the input sequence is unlabeled.
- **--no-fast** If you don't have enough memory, you can add this option, which will reduce memory usage, but will increase the time consumption.
- **--low** If you need to cluster very large files (100 million+ sequences), you can add this option, which will enable the lowest memory usage mode, which will default to --no-tag and --no-fast and will output multiple parallel output files in multiple processes, which you can merge yourself.
- **--align** Adding this option will enable the global comparison feature. Turning it on will improve the clustering effect, but will seriously slow down the clustering speed. We allow you to customize the global matching algorithm, you just need to replace the align.py file and change 'now_align_alg' to true in config.json
- **--stat** Statistical mode, we allow to turn on statistical mode in --no-tag mode. After the clustering is finished, the statistics mode will give some feature statistics of the input file. We will provide the use of statistics mode after the paper is accepted

*Startup argvs will override the config file

### Customize Config
You can also modify the config_dict in load_config.py to get more custom options.
- read_len
  - number
  - Length of read
  - Default: 152
- end_tree_len
  - number
  - Depth of the tree at both ends
  - Default: 15
- other_tree_len
  - number
  - Depth of other trees
  - Defalut: 15
- other_tree_nums
  - number
  - Number of other trees
  - Default: 2
- thd_tree_loc
  - number
  - Location of the third tree
  - 40
- four_tree_loc
  - number
  - Location of the fourth tree
  - 40
- Vertical_drift
  - number
  - Vertical drift value
  - 2
- Horizontal_drift
  - number
  - Horizontal drift value
  - 3
- tree_threshold
  - number
  - Number of tree drift retrievals
  - 10
- now_clust_threshold
  - number
  - Drift threshold for new clusters
  - 8
- tag_nums
  - number
  - Number of clusters
  - 1
- processes_nums
  - number
  - Process mode
  - 0
- Cluster_size_threshold
  - number
  - Minimum number of elements in a cluster
  - 1
- h_index_nums
  - number
  - Number of bases in the leading primer
  - 0
- e_index_nums
  - number
  - Number of bases of back-end primers
  - 0
- read_len_min
  - number
  - Minimum processing length of read
  - 0
- align_fuc
  - boolean
  - Global Matching Mode
  - Default: false
- mmr_mode
  - boolean
  - Minimum memory usage mode
  - Default: false
- Virtual_mode
  - boolean
  - Input file with tags or not
  - Default: true
- fast_mode
  - boolean
  - High-speed mode (in this mode, the program will read the file into memory first)
  - Default: true
- tag_mode
  - boolean
  - Whether to enter a tag
  - Default: false
- Statistical_model
  - boolean
  - Feature statistics model
  - Default: false
- same_tree_len
  - boolean
  - The tree is the same length. If you want to customize the tree at the middle end, please modify this parameter.
  - Default: false
- now_align_alg
  - boolean
  - Whether to replace the global comparison algorithm, if so please modify this parameter.
  - Default: false

### Customize Align
You can customize align.py and then modify the global matching algorithm. 
We recommend that you modify the global matching algorithm for better results before turning on the global matching feature.

The modified algorithm requires that the input is two sequences, returns a list with elements in tuple format, each tuple contains two elements, the position that does not match, and the base at that position in read_2.

Note: You need to change the now_align_alg in config_dict in load_config.py to True after modifying the global matching algorithm.

### Customize Tree
We allow you to modify tree.py for more customization. Among other things you can modify the dna_dict in the trie class to allow Clover to handle DNA sequences that are not composed of ATGC. 
The format of the dictionary requires the key to be the type of base and the value to be a natural number starting from 0. There is no restriction on the exact order.

## License

Clover is licensed under the GNU General Public License, for more information read the LICENSE file or refer to:

http://www.gnu.org/licenses/

## Citation

Qu G, Yan Z, Wu H. Clover: tree structure-based efficient DNA clustering for DNA-based data storage[J]. Briefings in Bioinformatics, 2022, 23(5): bbac336.
