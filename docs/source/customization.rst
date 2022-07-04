Customization
=============

You can find more custom to suit your needs in this section.

Startup Argvs
-------------

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


Customize Config
----------------

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

Customize Align
----------------

You can customize align.py and then modify the global matching algorithm. 
We recommend that you modify the global matching algorithm for better results before turning on the global matching feature.

The modified algorithm requires that the input is two sequences, returns a list with elements in tuple format, each tuple contains two elements, the position that does not match, and the base at that position in read_2.

Note: You need to change the now_align_alg in config_dict in load_config.py to True after modifying the global matching algorithm.

Customize Tree
--------------

We allow you to modify tree.py for more customization. Among other things you can modify the dna_dict in the trie class to allow Clover to handle DNA sequences that are not composed of ATGC. 
The format of the dictionary requires the key to be the type of base and the value to be a natural number starting from 0. There is no restriction on the exact order.