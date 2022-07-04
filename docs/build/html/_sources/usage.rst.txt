Usage
=====

This section is responsible for quickly using the Clover

Clustering labeled data
-------------------------

We provide an `example dataset < >` where the first column is label and the second column is sequence.

You can execute the following command:

.. code-block:: bash

   python -m clover.main -I example_tag_data.txt -L 150 -T 3

- **-I [input_file]** Input files, as the paper is currently under review, we only allow input in txt format.

- **-L [int]** Length of read.

- **-T [int]** The total number of clusters in the file, this selection can also be left out, although it will result in no coverage in the final statistics.


Clustering unlabeled data
-------------------------

We provide an `example dataset < >` where the first column is index and the second column is sequence.

You can execute the following command:

.. code-block:: bash

   python -m clover.main -I example_index_data.txt -O output_file -L 150 -P 0 --no-tag

- **-I [input_file]** Input file, as the paper is currently under review, we currently only allow input in txt,fasta,fastq

- **-O [output_file_name]** Output file name,Output file name.The output file will be located in the Clover folder

- **-P [int]** Select the process mode, if p=0 means single process operation. p=1, 2 means 4 processes, 16 processes, and so on,respectively. We do not recommend that p is greater than 2.