Package Structure
=================


The following figure shows how the package and its algorithms are structured:

.. code-block:: html

    ├── clover                                     // Source codes of Clvoer
    │    ├── __init__.py                           // Exhibition of class and method calls
    │    ├── align.py                              // Global Matching Module
    │    │    ├── global_align                     // Global comparison function
    │    ├── load_config.py                        // Import Parameters Module
    │    │    ├── load_json                        // Read parameter function
    │    │    ├── generate_vertical_drifts_list    // generate vertical drifts function
    │    │    ├── out_put_config                   // Parameter assignment function
    │    ├── main.py                               // Clover Main Module
    │    │    ├── MyProcess                        // Multi-process module
    │    ├── tree.py                               // Tree Structure Module
    │    │    ├── Trie                             // Tree structure class
    ├── example                                    // Example dataset for Clover
    │    ├── example_index_data.txt                // Example dataset without labels
    │    ├── example_tag_data.txt                  // Example dataset with labels
    ├── experiment result                          // Clover's experimental results
    │    ├── raw_data.xlsx                         // Experimental results table
    ├── tests                                      // Test module of source codes
    │    ├── test_align.py                         // Unit test for Global Matching
    │    ├── test_clust.py                         // Unit test for Clustering unit
    │    ├── test_tree.py                          // Unit test for Tree Structure
    ├── README.md                                  // Description document of library

The installation process using 'pip' only includes folder 'clover' and 'tests'.
