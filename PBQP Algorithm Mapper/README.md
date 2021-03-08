# DYNAMAP
Dynamic Algorithm Mapping for CNN inference

## Software Algorithm Mapping

- System Requirement: See README in the root directory

### Software Algorithm Mapper User Guide

1. Run the CNN graph construction script:

```
    python cnn_const.py <Your_input_file>.in <Your_output_file>.in
```

Note that the <Your_input_file>.in file is an input config file where you can specify your CNN model metadata,
including the layer parameters (we use the notions in accordance with our paper, which is shown in the table below) and edge connections (i.e. layer ordering). 
| Feature map size      | Kernel size | # input channels |# output channels|
| :----:         |    :----:   |      :----:  | :----:  |
| H1,H2      | K1,K2       |  C_in  |  C_out  |

For each layer, specify the layer numbering (no requirement on specific order), the type of each layer (the program only considers the convolution layers "c" and depth/filter concatenation layers "d"), and two optional parameters: whether it is branching("fork")/merging("join") layer, and the outdegree(indegree) for the branching(merging) layers. 
Fpr each edge, specify the end-node numbers in a new line.
Our program performs node splitting and edge-reordering to automatically construct the PBQP graph and use the metadata to populate the node/edge cost vectors.

An example input file for GoogleNet inception module is included [incep_gn.in].

2. If everything is ok, the program will display "CNN GRAPH CONSTRUCTION COMPLETE" and drops an output <Your_output_file>.in file in the pbqp/testcases folder.

3. Go to directory pbqp/src and type

```
    make
    make install
    cd ..
```

4. Run the PBQP solver

```
    ./run_test
```

5. Dump files for the reduction (including topology graphs for intermediate steps) can be produced by 
```
    ./build_dump input_file
```
where input_file is the input file(s) in the test case directory. The dump can be found in the pbqp/dump folder. Open index.html to view the complete results



