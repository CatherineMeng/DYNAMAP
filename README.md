# DYNAMAP
Dynamic Algorithm Mapping for CNN inference

## Software Algorithm Mapping

- System Requirement:
	- GCC compiler under linux
- Dependency:
	- Python 3.7
    - [Perl HTML Tag-Reader](http://pepper.linuxfocus.org/~guido/index.html#TagReader)
    - TeX system
    - [Tex to GIF](http://www.fourmilab.ch/webtools/textogif/textogif.html)
    - [Graphviz](http://www.graphviz.org/)
    - The core graph reduction functions are from the open-source [C-based Solver](http://www.complang.tuwien.ac.at/scholz/pbqp.html). Make sure to download the modified src in the Algorithm Mapper directory for complete dump file display.

### Software Algorithm Mapper User Guide

1. Download the Algorithm Mapper

2. Go to directory pbqp/src and type

```
    make
    make install
    cd ..
```
   
3. Run the CNN graph construction script:


```
    python cnn_const.py
```

- The program will ask for inputs about device capabilities, model name and CNN metadata. We provide preset config files for Inception modules.
- A config file with model name will be created in the "pbqp/testcases" folder.

4. Run the PBQP solver

```
    ./run_test
```

5. Dump files for the reduction (including topology graphs for intermediate steps) can be produced by 

```
    ./build_dump input_file
```

where input_file is the input file(s) in the test case directory. The dump can be found in the pbqp/dump folder. Open index.html to view the complete results

## Hardware Generation

- Dependency: 
	- VITIS HLS
	- Vivado 

