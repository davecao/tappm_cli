TAPPM package
=============

A HMM-based predictor for TA proteins, implemented in python.

 TAPPM predicts TA proteins solely from target amino acid sequences according to 
 the knowledge of the sequence features of TMDs and the peripheral regions of 
 TA proteins. The hidden markov models of TA proteins as well as three different
 types of transmembrane proteins with similar structures and compared their 
 likelihoods as TA proteins. Using these models, TAPPM achieved high prediction 
 accuracy (area under the receiver operator curve values reaching 0.963). 

The web application of TAPPM is also available freely at 
http://tenuto.bi.a.u-tokyo.ac.jp/tapp/

Prerequisites
-------------

1. Numpy
2. jinja2

Download and Installation
-------------------------
1. Numpy

    git clone git://github.com/numpy/numpy.git numpy  
    cd numpy  
    python setup.py build install  

2. tappm

    cd tappm  
    python setup.py build install_dist  

Options
--------------
Usage: tappm_cli.py [optioins] -i xx.fasta

Options:
    --version             show program's version number and exit
    -h, --help            show this help message and exit

  General options:
    --log=LOGFILENAME   The name of a log file. If not specified, the name
                        will be tappm_cli.log
    --log-level=LOG_LEVEL
                        The name of a log file. If not specified, the name
                        will be tappm_cli.log
    -v, --verbose       Show verbose info

  Input/output settings:
    -i FASTAFILE, --inputfile=FASTAFILE
                        The input file in fasta format[REQUIRED].
    --fmt=DBFORMAT      The format of the input file, default 'free'.
    --outfmt=OUTFMT     Specify the format of output. xml', tabular' and
                        'text' are supported now, and corresponding suffixes
                        are xml, tab and text, respectively.
    --outdir=OUTDIR     Specify the output directory, default is '.'.
    --out=OUTFILE       Specify the output file name, default is
                        'tappm_report.tabular'.
    -t THRESHOLD, --threshold=THRESHOLD
                        The threshold value used forfinal decision to predict
                        TA or not. default is -0.0167222981.

  Multithreads options for loading data:
    Option arguments for multithreads.

    --mcpu=MCPU         The number of threads, default is the number of cores.
