#!/usr/bin/env python
"""
 A HMM-based predictor for TA proteins, implemented in python.

TAPPM predicts TA proteins solely from target amino acid sequences according to
the knowledge of the sequence features of TMDs and the peripheral regions of
TA proteins. The hidden markov models of TA proteins as well as three different
types of transmembrane proteins with similar structures and compared their
likelihoods as TA proteins. Using these models, TAPPM achieved high prediction
accuracy (area under the receiver operator curve values reaching 0.963).

A web application of TAPPM is also available freely at
http://tenuto.bi.a.u-tokyo.ac.jp/tapp/

Contact info
---------------------
By Post mail:

Bioinformation Engineering Laboratory
Graduate School of Agricultural and Life Science
University of Tokyo.
No.1-1-1 Yayoi Bukyo-ku, Tokyo, Japan. ZIP:113-8657


By email:

Mr. Shunsuke Shigemitsu
choge@choge.net
or
Dr. Wei CAO
davecao@bi.a.u-tokyo.ac.jp

"""
import os
import sys
# import json
# import numpy as np

from signal import signal, SIGPIPE, SIG_DFL
from multiprocessing import cpu_count
from optparse import OptionParser, OptionGroup

from tappm import FastaReader, FastaBuilder, MyHmmPredictor, MODELPATH
from tappm.fasta import SwissProt, TrEMBLE, GenBank_refseq, BasicProteinFasta
from tappm.utils import Console, IndentedHelpFormatterWithNL,\
                        convert_numpy_types
from tappm.io import report

__date__ = "2016/01/16"
__version__ = "1.0.0"

# reset
signal(SIGPIPE, SIG_DFL)
# global
LOGGER = None

FORMAT = {
    'swissprot': SwissProt,
    'tremble': TrEMBLE,
    'genebank': GenBank_refseq,
    'free': BasicProteinFasta
}

MODELS = {
    'TA': MODELPATH + os.sep + 'ta4.xml',
    'MP': MODELPATH + os.sep + 'mp.xml',
    'TACODE': 'TTHHHHHHHHHHHHHHHHHHHHHHHHHCCCCCGTT',
    'MPCODE': 'SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSG'
              'LLLLLLLLLLLLLLLLLLLLCCCCCHHHHHHHHHHHHHHHHHHHHHHHHH'
}


def parse_cmd(argv):
    """
        Parse command line arguments
    """
    usage = 'usage: %prog [optioins] -i xx.fasta'
    parser = OptionParser(
                formatter=IndentedHelpFormatterWithNL(),
                add_help_option=True, usage=usage, version=__version__)
    # --- 1. general options ---
    general_opts = OptionGroup(parser, "General options")
    general_opts.add_option(
        "--log", dest='logfilename', default='tappm_cli',
        type='string',
        help="The name of a log file. If not specified, the name will" +
             " be tappm_cli.log"
    )
    general_opts.add_option(
        "--log-level", dest='log_level', default='info', type='choice',
        choices=['debug', 'info', 'warnings', 'error', 'critical', 'none'],
        help="The name of a log file. If not specified, the name will" +
             " be tappm_cli.log"
    )
    general_opts.add_option(
        "-v", "--verbose",
        action="store_true", dest="verbose", default=False,
        help="Show verbose info"
    )
    # --- 1. input and output ---
    inout_opts = OptionGroup(parser, "Input/output settings")
    inout_opts.add_option(
        "-i", "--inputfile", dest="fastafile",
        help="The input file in fasta format[REQUIRED]."
    )
    inout_opts.add_option(
        "--fmt", dest="dbformat", default="free", type='choice',
        choices=['swissprot', 'tremble', 'genebank', 'free'],
        help="The format of the input file, default 'free'."
    )
    # inout_opts.add_option(
    #     "-m", "--modelfile", dest="modelfile",
    #     help="The HMM model file[REQUIRED]."
    # )

    inout_opts.add_option(
        "--outfmt", dest="outfmt", type='choice',
        choices=['tabular', 'text', 'xml'], default='tabular',
        help="Specify the format of output. xml', tabular' and 'text'"
             " are supported now, and corresponding suffixes are xml, tab "
             "and text, respectively."
    )
    inout_opts.add_option(
        "--outdir", dest="outdir", type='string',
        default='.',
        help="Specify the output directory, default is '.'."
    )
    inout_opts.add_option(
        "--out", dest="outfile", type='string',
        default='tappm_report',
        help="Specify the output file name, default is 'tappm_report.tabular'."
    )

    inout_opts.add_option(
        "-t", "--threshold", dest="threshold", type='float',
        default=-0.0167222981,
        help="The threshold value used for"
             "final decision to predict TA or not. default is -0.0167222981."
    )
    # --- 2. Multithreads options ---
    multihtreads_opts = OptionGroup(
         parser,
         "Multithreads options for loading data",
         "Option arguments for multithreads."
    )

    multihtreads_opts.add_option(
        "--mcpu", dest="mcpu", type='int', default=None,
        help="The number of threads, default is the number of cores."
    )
    parser.add_option_group(general_opts)
    parser.add_option_group(inout_opts)
    parser.add_option_group(multihtreads_opts)

    options, arguments = parser.parse_args(argv)
    if arguments == 0:
        print ("Error: no arguments found")
        parser.print_help()
        sys.exit(1)

    # check input fasta files
    if not options.fastafile:
        print ("Error: do not specify an input file")
        parser.print_help()
        sys.exit(1)

    return options


def main(argv):

    global LOGGER
    # parse command line arguments
    cmds = ' '.join(argv)
    opt = parse_cmd(argv)
    # general opts
    verbose = opt.verbose
    # input and output
    inputfile = opt.fastafile
    dbformat = FORMAT[opt.dbformat]
    threshold = opt.threshold
    outfmt = opt.outfmt
    outputfile = opt.outdir + os.sep + opt.outfile + '.' + outfmt

    # input file path, name, and ext
    inputfilename = os.path.basename(inputfile)
    mcpu = opt.mcpu if opt.mcpu else cpu_count()
    # Logger settings
    logfile = opt.logfilename
    log_label = inputfilename
    LOGGER = Console(log_label, prefix='@>', console=opt.log_level)
    # LOGGER = Console(log_label, prefix='@>')
    LOGGER.start(logfile)

    # Load fasta file
    reader = FastaReader(FastaBuilder(dbformat), protein=True)
    fasta_list = reader.parse_file(inputfile)
    if verbose:
        print("Input file: {}".format(inputfile))
        print("Read {} sequences".format(len(fasta_list)))
        print("Paramters:")
        print("  threshold:{:10.6f}".format(threshold))
        print("  nCPU: {}".format(mcpu))
        print("  Model:{}".format(MODELPATH))
    LOGGER.info("cmd: {}".format(cmds))
    LOGGER.info("Input file: {}".format(inputfile))
    LOGGER.info("Read {} sequences".format(len(fasta_list)))
    LOGGER.info("Paramters:")
    LOGGER.info("  threshold:{:10.6f}".format(threshold))
    LOGGER.info("  nCPU: {}".format(mcpu))
    # MyHmmPredictor
    # TA prediction
    ta_predictor = MyHmmPredictor(filename=MODELS['TA'], cpus=mcpu)
    ta_predictor.set_decoder(MODELS['TACODE'])
    # MP prediction
    mp_predictor = MyHmmPredictor(filename=MODELS['MP'], cpus=mcpu)
    mp_predictor.set_decoder(MODELS['MPCODE'])
    # print MODELS['MPCODE']
    LOGGER.info("Start to scan sequences.")
    LOGGER.timeit(label='prediction')
    prediction_ta = ta_predictor.predict(fasta_list, reverse=True)
    prediction_mp = mp_predictor.predict(fasta_list)
    resultList = convert_numpy_types(
                    prediction_ta,
                    prediction_mp,
                    threshold,
                    outfmt,
                    fasta_list)
    LOGGER.report(msg='Completed in %.2fs', label='prediction')
    if verbose:
        print("Prepare output")
    LOGGER.info("Prepare output")
    LOGGER.timeit(label='output')
    report(resultList, outputfile, fmt=outfmt, backupCount=5,
           threshold=threshold, totalSeq=len(fasta_list))
    LOGGER.report(msg='Completed in %.2fs', label='output')
    if verbose:
        print("Finished")


if __name__ == '__main__':
    main(sys.argv)
