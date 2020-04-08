################################################################
# Project: Illumina reads evaluation, simulation and correction.
# Function: parameters and reads simulation
# Author: Peng Jia
# Email: pengjia@stu.xjtu.edu.cn
################################################################
import argparse
import os
import random
import subprocess
import numpy as np
import pandas as pd
import pysam
import uuid
from multiprocessing import Pool
from var import *
from sim import *
from eval import *

global paras


def arguments():
    """
    Deal with the command input and paramenter extraction.
    :return: class of parser containing input parameters.
    """
    commands = ["sim", "eval","var"]
    parser = argparse.ArgumentParser(description='PQ Tools.')
    parser.usage = " NGSsim.py <command> [options]"

    subparsers = parser.add_subparsers(title="command", metavar="", dest='command')
    parser_eval = subparsers.add_parser('eval', help='Evaluate the next generation sequencing data.')
    parser_eval.description = 'Evaluate the next generation sequencing data.'
    parser_eval.add_argument('-i', '--bam_file', required=True, type=str, nargs=1,
                             help="next generation sequencing bam file with samtools index! [required]")
    parser_eval.add_argument('-m', '--microsatellites', required=True, type=str, nargs=1,
                             help="microsatellites information in csv format. e.g. path/to.microsatellite.csv  [required]")
    parser_eval.add_argument('-o', '--output', required=True, type=str, nargs=1,
                             default=False, help="output file.")
    parser_eval.add_argument('-r', '--reference', required=True, type=str, nargs=1,
                             default=False, help="reference file.")
    parser_eval.add_argument('-oh', '--only_homopolymer', required=False, type=bool, nargs=1,
                             default=False, help="only evaluate polymerase slippage in homopolymer regions.")
    parser_sim = subparsers.add_parser('sim', help='Next generation sequencing reads simulated.')
    parser_sim.description = 'NGS simulation!'

    parser_sim.add_argument('-fa', '--fasta', required=True, type=str, nargs=1,
                            help="input genome.")
    parser_sim.add_argument("-fq ", "--fastq", required=False, type=str, nargs=1,
                            help="output fastq file.")
    parser_sim.add_argument("-r1 ", "--read1", required=False, type=int, nargs=1, default=[150],
                            help="first read length.")
    parser_sim.add_argument("-r2 ", "--read2", required=False, type=int, nargs=1, default=[150],
                            help="Second read length.")
    parser_sim.add_argument("-d ", "--depth", required=False, type=int, nargs=1, default=[30],
                            help="Second read length.")
    parser_sim.add_argument("-ismean", "--insertSzieMean", required=False, type=int, nargs=1, default=[500],
                            help="mean of the insert size")
    parser_sim.add_argument("-isstd", "--insertSzieStd", required=False, type=int, nargs=1, default=[100],
                            help="standard deviation of the insert size")
    parser_sim.add_argument("-cs", "--clusterSize", required=False, type=int, nargs=1, default=[100000],
                            help="molecular number of one sequencing cluster.")
    parser_sim.add_argument("-bs", "--bufSize", required=False, type=int, nargs=1, default=[1000],
                            help="buf size to write to the output of multi threads.")

    parser_sim.add_argument("-t ", "--threads", required=False, type=int, nargs=1, default=[1],
                            help="threads.")
    parser_var = subparsers.add_parser("var", help='Add variants to a genome.')
    parser_var.description = 'Add variants to a genome.'
    parser_var.add_argument('-i', "--input", required=True, type=str, nargs=1,
                            help="input genome (*.fa/*.fasta).")
    parser_var.add_argument("-o", "--output", required=True, type=str, nargs=1,
                            help="output genome(*.fa/*.fasta).")
    parser_var.add_argument("-vc", "--var_conf", required=True, type=str, nargs=1,
                            help="variants configure file (*.yaml).")
    parser_var.add_argument("-v", "--vcf", required=True, type=str, nargs=1, default=[],
                            help="output vcf file (*.vcf).")
    parser_var.add_argument("-c", "--conf", required=False, type=str, nargs=1,default=[],
                            help="configure file of the program (*.yaml).")

    if len(os.sys.argv) == 1 or os.sys.argv[1] not in commands:
        parser.print_help()
        parser.parse_args()
        return False
    # print("hhh")
    return parser




def initEval(parase):
    print()


def main():
    """
    Main function.
    :return:
    """
    arg = arguments()
    # print("here")
    # print(arg.parse_args())
    if arg:
        parase = arg.parse_args()
        if parase.command == "sim":
            sim(parase)
        elif parase.command == "eval":
            eval(parase)
        elif parase.command  == "var":
            var(parase)
            # eval(parase)


if __name__ == "__main__":
    main()
