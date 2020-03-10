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

global paras


def arguments():
    """
    Deal with the command input and paramenter extraction.
    :return: class of parser containing input parameters.
    """
    commands = ["sim", "eval"]
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

    parser_sim.add_argument("-t ", "--threads", required=False, type=int, nargs=1,default=[1],
                            help="threads.")

    if len(os.sys.argv) == 1 or os.sys.argv[1] not in commands:
        parser.print_help()
        parser.parse_args()
        return False
    return (parser)

def eval(parase):
    """
    Evaluate the reads quality .
    :param parase:
    :return:
    """
    print()


def initEval(parase):
    print()


def initSimSeq(parase):
    """
    Initialize the parameters for reads simulations
    :param parase:
    :return: dict of simulations parameters
    """
    global paras
    paras = {}
    paras["fapath"] = parase.fasta[0]
    paras["fqpath"] = parase.fastq[0]
    paras["r1Len"] = parase.read1[0]
    paras["r2Len"] = parase.read2[0]
    paras["insertSzieMean"] = parase.insertSzieMean[0]
    paras["insertSzieStd"] = parase.insertSzieStd[0]
    paras["clusterSize"] = parase.clusterSize[0]
    paras["bufSize"] = parase.bufSize[0]
    paras["threads"] = parase.threads[0]
    paras["depth"] = parase.depth[0]

    paras["polymerase"]={}
    paras["polymerase"][paras["r1Len"]]=[ ((i+1)/(paras["r1Len"]))**2 for i in range(paras["r1Len"])]
    paras["polymerase"][paras["r2Len"]]=[ ((i+1)/(paras["r2Len"]))**2 for i in range(paras["r2Len"])]
    paras["polymerase"][paras["r1Len"]]=[ (1.1**(i+1))/(1.1**paras["r1Len"]) for i in range(paras["r1Len"])]
    paras["polymerase"][paras["r2Len"]]=[ (1.1**(i+1))/(1.1**paras["r2Len"]) for i in range(paras["r2Len"])]

    pa = 0.00001
    pc = 0.00001
    pg = 0.00001
    pt = 0.00001
    qa = 0.000001
    qc = 0.000001
    qg = 0.000001
    qt = 0.000001
    proP = {"A": pa, "C": pc, "G": pg, "T": pt}
    proQ = {"A": qa, "C": qc, "G": qg, "T": qt}
    paras["p"] = proP
    paras["q"] = proQ
    return paras


def getReverseComplematary(read):
    """
    Get the reverse complementary sequencing of given reads.
    :param read:
    :return:
    """
    reverseComplematary = ""
    for ipos in read[::-1]:
        if ipos == "A":
            reverseComplematary += "T"
        elif ipos == "T":
            reverseComplematary += "A"
        elif ipos == "G":
            reverseComplematary += "C"
        elif ipos == "C":
            reverseComplematary += "G"
    return reverseComplematary


def sequencing(read, readLen=150):
    """
    Read sequencing simulation.
    :param read: True read sequence (template)
    :param readLen: read length
    :param para: parameter of simulation.
    :return: sequencing read and qualities
    """
    # print("read",len(read))
    read = read + "N"
    readstrindexmax = len(read)
    curent = np.zeros(paras["clusterSize"])
    curent = pd.Series(curent, dtype="int")
    query = ""
    qualList = []
    # print(curentSet)
    for pos in range(readLen):
        curentSet = curent.value_counts()
        qualPro = curent.value_counts(normalize=True)
        base = read[np.random.choice(qualPro.index, p=qualPro.iloc[:].ravel())]
        query += base
        qualList.append(qualPro.max())
        # print("curentSet",curentSet)
        curent_new = np.array([])
        for ccc in curentSet.index:
            cccNum = curentSet[ccc]
            basenow = read[ccc]
            pdel = paras["p"][basenow]*paras["polymerase"][readLen][pos]
            pins = paras["q"][basenow]*paras["polymerase"][readLen][pos]
            pnorm = 1 - pdel - pins
            psteps = np.array([pdel, pnorm, pins])
            steps = np.random.choice([0, 1, 2], cccNum, p=psteps.ravel())
            curent_new = np.hstack([curent_new, steps + ccc])
        curent = pd.Series(curent_new, dtype="int")
        curent[curent > readstrindexmax] = readstrindexmax

    def _pro2basequality(qual):
        """
        Converse sequencing quality from float to phase33
        :param qual:
        :return:
        """
        if qual > 0.99999:
            qual = 0.99999
        return chr(int(-10 * np.log10(1 - qual)) + 33)

    quals = "".join(list(map(_pro2basequality, qualList)))
    return query, quals


def _processFragments(fragment):
    """
    Process Fragments for paired read sequencing
    :param contig:
    :param fragment:
    :param fafile:
    :return:
    """
    r1Len = paras["r1Len"]
    r2Len = paras["r2Len"]
    fafile = pysam.FastaFile(paras["fapath"])
    contig = fragment[-1]
    if fragment[2] == 1:
        leftstart = fragment[0]
        # leftend=fragment[0]+r1Len
        # rightstart=fragment[1]-r2Len
        rightend = fragment[1]
        read1 = fafile.fetch(contig, leftstart, leftstart + 2 * r1Len).upper()
        read2 = getReverseComplematary(fafile.fetch(contig, rightend - 2 * r2Len, rightend).upper())
    else:
        leftstart = fragment[0]
        # leftend = fragment[0] + r2Len
        # rightstart = fragment[1] - r1Len
        rightend = fragment[1]
        read2 = fafile.fetch(contig, leftstart, leftstart + 2 * r2Len).upper()
        read1 = getReverseComplematary(fafile.fetch(contig, rightend - 2 * r1Len, rightend).upper())
    fq1str = ""
    fq2str = ""
    if _checkstr(read1) and _checkstr(read2):
        query1, qual1 = sequencing(read1, r1Len)
        query2, qual2 = sequencing(read2, r2Len)
        readname = uuid.uuid1()
        # fq1.write("@" + str(readname) + "\n" + query1 + "\n+" + "\n" + qual1 + "\n")
        # fq2.write("@" + str(readname) + "\n" + query2 + "\n+" + "\n" + qual2 + "\n")
        fq1str = "@" + str(readname) + "\n" + query1 + "\n+" + "\n" + qual1 + "\n"
        fq2str = "@" + str(readname) + "\n" + query2 + "\n+" + "\n" + qual2 + "\n"
    # print(fq1str)
    return (fq1str, fq2str)


def _checkstr(read):
    """
    Remove sequence with N
    :param read:
    :return:
    """
    readChr = set(read)
    chrNum = 0
    for chr in ["A", "T", "C", "G"]:
        if chr in readChr:
            chrNum += 1
    if chrNum < len(readChr):
        return False
    else:
        return True


def buildLib(fa, depth=30, r1Len=150, r2Len=150, insertsizeMean=500,
             inserSizeStd=100, fq1="", fq2=""):
    """
    :param fa: reference
    :param depth: sequencing depth (e.g. 30x)
    :param r1Len: R1 length
    :param r2Len: R2 length
    :param insersizeMean:  mean of insertszie
    :param inserSizeStd:  standard divation of insertsize
    :param paras: simulation parameters
    :param fq1: output R1
    :param fq2: output R2
    :return:
    """
    fafile = pysam.FastaFile(fa)
    contigs = fafile.references

    def _processOneContig(contig):
        """
        Process framgents of one contig in reference.
        :param contig:
        :return:
        """
        fafile = pysam.FastaFile(fa)
        contigLen = fafile.get_reference_length(contig)
        fragmentsNum = int(contigLen * depth / (r1Len + r2Len))
        # print("fra", fragmentsNum)
        fragments = []
        fragmentsCureentNum = 0
        while True:
            insertSize = random.normalvariate(insertsizeMean, inserSizeStd)
            position = np.random.randint(0, contigLen - 1)
            direction = np.random.choice([-1, 1])
            left = position
            right = position + insertSize

            if right > contigLen - 1:
                continue
            elif left < 0:
                continue
            else:
                fragmentsCureentNum += 1
                fragments.append([left, right, direction, contig])
            if fragmentsCureentNum == 200000:
                break
        return fragments

    totalfragments = {}
    for contig in contigs:
        # print("Process Contig " + contig)
        totalfragments[contig] = _processOneContig(contig)
    pool = Pool(paras["threads"])
    for contig in totalfragments:
        splitnum = paras["bufSize"]
        fragmentNum = 0
        totalfragmentsNum = len(totalfragments[contig])
        for partfrangments in [totalfragments[contig][i:i + splitnum] for i in
                               range(0, len(totalfragments[contig]), splitnum)]:
            fragmentNum += len(partfrangments)

            partresult = pool.map(_processFragments, partfrangments)
            # print(fragmentNum)
            # print(partresult)
            print("[Info] Processing: " + contig, " Total: " + str(totalfragmentsNum), "Finish: " + str(fragmentNum))
            for readstr in partresult:
                fq1.write(readstr[0])
                fq2.write(readstr[1])


def sim(parase):
    """
    Simulation of reads
    :param parase:
    :return:
    """
    # parase.print_help()
    initSimSeq(parase)
    # reffa = parase.fasta[0]
    # print(parase.fastq)
    fq1 = open(parase.fastq[0] + "_R1.fq", "w")
    fq2 = open(parase.fastq[0] + "_R2.fq", "w")
    buildLib(fa=paras["fapath"], depth=paras["depth"], r1Len=paras["r1Len"], r2Len=paras["r2Len"],
             insertsizeMean=paras["insertSzieMean"],
             inserSizeStd=paras["insertSzieStd"], fq1=fq1, fq2=fq2)
    fq1.close()
    fq2.close()


def main():
    """
    Main function.
    :return:
    """
    arg = arguments()
    if arg:
        parase = arg.parse_args()
        if parase.command == "sim":
            sim(parase)
        else:
            eval(parase)


if __name__ == "__main__":
    main()
