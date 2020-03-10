import argparse
import os
import random
import subprocess

import numpy as np
import pandas as pd
import pysam
import uuid


def arguments():
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
    parser_sim.add_argument("-t ", "--threads", required=False, type=str, nargs=1,
                            help="threads.")

    # print(parser.error())
    # print(parser_eval)
    if len(os.sys.argv) == 1:
        parser.print_help()
        return False
    return (parser)

#
# def bamEvalue(bampath, bedpath, refpath, outpath, thread=1, tmpPath=os.getcwd(), prefix=15, suffix=15):
#     bedfile = open(bedpath, "r")
#     regions = [line[:-1].split("\t") for line in bedfile]
#     ref = pysam.FastaFile(refpath)
#     samfile = pysam.AlignmentFile(bampath, "rb")
#     try:
#         samfile.check_index()
#     except Exception as msg:
#         print(msg)
#         print("[WARNING] Bam file has no index, doing index now ...")
#         subprocess.call('samtools index -@ ' + str(thread) + ' ' + bampath, shell=True)
#     else:
#         print("[INFO] Check bam index: OK!")
#     file = open(outpath, "w")
#     file.write(",".join(
#         ["chrom", "start", "end"] + [i + "_" + j + "_" + k for i in ["prefix", "MS", "suffix"] for j in
#                                      ["forward", "reverse"] for k in
#                                      ["base", "mismatch", "deletion", "insertion"]]) + "\n")
#     rnum = 0
#     total = len(regions)
#     for region in regions:
#         rnum += 1
#         print(rnum, total, region)
#         myregion = region
#         chr = str(region[0])
#         start = int(region[1])
#         end = int(region[2])
#         prefixPos = start - prefix
#         suffixPos = end + suffix
#         ThisReads = samfile.fetch(chr, start, end)
#         tmpBamPathForward = tmpPath + "/" + chr + "_" + str(start) + "_" + str(end) + "_Forward.bam"
#         tmpBamPathReverse = tmpPath + "/" + chr + "_" + str(start) + "_" + str(end) + "_Reversed.bam"
#         tmpFarward = pysam.AlignmentFile(tmpBamPathForward, "wb", template=samfile)
#         tmpReverse = pysam.AlignmentFile(tmpBamPathReverse, "wb", template=samfile)
#         for read in ThisReads:
#             if read.is_reverse:
#                 tmpReverse.write(read)
#             else:
#                 tmpFarward.write(read)
#         tmpFarward.close()
#         tmpReverse.close()
#         subprocess.call('samtools index -@ ' + str(thread) + ' ' + tmpBamPathForward, shell=True)
#         subprocess.call('samtools index -@ ' + str(thread) + ' ' + tmpBamPathReverse, shell=True)
#         tmpFarward = pysam.AlignmentFile(tmpBamPathForward, "rb")
#         tmpReverse = pysam.AlignmentFile(tmpBamPathReverse, "rb")
#         # print(tmpReverse.mapped)
#         if tmpReverse.mapped < 10 or tmpFarward.mapped < 10:
#             tmpFarward.close()
#             tmpReverse.close()
#             subprocess.call('rm ' + tmpBamPathForward + "*", shell=True)
#             subprocess.call('rm ' + tmpBamPathReverse + "*", shell=True)
#             continue
#         tmpsam = {"forward": tmpFarward, "reverse": tmpReverse}
#         # tmpsam={"reverse":tmpReverse}
#         thisRegions = {"prefix": [prefixPos, start], "ms": [start, end], "suffix": [end, suffixPos]}
#
#         print(tmpFarward)
#         # for reg in thisRegions:
#         #     startP = thisRegions[reg][0]
#         #     endP = thisRegions[reg][1]
#         #     refSeq = ref.fetch(chr)[startP:endP].upper()
#         #     for direction in tmpsam:
#         #         pos = 0;
#         #         deletion = 0;
#         #         insertion = 0;
#         #         mismatch = 0
#         #         base = 0
#         #         # pileups = tmpsam[direction].pileup(chr, startP, endP, truncate=True, stepper="all", min_base_quality=0)
#         #         # for pileupPos in pileups:
#         #         #     pos += 1
#         #         #     query = [ch.upper() for ch in pileupPos.get_query_sequences(add_indels=True)]
#         #         #     base += len(query)
#         #         #     refPos = refSeq[pos - 1]
#         #         #     for posQuery in query:
#         #         #         if "*" == posQuery[0]:
#         #         #             continue
#         #         #         # else:
#         #         #         #     if posQuery[0] != refPos:
#         #         #         #         mismatch += 1
#         #         #         if len(posQuery) > 3:
#         #         #             if posQuery[1] == "+":
#         #         #                 insertion += int(posQuery[2])
#         #         #             elif posQuery[1] == "-":
#         #         #                 deletion += int(posQuery[2])
#         #         myregion.append(base)
#         #         myregion.append(mismatch)
#         #         myregion.append(deletion)
#         #         myregion.append(insertion)
#         tmpFarward.close()
#         tmpReverse.close()
#         subprocess.call('rm ' + tmpBamPathForward + "*", shell=True)
#         subprocess.call('rm ' + tmpBamPathReverse + "*", shell=True)
#         file.write(",".join(map(str, myregion)) + "\n")
#     file.close()
#     samfile.close()
#
#
# def argumentProcress():
#     """
#     argument procress
#     """
#     global arguments
#     parser = argparse.ArgumentParser(description='MSIHunter: a Microsatellite Instability(MSI)'
#                                                  ' detection tools using only tumor sequencing data!\n'
#                                                  'You can test multiple sample one time in this tools')
#     parser.add_argument('-i', '--input_configure', required=True, type=str, nargs=1,
#                         help="The path of input configure files [required]")
#     parser.add_argument('-o', '--workspace', required=True, type=str, nargs=1, default=["./workspace"],
#                         help="prefix of the output [required]")
#     parser.add_argument('-mc', '--microsatellites_configure', required=True, type=str, nargs=1, default=["NA"],
#                         help="path of the microsatellites configure files [required]")
#     parser.add_argument('-t', '--threads', type=int, nargs=1, default=[4],
#                         help="mumber of additional threads to use [default:4]")
#     parser.add_argument('-q', '--minimum_mapping_quality', type=int, nargs=1, default=[20],
#                         help="minimum mapping quality of read [default:20]")
#     parser.add_argument('-s', '--minimum_support_reads', type=int, nargs=1, default=[20],
#                         help="minimum support reads of an available microsatellite [default:20]")
#     args = parser.parse_args()
#     arguments = {}
#     arguments["input"] = args.input_configure[0]
#     arguments["workspace"] = args.workspace[0] if args.workspace[0][-1] == "/" else args.workspace[0] + "/"
#     arguments["Microsatellite"] = args.microsatellites_configure[0]
#     arguments["threads"] = args.threads[0]
#     arguments["minimum_support_reads"] = args.minimum_support_reads[0]
#     arguments["minimum_mapping_quality"] = args.minimum_mapping_quality[0]
#     ErrorStat = False
#     if os.path.isfile(arguments["input"]):
#         print("[INFO] The train set is : " + arguments["input"])
#     else:
#         print('[ERROR] The train set "' + arguments["input"] + '" is not exist, please check again')
#         ErrorStat = True
#     if os.path.isfile(arguments["Microsatellite"]):
#         print("[INFO] The Microsatellites file  is : " + arguments["Microsatellite"])
#     else:
#         print('[ERROR] The Microsatellites file "'
#               + arguments["Microsatellite"]
#               + '" is not exist, please check again')
#         ErrorStat = True
#     if os.path.exists(arguments["workspace"]):
#         print(
#             '[ERROR] The workspace is still exist! in case of overwrite files in this workspace, '
#             'please give a new work space')
#         ErrorStat = True
#         if ErrorStat: return False
#     else:
#         if ErrorStat: return False
#         os.mkdir(arguments["workspace"])
#         os.mkdir(arguments["workspace"] + "detailInfo/")
#         print("[INFO] The workspace is : " + arguments["workspace"])
#     return True
#
#
# def loadMicroSatellite(microsatellitesFile):
#     """
#     :return:
#     """
#     global dfMicroSatellites, arguments
#     dfMicroSatellites = {}
#     if not os.path.isfile(microsatellitesFile):
#         print("[ERROR]: The microsatellite file '" + microsatellitesFile + "' is not exist!")
#     return False
#     microsatele = pd.read_csv(microsatellitesFile, index_col=0)
#
#     columns = ["chr", "pos", "motif", "motifLen", "repeatTimes", "prefix", "suffix"]
#     errorItem = []
#     for item in columns:
#         if item not in dfMicroSatellites[cancer].columns:
#             errorItem.append(item)
#     if len(errorItem) > 0:
#         print("[ERROR] The item " + ",".join(
#             errorItem) + " of the Micorsatellites file is not exist! Please check again!")
#         return False
#     return microsatele
#
#
# def getRepeatTimes(alignment, motif, motifLen, prefix, suffix):
#     """
#     :param alignment:
#     :param motif:
#     :param motifLen:
#     :param prefix:
#     :param suffix:
#     :return:
#     """
#     global arguments
#     if alignment.mapping_quality < arguments["minimum_mapping_quality"]:
#         return -1
#     readString = alignment.query
#     prefixState = readString.find(prefix)
#     if prefixState < 0: return -1
#     suffixState = readString.rfind(suffix)
#     if suffixState < 0: return -3
#     if prefixState + 5 >= suffixState: return -2
#     while prefixState >= 0:
#         count = 0
#         start = prefixState + 5
#         while start == readString.find(motif, start):
#             count += 1
#             start = readString.find(motif, start) + motifLen
#         if (motifLen == 1 and count >= 1) or (motifLen > 1 and count >= 1):
#             if start == readString.find(suffix, start):
#                 # print(count, "    ", prefix,motif, suffix,repeat)
#                 return count
#         prefixState = readString.find(prefix, prefixState + 1)
#     return -4
#
#
# def calcuShiftProbability(disDict, refRepeatTimes):
#     """
#     :param disDict:
#     :param refRepeatTimes:
#     :return:
#     """
#     insShfit = 0;
#     delShfit = 0;
#     normal = 0
#     for rpt in disDict:
#         if rpt - refRepeatTimes > 0:
#             insShfit = insShfit + (rpt - refRepeatTimes) * disDict[rpt]
#             normal = normal + rpt * disDict[rpt] - (rpt - refRepeatTimes) * disDict[rpt]
#         else:
#             delShfit = delShfit + (refRepeatTimes - rpt) * disDict[rpt]
#             normal = normal + rpt * disDict[rpt] - (refRepeatTimes - rpt) * disDict[rpt]
#     return round(delShfit / (insShfit + delShfit + normal), 4), round(insShfit / (insShfit + delShfit + normal), 4)
#
#
# def procressInputConfigure():
#     """
#     :return:
#     """
#     global arguments, dfMicroSatellites, caseInfo, Result
#
#     print("[INFO] Loading the input configure file from " + arguments["input"] + " ...")
#     caseInfo = pd.read_csv(arguments["input"], index_col=0, dtype="str")
#     columns = ["bamPath", "cancerType"]
#     errorItem = []
#     for item in columns:
#         if item not in caseInfo.columns:
#             errorItem.append(item)
#     if len(errorItem) > 0:
#         print("[ERROR] The item " + ",".join(
#             errorItem) + " of the input configure file is not exist! Please check again!")
#         return False
#     if caseInfo.isnull().any().any():
#         print("[ERROR] The information of your input configure file is incomplete,please check again ")
#         return False
#
#     cancerTypeList = set(caseInfo["cancerType"])
#     for cancer in cancerTypeList:
#
#         if cancer not in dfMicroSatellites:
#             print("[ERROR] There is no microsatellite file for " + cancer + ",please check again ")
#     for id, info in caseInfo.iterrows():
#         if not os.path.isfile(info['bamPath']):
#             print("[ERROR] This file " + info["bamPath"] + "is not exist! Please check it again!")
#             return False
#     Result = caseInfo
#     return True
#

# def procressOneBam(caseid):
#     """
#     :param caseid:
#     :return:
#     """
#     global caseInfo, arguments, dfMicroSatellites
#     try:
#         infoTest = caseInfo.loc[caseid, :]
#         lociRes = pd.DataFrame()
#         # print(infoTrain)
#         print("[INFO] procressing " + str(caseid) + "...")
#         os.mkdir(arguments["workspace"] + "detailInfo/" + str(caseid))
#         file = open(arguments["workspace"] + "detailInfo/" + str(caseid) + "/" + str(caseid) + ".dis", "w")
#         file.close()
#         cancer = infoTest.loc["cancerType"]
#         bam = pysam.AlignmentFile(infoTest["bamPath"], "rb")
#         Distribution = {}
#         curentMSNum = 0
#         for id, info in dfMicroSatellites[cancer].iterrows():
#             chrId = info["chr"]
#             posStart = info["pos"]
#             motifLen = int(info["motifLen"])
#             motif = info["motif"]
#             repeatTimes = int(info["repeatTimes"])
#             prefix = info['prefix']
#             suffix = info['suffix']
#             posEnd = posStart + motifLen * repeatTimes
#             queryStrat = posStart - 5
#             queryEnd = posEnd + 5
#             alignmentList = []
#             for alignment in bam.fetch(chrId, queryStrat,
#                                        queryEnd):
#                 # (refName,start,end): read which at least has a base between  start+1 and end-1
#                 alignmentList.append(alignment)
#             if len(alignmentList) < arguments["minimum_support_reads"]: continue
#             repeatTimesDict = {}
#             for alignment in alignmentList:
#                 if alignment.is_unmapped: continue
#                 thisRepeatTimes = getRepeatTimes(alignment, motif, motifLen, prefix, suffix)
#                 if thisRepeatTimes < 0: continue
#                 if thisRepeatTimes not in repeatTimesDict: repeatTimesDict[thisRepeatTimes] = 0
#                 repeatTimesDict[thisRepeatTimes] += 1
#             if sum(repeatTimesDict.values()) < arguments["minimum_support_reads"]:
#                 continue
#             else:
#                 curentMSNum += 1
#                 Distribution[id] = repeatTimesDict
#                 series = dfMicroSatellites[cancer].loc[id]
#                 series["threshold"] = round(series["threshold"], 4)
#                 series["p"] = round(calcuShiftProbability(repeatTimesDict, repeatTimes)[0], 4)
#                 lociRes = lociRes.append(series)
#             if curentMSNum % 1000 == 0:
#                 file = open(arguments["workspace"] + "detailInfo/" + str(caseid) + "/" + str(caseid) + ".dis", "a")
#                 for id, dis in Distribution.items():
#                     file.write(str(id) + "\n" +
#                                " ".join([str(key) + ":" + str(dis[key]) for key in sorted(list(dis.keys()))]) + "\n")
#                 file.close()
#                 Distribution = {}
#         if len(Distribution) > 0:
#             file = open(arguments["workspace"] + "detailInfo/" + str(caseid) + "/" + str(caseid) + ".dis", "a")
#             for id, dis in Distribution.items():
#                 file.write(str(id) + "\n" +
#                            " ".join([str(key) + ":" + str(dis[key]) for key in sorted(list(dis.keys()))]) + "\n")
#             file.close()
#         lociRes = lociRes[
#             ["chr", "pos", "motif", "motifLen", "repeatTimes", "prefix", "suffix", "threshold", "p"]]
#         lociRes.index.name = "id"
#         fileRes = open(arguments["workspace"] + "detailInfo/" + str(caseid) + "/" + str(caseid), "w")
#         totalNum = len(lociRes)
#         InstableNum = len(lociRes[lociRes["p"] > lociRes["threshold"]])
#         InstablePercentage = InstableNum / totalNum * 100
#
#         # Result.loc[caseid,"percentageOfInstableMicosatellites"]=round(InstablePercentage,4)
#         fileRes.write(
#             "CaseID: " + str(caseid) + "\n"
#                                        "TotalNumberOfMicrosatellites: " + str(totalNum) + "\n"
#                                                                                           "NumberOfInstableMicrosatellites: " + str(
#                 InstableNum) + "\n"
#                                "PercentageOfInstableMicrosatellite: " + str(round(InstablePercentage, 2)) + " %\n"
#         )
#         fileRes.close()
#         lociRes.to_csv(arguments["workspace"] + "detailInfo/" + str(caseid) + "/" + str(caseid) + ".pro")
#         return True
#     except:
#         print("[INFO] procressing " + str(caseid) + " ERROR!!")
#         return False
#
#
# def resultOut():
#     global caseInfo, arguments
#     try:
#         for caseid in caseInfo.index:
#             caseInfo.loc[caseid, "MSIScore(%)"] = \
#                 open(arguments["workspace"] + "detailInfo/" + str(caseid) + "/" + str(caseid)).readlines()[-1].split(
#                     ": ")[
#                     -1][:-2]
#         caseInfo.to_csv(arguments["workspace"] + "Result.csv")
#         return True
#     except:
#         print("[ERROR] error when processing resultOut()")
#         return False


# def main():
#     global caseInfo, Result, arguments
#     if not argumentProcress():
#         return 1
#     if not loadMicroSatellite():
#         return 2
#     if not procressInputConfigure():
#         return 3
#     pool = Pool(arguments["threads"])
#     pool.map(procressOneBam, list(caseInfo.index))
#     pool.close()
#     if not resultOut():
#         return 4
#     return 0
#
def eval(parase):
    print()
# def eval(parase):
#     arg = arguments()
#     if arg:
#         print(arg.parse_args())
#         parase = arg.parse_args()
#         bamPath = parase.bam_file[0]
#         print(bamPath)
#         bedPath = parase.microsatellites[0]
#         ref = parase.reference[0]
#         output = parase.output[0]
#         microsate = parase.microsatellites[0]
#         microsateDF = loadMicroSatellite(microsate)
#         if not microsateDF:
#             return
#
#         bamEvalue(bamPath, bedpath=bedPath, refpath=ref,
#                   outpath=output)
#

def initEval(parase):
    print()


def initSimSeq(parase):
    pa = 0.0001
    pc = 0.0001
    pg = 0.0001
    pt = 0.0001
    qa = 0.00001
    qc = 0.00001
    qg = 0.00001
    qt = 0.00001
    proP = {"A": pa, "C": pc, "G": pg, "T": pt}
    proQ = {"A": qa, "C": qc, "G": qg, "T": qt}
    clusternum = 100000
    return {"p":proP, "q":proQ, "cluster":clusternum}
def getReverseComplematary(read):
    reverseComplematary=""
    for ipos in read[::-1]:
        if ipos=="A":
            reverseComplematary += "T"
        elif ipos =="T":
            reverseComplematary += "A"
        elif ipos =="G":
            reverseComplematary += "C"
        elif ipos =="C":
            reverseComplematary += "G"
    return reverseComplematary

def sequencing(read,readLen=150,para=""):
    curent=np.zeros(para["cluster"])
    # print(curent)
    query=""
    qualList=[]
    for pos in range(readLen):
        posStr=[]
        for ccc in range(len(curent)):
            cc=curent[ccc]
            posStr.append(read[int(cc)])
            pdel = para["p"][posStr[ccc]]
            pins = para["q"][posStr[ccc]]
            pnorm = 1 - pdel - pins
            psteps = pd.array([pdel, pnorm, pins])
            step = np.random.choice([0, 1, 2], p=psteps.ravel())
            curent[ccc] += step

        se=pd.Series(posStr)
        see= se.value_counts(normalize=True)
        index=see.index
        pro=see.iloc[:]
        base=np.random.choice(index, p=pro.ravel())
        print(pro)
        qualityP=np.max(pro)
        query+=base
        qualList.append(qualityP)

        # for ccc in range(len(curent)):
        #     pdel=para["p"][posStr[ccc]]
        #     pins=para["q"][posStr[ccc]]
        #     pnorm=1-pdel-pins
        #     psteps=pd.array([pdel,pnorm,pins])
        #     step=np.random.choice([0,1,2],p=psteps.ravel())

    print(qualList)
    def _pro2basequality(qual):
        if qual>0.99999:
            qual=0.99999
        return chr(int(-10*np.log10(1-qual))+33)

    quals="".join(list(map(_pro2basequality,qualList)))
    return query, quals





        # print(quality)

        # print(index,pro)


        # print(read[i])


    #
    # print(read)
    # print(para)



def buildLib(fa, depth=30, r1Len=150, r2Len=150, insersizeMean=500, inserSizeStd=100,paras="",fq1="",fq2=""):
    fafile = pysam.FastaFile(fa)
    contigs = fafile.references
    # fafile.close()

    def _checkstr(read):
        readChr = set(read)
        chrNum = 0
        for chr in ["A", "T", "C", "G"]:
            if chr in readChr:
                chrNum += 1
        if chrNum < len(readChr):
            return False
        else:
            return True

        # print(read)

    def _processFragments(contig, fragment,fafile):
        if fragment[2]==1:
            leftstart=fragment[0]
            # leftend=fragment[0]+r1Len
            # rightstart=fragment[1]-r2Len
            rightend=fragment[1]
            read1=fafile.fetch(contig,leftstart,leftstart+2*r1Len).upper()
            read2=getReverseComplematary(fafile.fetch(contig,rightend-2*r2Len,rightend).upper())
        else:
            leftstart = fragment[0]
            # leftend = fragment[0] + r2Len
            # rightstart = fragment[1] - r1Len
            rightend = fragment[1]
            read2 = fafile.fetch(contig, leftstart, leftstart+2*r2Len).upper()
            read1 = getReverseComplematary(fafile.fetch(contig, rightend-2*r1Len, rightend).upper())
        if _checkstr(read1) and _checkstr(read2):
            query1,qual1=sequencing(read1,r1Len,paras)
            query2,qual2=sequencing(read2,r2Len,paras)
            readname=uuid.uuid1()
            fq1.write("@"+str(readname)+"\n"+query1+"\n+"+"\n"+qual1+"\n")
            fq2.write("@"+str(readname)+"\n"+query2+"\n+"+"\n"+qual2+"\n")


    def _processOneContig(contig):
        fafile = pysam.FastaFile(fa)
        contigLen = fafile.get_reference_length(contig)
        fragmentsNum = int(contigLen * depth / (r1Len + r2Len))
        print("fra", fragmentsNum)
        fragments = []
        fragmentsCureentNum = 0
        while True:
            insertSize = random.normalvariate(insersizeMean, inserSizeStd)
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
                fragments.append([left, right, direction])
            # if fragmentsCureentNum==fragmentsNum:

            if fragmentsCureentNum == 2000:
                break
        return fragments

    totalfragments = {}
    for contig in contigs:
        print("Process Contig" + contig)
        totalfragments[contig] = _processOneContig(contig)
    fragmentNum=0
    for contig in totalfragments:
        for onefragments in totalfragments[contig]:
            fragmentNum +=1
            print(fragmentNum)
            _processFragments(contig, onefragments,fafile)

    # falen = fafile.get_reference_length(contigs[0])


def sim(parase):
    seqparas = initSimSeq(parase)
    reffa = parase.fasta[0]
    print(parase.fastq)
    fq1=open(parase.fastq[0]+"_R1.fq","w")
    fq2=open(parase.fastq[0]+"_R2.fq","w")
    buildLib(reffa,paras=seqparas,fq1=fq1,fq2=fq2)
    fq1.close()
    fq2.close()

def main():
    arg = arguments()
    if arg:
        parase = arg.parse_args()
        if parase.command == "sim":
            sim(parase)
        else:
            eval(parase)



main()

