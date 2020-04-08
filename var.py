import pysam
import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# fafile = ""
global thisContigVarDf, startEndDf, seqStat, regionStat


def initVarArg(parase):
    """
    Initialize the parameters for reads simulations
    :param parase:
    :return: dict of simulations parameters
    """
    global paras
    paras = {}
    paras["input"] = parase.input[0]
    paras["output"] = parase.output[0]
    if len(parase.var_conf[0]) > 0:
        paras["var_conf"] = parase.var_conf[0]
    else:
        paras["var_conf"] = False
    if len(parase.vcf[0]) > 0:
        paras["vcf"] = parase.vcf[0]
    else:
        paras["vcf"] = False
    if len(parase.conf[0]) > 0:
        paras["conf"] = parase.conf[0]
    else:
        paras["conf"] = False

    return paras


def arrangeVar2Contig():
    # global startEndDf
    varConfPerContig = {}
    yamlLoader = open(paras["var_conf"], 'r', encoding='utf-8')
    varconf = yamlLoader.read()
    varconfdict = yaml.load(varconf)
    outfa = open(paras["output"], "w")
    fafile = pysam.FastaFile(paras["input"])
    contigs = fafile.references
    contigsLengths = fafile.lengths
    # fafile.close()
    lenPerContig = {}
    totalLen = 0
    for contig, contigLen in zip(contigs, contigsLengths):
        lenPerContig[contig] = contigLen
        totalLen += contigLen
    for contig in lenPerContig:
        varConfPerContig[contig] = {}
        for vartype in varconfdict:
            varConfPerContig[contig][vartype] = {}
            for item in varconfdict[vartype]:
                if item == "num":
                    varConfPerContig[contig][vartype][item] = \
                        lenPerContig[contig] / (totalLen) * varconfdict[vartype][item]
                else:
                    varConfPerContig[contig][vartype][item] = varconfdict[vartype][item]

    def _getDistribution(minLen, maxLen, num, disType="pwoerlow"):
        if disType == "uniform":
            lens = np.random.uniform(minLen, maxLen, num).round()
            lens.sort()
            return lens[::-1]
        if disType == "pwoerlow":
            lens = np.random.exponential(scale=(maxLen - minLen) / 10, size=num).round() + minLen - 1
            lens.sort()
            return lens[::-1]

    def _getOneContigLenDis(varConf):
        thisContigVarDf = pd.DataFrame(columns=["chr", "genomeSpan", "varLen", "type"])
        thisContigVarDf.index.name = "id"

        for vartype in varConf:
            if vartype == "DEL" or vartype == "INV" or vartype == "TD":
                varLen = _getDistribution(minLen=varConf[vartype]["minLen"],
                                          maxLen=varConf[vartype]["maxLen"],
                                          num=int(varConf[vartype]["num"]),
                                          )
                thisTypeVarNum = len(varLen)
                thisContigThisTypeVarDf = pd.DataFrame({"genomeSpan": varLen,
                                                        "varLen": varLen},
                                                       index=["_".join([contig, vartype, str(i)]) for i in
                                                              range(thisTypeVarNum)])
                thisContigThisTypeVarDf["chr"] = contig
                thisContigThisTypeVarDf["type"] = vartype
                thisContigVarDf = pd.concat([thisContigVarDf, thisContigThisTypeVarDf])



            elif vartype == "INS" or vartype == "SNP":
                varLen = np.ones(int(varConf[vartype]["num"]))
                thisTypeVarNum = len(varLen)
                thisContigThisTypeVarDf = pd.DataFrame({"varLen": varLen},
                                                       index=["_".join([contig, vartype, str(i)]) for i in
                                                              range(thisTypeVarNum)])
                thisContigThisTypeVarDf["genomeSpan"] = 1
                thisContigThisTypeVarDf["chr"] = contig
                thisContigThisTypeVarDf["type"] = vartype
                thisContigVarDf = pd.concat([thisContigVarDf, thisContigThisTypeVarDf])

            # varLenDisContig[vartype] = varLen
        return thisContigVarDf

    def arrangePos4Vars(contig, contigLen):
        global startEndDf
        markedlen = 10 ** len(str(int(thisContigVarDf["genomeSpan"].max()))) + 1
        searchRegion = np.ones(contigLen)
        searchSet = np.arange(contigLen)
        for s, e in zip(start, end):
            if e >= contigLen - markedlen:
                searchRegion[s:] = 0
                break
            searchRegion[s:e + markedlen] = 0

        for id, info in thisContigVarDf.iterrows():
            # print(thisContigVarDf)
            if info["genomeSpan"] * 10 < markedlen:
                markedlen = markedlen // 10 + 1
                searchRegion = np.ones(contigLen)
                searchSet = np.arange(contigLen)
                for s, e in zip(start, end):
                    if e >= contigLen - markedlen:
                        searchRegion[e:] = 0
                        break
                    searchRegion[s:e + markedlen] = 0
            thisend = int(np.random.choice(searchSet, p=searchRegion / searchRegion.sum()))
            thisstart = int(thisend - info["genomeSpan"])
            start.append(thisstart)
            end.append(thisend)
            searchRegion[thisstart:thisend + markedlen] = 0
            # print(pd.Series({"start":thisstart,"end":thisend},name=id))
            startEndDf = startEndDf.append(pd.Series({"start": thisstart, "end": thisend}, name=id))
            # print(startEndDf)

            # print(id, info, thisend)
        return

    def _write2fa(contig):
        global thisContigVarDf
        print("write2fa")
        snpDict = {"A": [["T", "C", "G"], [0.4, 0.3, 0.3]], "T": [["A", "C", "G"], [0.4, 0.3, 0.3]],
                   "C": [["T", "A", "G"], [0.3, 0.3, 0.4]], "G": [["T", "C", "A"], [0.3, 0.4, 0.3]]}

        seq = list(fafile.fetch(contig))
        thisContigVarDf["start"] = thisContigVarDf["start"].apply(int)
        thisContigVarDf["end"] = thisContigVarDf["end"].apply(int)
        thisContigVarDf["varLen"] = thisContigVarDf["varLen"].apply(int)
        thisContigVarDf["genomeSpan"] = thisContigVarDf["genomeSpan"].apply(int)
        print(thisContigVarDf)

        def _getReverseComplematary(read):
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

        def _randomSeq(seqLen):
            seq = []
            for i in range(seqLen):
                seq.append(np.random.choice(["A", "C", "G", "T"]))
            return "".join(seq)

        for id, info in thisContigVarDf.iterrows():
            thisstart = info["start"]
            thisend = info["end"]
            vartype = info["type"]
            if vartype == "SNP":

                seq[thisstart] = np.random.choice(snpDict[seq[thisstart].upper()][0],
                                                  p=snpDict[seq[thisstart].upper()][1])
            elif vartype == "DEL":
                for delpos in range(thisstart, thisend):
                    seq[delpos] = "&"
            elif vartype == "INS":
                seq[thisstart] == seq[thisstart] + _randomSeq(info["varLen"])
            elif vartype == "DUP":
                seq[thisend] == seq[thisstart:thisend] + seq[thisend]
            elif vartype == "INV":
                seq[thisstart:thisend] = list(_getReverseComplematary("".join(seq[thisstart:thisend])))
        outfa.write(">"+contig+"\n")
        # seqstr="".join(seq)
        # seqstr.replace()
        outfa.write("".join(seq).replace("&","")+"\n")



        # startEndDf=pd.DataFrame()

    # outfa=pysam.fasta()
    for contig, contigLen in zip(varConfPerContig, contigsLengths):
        global thisContigVarDf, seqStat, regionStat, startEndDf
        startEndDf = pd.DataFrame()
        thisContigVarDf = _getOneContigLenDis(varConfPerContig[contig])
        thisContigVarDf.sort_values(inplace=True, ascending=False, by=["genomeSpan"])
        seqStat = np.ones(contigLen)
        pos = -1  # from 1
        start, end = [], []
        lastpos = "A"
        for contigStr in "A" + fafile.fetch(contig):
            if contigStr not in ["A", "T", "C", "G", "a", "t", "c", "g"]:
                seqStat[pos] = 0
                if lastpos in ["A", "T", "C", "G", "a", "t", "c", "g"]:
                    print("start", pos)
                    start.append(pos)
            if contigStr in ["A", "T", "C", "G", "a", "t", "c", "g"]:
                seqStat[pos] = 0
                if lastpos not in ["A", "T", "C", "G", "a", "t", "c", "g"]:
                    print("end", pos)
                    end.append(pos)
            lastpos = contigStr
            pos += 1
        print("arrangePos4Vars", contig, contigLen)
        arrangePos4Vars(contig, contigLen)
        thisContigVarDf = pd.concat([thisContigVarDf, startEndDf], axis=1)
        thisContigVarDf.sort_values(inplace=True, ascending=False, by=["end"])
        _write2fa(contig)
    fafile.close()
    outfa.close()


def var(parase):
    # print(parase)
    initVarArg(parase)
    arrangeVar2Contig()
    print("hhhhhh")
