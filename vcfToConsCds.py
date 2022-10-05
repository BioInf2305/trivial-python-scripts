import sys
import argparse
import re
import random
from pyfasta import Fasta
from pysam import VariantFile
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def writeFasta(vcfIn, gffIn, fastaIn, transcriptFIn, sampleMap, sampleMode, popMode):

    sampleToPopDict = sampleMapToPopDict(sampleMap)
    transcriptList = transcriptFtoList(transcriptFIn)
    transcriptCdsDict, transcriptStrandDict = transcriptListToDict(transcriptList, gffIn)
    refFastaDict = fastaInToDict(fastaIn, transcriptCdsDict)
    sampleAltFastaDict = vcfToCons(vcfIn, refFastaDict, sampleMode)
    allSampleConsDict = {}
    ambiNuclDict = {"R":"AG","Y":"CT","S":"GC","W":"AT","K":"GT","M":"AC"}
    for transcript in sampleAltFastaDict:
        allSampleConsDict = {}
        samples = sampleAltFastaDict[transcript]
        with open(transcript+".indi.consensus.txt","w") as dest:
            for sample in samples:
                dest.write(">"+sample+"\n")
                sampleSeq = Seq("".join(list(samples[sample].values())),IUPAC.ambiguous_dna)
                if transcriptStrandDict[transcript] == "-":
                    sampleSeq = sampleSeq.reverse_complement()
                sampleSeq = str(sampleSeq)
                sampleSeqSubList = [sampleSeq[i:i+60] for i in range(0,len(sampleSeq),60)]
                for seq in sampleSeqSubList:
                    dest.write("".join(seq))
                    dest.write("\n")
                for position in samples[sample]:
                    if sampleToPopDict[sample] not in allSampleConsDict:
                        allSampleConsDict[sampleToPopDict[sample]] = {}
                    if position not in allSampleConsDict[sampleToPopDict[sample]]:
                        allSampleConsDict[sampleToPopDict[sample]][position] = []
                    allSampleConsDict[sampleToPopDict[sample]][position].append(samples[sample][position])
            with open(transcript+".pop.consensus.txt","w") as dest1:
                for pop in allSampleConsDict:
                    dest1.write(">"+pop+"\n")
                    popSeq=[]
                    popConsDict = allSampleConsDict[pop]
                    for position in popConsDict:
                        uniqAlleles = list(set(popConsDict[position]))
                        if len(uniqAlleles)>3:
                            print("ERROR: multi-allelic position "+position)
                            sys.exit(1)
                        elif len(uniqAlleles) > 1 and len(uniqAlleles) <= 3:
                            alleleCount = [popConsDict[position].count(allele) for allele in uniqAlleles]
                            base = uniqAlleles[alleleCount.index(max(alleleCount))]
                        else:
                            base = uniqAlleles[0]
                        if base in ambiNuclDict:
                            if popMode == "4":
                                randIdx = random.randint(0,1)
                                base = ambiNuclDict[base][randIdx]
                        popSeq.append(base)
                    if transcriptStrandDict[transcript] == "-":
                        popSeq = str(Seq("".join(popSeq),IUPAC.ambiguous_dna).reverse_complement())
                    popSeqSubList = [popSeq[i:i+60] for i in range(0,len(popSeq),60)]
                    for seq in popSeqSubList:
                        dest1.write("".join(seq))
                        dest1.write("\n")

def sampleMapToPopDict(sampleMap):
    sampleToPopDict = {}
    with open(sampleMap) as source:
        for line in source:
            line = line.rstrip()
            line = line.split()
            sampleToPopDict[line[0]] = line[1]
    return sampleToPopDict

def transcriptFtoList(transcriptFIn):

    transcriptList = []
    with open(transcriptFIn) as source:
        for line in source:
            line = line.rstrip().split()
            transcriptList.append(line[0])
    return transcriptList


def transcriptListToDict(transcriptList, gffIn):

    transcriptCdsDict = {}
    transcriptStrandDict = {}
    with open(gffIn) as source:
        for lineC in source:
            line = lineC.rstrip().split("\t")
            if not lineC.startswith("#"):
                if line[2] == "CDS":
                    pattern = re.compile('Parent=transcript:([^;]+)')
                    match = re.findall(pattern, line[8])
                    if len(match) == 0 or len(match)>1:
                        print("the cds id is not captured by regex")
                        sys.exit(1)
                    else:
                        if(match[0] in transcriptList):
                            cdsList = [line[0],int(line[3]),int(line[4])]
                            if match[0] not in transcriptCdsDict:
                                transcriptCdsDict[match[0]] = []
                            transcriptCdsDict[match[0]].append(cdsList)
                            transcriptStrandDict[match[0]] = line[6]
    return transcriptCdsDict, transcriptStrandDict

def fastaInToDict(fastaIn, transcriptCdsDict):

    refFasta = Fasta(fastaIn)
    sampleFastaDict = {}
    for transcript in transcriptCdsDict:
        sampleFastaDict[transcript] = {}
        cdsList = transcriptCdsDict[transcript]
        for cds in cdsList:
            if cds[0] not in sampleFastaDict[transcript]:
                sampleFastaDict[transcript][cds[0]] = {}
            tmpSeq = refFasta[cds[0]][cds[1]-1:cds[2]]
            genomePos = cds[1]
            for base in tmpSeq:
                sampleFastaDict[transcript][cds[0]][genomePos]=base
                genomePos +=1
    return sampleFastaDict

def vcfToCons(vcfIn, refFastaDict, mode):

    vcfF = VariantFile(vcfIn)
    sampleAltConsDict = {}
    totalSamples = list((vcfF.header.samples))
    ambiCodes = {"AT":"W","TA":"W","AC":"M","CA":"M","AG":"R","GA":"R","TC":"Y","CT":"Y","TG":"K","GT":"K","GC":"S","CG":"S","AA":"A","TT":"T","CC":"C","GG":"G"}
    for transcript in refFastaDict:
        tmpDict = refFastaDict[transcript]
        chrom = list(tmpDict.keys())[0]
        startPos = list(tmpDict[chrom].keys())[0]
        endPos = list(tmpDict[chrom].keys())[-1]
        sampleAltConsDict[transcript] = {}
        for sample in totalSamples:
            sampleAltConsDict[transcript][sample] = refFastaDict[transcript][chrom].copy()
        for rec in vcfF.fetch(chrom,startPos, endPos):
            if rec.pos in refFastaDict[transcript][chrom]:
                if len(rec.alts)>1 or len(rec.ref)>1:
                    print("WARNING: MULTIALLELIC VARIANT, SKIPPING "+ str(rec.chrom)+" "+str(rec.pos)+"\n")
                elif len(rec.alts[0])>1 or len(rec.ref[0])>1:
                    print("WARNING: MULTIALLELIC VARIANT, SKIPPING "+ str(rec.chrom)+" "+str(rec.pos)+"\n")
                elif (rec.ref[0]).lower() != refFastaDict[transcript][chrom][rec.pos].lower():
                    print("ERROR at position "+ str(chrom)+ " "+str(rec.pos)+" reference allele in fasta and in vcf does not match")
                    sys.exit(1)
                else:
                    for sample in totalSamples:
                        gt = rec.samples[sample]["GT"]
                        if gt == (0,1) or gt == (1,1):
                            if gt == (0,1):
                                geno = rec.ref[0] + rec.alts[0]
                            elif gt == (1,1):
                                geno = rec.alts[0] + rec.alts[0]
                            else:pass
                            if mode == "1":
                                base = geno[1]
                            elif mode == "2":
                                base = ambiCodes[geno]
                            elif mode == "3":
                                randNum = random.randint(0,1)
                                print(geno)
                                print("here is the random number ", randNum)
                                base = geno[randNum]
                            else:
                                print("ERROR: the mode should be either 1, 2 or 3, for more info run help argument")
                                sys.exit(1)
                            sampleAltConsDict[transcript][sample][rec.pos] = base
    return sampleAltConsDict

if __name__ =="__main__":
    parser = argparse.ArgumentParser(description = "\
    This python script will make consensus for each sample, \
    each population in vcf file given gff file, reference sequence and file of transcript ids",epilog = "author: Maulik Upadhyay (Upadhyay.maulik@gmail.com)")
    parser.add_argument('-v',"--varFile",metavar = "File", help ="compressed and indexed vcf or bcf file",required = True)
    parser.add_argument('-g',"--gffFile",metavar = "File", help ="gff file",required = True)
    parser.add_argument('-f',"--fasFile",metavar = "File", help = "reference sequence in fasta format (used to generated vcf file)", required = True)
    parser.add_argument('-t',"--transcriptFile",metavar = "File",help = "file containing transcript ids in each new line; note that transcript should \
    be present in the gff file", required = True)
    parser.add_argument('-s',"--sampleMapFile",metavar = "File",help = "file containing two co: sampleid(first column) and popId (second column)",required = True)
    parser.add_argument('-m',"--sampleMode",metavar = "Str", help ="rule to create consesnsus, \
    there are four modes:,\
    1:replace 1/0 or 1/1 in sample with only alternative allele\
    2: replace 1/1 with alternative allele and 1/0 with IUPAC ambiguous base\
    3: in case of 1/1 or 1/0 select base randomly, default is 1",default = "1", required = False)
    parser.add_argument('-M',"--popMode",metavar = "String",help="rule to create population consensus, most frequent base will be considered,\
    if the value is any other than 4, ambiguous bases will be kept as is",default = "4", required = False)

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    else:
        writeFasta(args.varFile, args.gffFile, args.fasFile, args.transcriptFile, args.sampleMapFile, args.sampleMode, args.popMode)
