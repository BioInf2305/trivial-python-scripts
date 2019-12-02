import sys
from collections import OrderedDict
from Bio.Seq import Seq

##usage: the script can be used to reverse complement the sequences given its header 
##for example python3 ReverseComplementSeq.py Example1.fasta SequenceHeaders.txt

def RevCompFasta(Fasta,GeneInfo):
    dest=open(Fasta+".strandAdjusted","w")
    dest1=open(GeneInfo+".missingGeneInFastq","w")
    with open(GeneInfo) as InfoFile:
        InfoDict=OrderedDict()
        GeneList=[]
        for line in InfoFile:
            a=line.split()
            InfoDict[a[1]]=[a[0],a[2]]
    with open(Fasta) as FastaFile:
        ReverseComplement=False
        for line in FastaFile:
            if line.startswith(">"):
                line=line.strip()
                Ke=line.strip(">")
                if InfoDict[Ke][1]=="-":
                    ReverseComplement=True
                else:ReverseComplement=False
                GeneList.append(Ke)
                dest.write(line)
                dest.write("\n")
            elif ReverseComplement==True:
                print("reverse complementing the transcript ", Ke)
                dest.write(str(Seq(line.rstrip()).reverse_complement()))
                dest.write("\n")
            else:
                dest.write(str(line.rstrip()))
                dest.write("\n")
    GffGeneList=list(InfoDict.keys())
    for gene in GffGeneList:
        if gene not in GeneList:
            dest1.write(gene)
            dest1.write("\n")
    dest.close()
inp=sys.argv[1]
inp1=sys.argv[2]
RevCompFasta(inp,inp1)
