import sys
from collections import OrderedDict

##this script identify if the cds sequence of any gene has multiple stop codons, note that the script assume that sequences are on the forward strand.
##if you do not have the sequence on the forward strand, you can use ReverseComplementSeq.py to revese-complement your script and then use this script.
##usage python3 IdentifyCdsMultiStopCod.py <input_fasta_file> <output_file>
##the output file will contain the header of the sequences that are not divisible by three or that have multiple stop codons. 

def ReplaceStopCodon(fil,file2):
    SeqDict=MakeSeqDict(fil)
    WriteId(SeqDict,file2)

def MakeSeqDict(fil):
    SeqDict=OrderedDict()
    Seq=''
    with open(fil) as fasta:
        for line in fasta:
            if line.startswith(">"):
                if len(Seq)!=0:
                    SeqDict[Header]=Seq
                    Seq=''
                Header=line.strip()
            else:
                line=line.strip()
                Seq+=line
    return SeqDict

def WriteId(SeqDict,file2):
    dest=open(file2,'w')
    Stop=["TAG","TAA","TGA"]
    Keys=list(SeqDict.keys())
    for key in Keys:
        Seq=SeqDict[key]
        StopCodon=0
        Leng=len(Seq)%3
        if Leng!=0:
            dest.write(key)
            dest.write("\n")
        else:
            Triplets=[Seq[i:i+3] for i in range(0, len(Seq), 3)]
            for i in Triplets:
                if i in Stop:
                    StopCodon+=1
            if StopCodon>1:
                dest.write(key)
                dest.write("\n")
    dest.close()
inp=sys.argv[1]
inp1=sys.argv[2]
ReplaceStopCodon(inp,inp1)
