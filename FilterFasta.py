import sys


def FilterFasta(inFile,outFile,l):
    dest=open(outFile,"w")
    seq=""
    with open(inFile) as source:
        for line in source:
            if line.startswith(">"):
                if seq!="":
                    if ("*" not in seq) and len(seq)>=int(l):
                        dest.write(header+"\n")
                        tmpList=[seq[i:i+60] for i in range(0,len(seq),60)]
                        dest.write("\n".join(tmpList)+"\n")
                        del tmpList[:]
                seq=""
                header=line.strip()
            else:
                seq+=line.strip()
    if ("*" not in seq) and len(seq)>=int(l):
        dest.write(header+"\n")
        tmpList=[seq[i:i+60] for i in range(0,len(seq),60)]
        dest.write("\n".join(tmpList)+"\n")

if __name__=="__main__":
    FilterFasta(sys.argv[1],sys.argv[2],sys.argv[3])
