import sys
import math
import argparse
from collections import OrderedDict
from pysam import VariantFile

def filterSitesAlleleCount(bcfIn, minAlleles, maxAlleles, bcfOut):
    minAlleles=int(minAlleles)
    maxAlleles=int(maxAlleles)
    bcf_in=VariantFile(bcfIn)
    bcf_out = VariantFile(bcfOut, 'w', header=bcf_in.header)
    for rec in bcf_in.fetch():
        if rec.alts!=None:
            if len(rec.alts)>=minAlleles and len(rec.alts)<=maxAlleles and "*" not in rec.alts and len(rec.ref)==1:
                bcf_out.write(rec)

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="A small python script to filter the sites based on \
    user defined total number of alternative alleles identified ", epilog="author: Maulik Upadhyay (Upadhyay.maulik@gmail.com)")
    parser.add_argument("-v","--vcfF",metavar="File",help="compressed or uncompressed vcf or \
    bcf file",required=True)
    parser.add_argument("-m","--minAlt",metavar="Int",help="minimum number of alternative allle counts required to be \
    present at each record. OPTIONAL FLAG",default=1,required=False)
    parser.add_argument("-M","--maxAlt",metavar="Int",help="max. number of alternative allele counts required to be \
    present at each record. OPTIONAL FLAG", default=1, required=False)
    parser.add_argument("-o","--outputF",metavar="Str",help="bcf or vcf file name for the output", required=True)
    args = parser.parse_args()
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    else:
        filterSitesAlleleCount(args.vcfF, args.minAlt, args.maxAlt, args.outputF)
