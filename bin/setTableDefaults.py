#!/usr/bin/env python

import sys
import string

flankDict = {
    "tss":2000,
    "tes":2000,
    "genebody":2000,
    "exon":500,
    "cgi":500,
    "enhancer":1500,
    "dhs":1000,
}

pointLabDict = {
    "tss":"TSS",
    "tes":"TES",
    "genebody":"TSS,TES",
    "exon":"Acceptor,Donor",
    "cgi":"Left,Right",
    "enhancer":"Enhancer",
    "dhs":"Left,Right",
}

# FItype = {
#     "tss":"endsite",
#     "tes":"endsite",
#     "genebody":"datatype",
#     "exon":"subregion",
#     "cgi":"subregion",
#     "enhancer":"subregion",
#     "dhs":"subregion"
# }

in_f = file(sys.argv[1])
out_f = file(sys.argv[2], "w")
out_f.write("\t".join(["Genome", "DefaultDB", "Region", "DefaultFI1", \
                       "DefaultFI2", "DefaultFI3", "PointLab", "Flank"]) + \
            "\n")

# while(True):
for line in in_f:
    # line = in_f.readline()
    # if len(line) == 0:
    #     break;
    # line = string.strip(line)

    lineL = line.rstrip().split("\t")
    genome = lineL[0]
    defaultDB = lineL[1]
    region = lineL[2]

    if region == "cgi":
        fi_1 = "NA"
        fi_2 = "ProximalPromoter"
        fi_3 = "protein_coding"
    elif region == "dhs":
        fi_1 = "H1hesc"
        fi_2 = "ProximalPromoter"
        fi_3 = "protein_coding"
    elif region == "exon":
        fi_1 = "chipseq"
        fi_2 = "canonical"
        fi_3 = "protein_coding"
    elif region == "genebody":
        fi_1 = "chipseq"
        fi_2 = "NA"
        fi_3 = "protein_coding"
    elif (region == "enhancer") and (genome == "hg19"):
        fi_1 = "H1hesc"
        fi_2 = "genebody"
        fi_3 = "protein_coding"
    elif (region == "enhancer") and (genome == "mm9"):
        fi_1 = "mESC"
        fi_2 = "genebody"
        fi_3 = "protein_coding"
    out_f.write("\t".join([genome, defaultDB, region, fi_1, fi_2, fi_3, \
                pointLabDict[region], str(flankDict[region])]) + "\n")
    
    # Extra: TSS and TES if region = genebody.
    if region == "genebody":
        region = "tss"
        out_f.write("\t".join([genome, defaultDB, region, fi_1, fi_2, fi_3, \
                    pointLabDict[region], str(flankDict[region])]) + "\n")
        region = "tes"
        out_f.write("\t".join([genome, defaultDB, region, fi_1, fi_2, fi_3, \
                    pointLabDict[region], str(flankDict[region])]) + "\n")

in_f.close()
out_f.close()
