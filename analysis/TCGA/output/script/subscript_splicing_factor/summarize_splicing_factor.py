#! /usr/bin/env python

import sys, os, glob, re

# savnet_result = sys.argv[1]

key2ctype = {}

print "Cancer_Type" + '\t' + "SF_gene" + '\t' + "Splicing_Key" + '\t' + "Gene_Symbol" + '\t' + "Splicing_Class" + '\t' + "Minus_Log_PV"

SJ_black_list = glob.glob("../splicing_factor/*_SJ.sig.txt")
for bfile in sorted(SJ_black_list):
    ctype, sf_gene, sp_type = os.path.basename(bfile).replace(".sig.txt", "").split("_")
    with open(bfile, 'r') as hin:

        header2ind = {}
        header = hin.readline().rstrip('\n').split('\t')
        for (i, cname) in enumerate(header):
            header2ind[cname] = i

        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = F[header2ind["SJ_1"]].strip('"') + ':' + F[header2ind["SJ_2"]].strip('"') + '-' + F[header2ind["SJ_3"]].strip('"')

            #####
            # take gene name 
            genes = F[header2ind["Gene_1"]].split(';') + F[header2ind["Gene_2"]].split(';')
            genes = sorted(list(set(genes)))

            if "---" in genes: genes.remove("---")
            if len(genes) > 0:
                genes_nm = filter(lambda x: x.find("(NM_") > 0, genes)
                if len(genes_nm) > 0: genes = genes_nm

            if len(genes) > 0:
                genes_single = filter(lambda x: x.find("-") == -1, genes)
                if len(genes_single) > 0: genes = genes_single
 
            gene = genes[0]
            gene = re.sub(r"\(N[MR]_\d+\)", "", gene)
   
            F[header2ind["Splicing_Class"]] = re.sub(r'^intronic-', "", F[header2ind["Splicing_Class"]]) 
            print ctype + '\t' + sf_gene + '\t' + key + '\t' + gene + '\t' + F[header2ind["Splicing_Class"]] +  '\t' + str(round(float(F[header2ind["Minus_Log_PV"]]), 3)) # + '\t' + str(round(float(F[header2ind["Effect_Size"]]), 3))


IR_black_list = glob.glob("../splicing_factor/*_IR.sig.txt")
for bfile in sorted(IR_black_list):
    ctype, sf_gene, sp_type = os.path.basename(bfile).replace(".sig.txt", "").split("_")
    with open(bfile, 'r') as hin:

        header2ind = {}
        header = hin.readline().rstrip('\n').split('\t')
        for (i, cname) in enumerate(header):
            header2ind[cname] = i

        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = F[header2ind["Chr"]].strip('"') + ':' + F[header2ind["Boundary_Pos"]].strip('"') + '-' + F[header2ind["Boundary_Pos"]].strip('"')
        
            print ctype + '\t' + sf_gene + '\t' + key + '\t' + F[header2ind["Gene_Symbol"]] + '\t' + "intron-retention" + '\t' + str(round(float(F[header2ind["Minus_Log_PV"]]), 3)) # + '\t' + str(round(float(F[header2ind["Effect_Size"]]), 3))
 

