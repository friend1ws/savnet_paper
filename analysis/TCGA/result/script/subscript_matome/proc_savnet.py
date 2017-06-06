#! /usr/bin/env python

import sys

input_file = sys.argv[1]
rescue_file = sys.argv[2]
output_file = sys.argv[3]

sptype = {"exon-skip": "Exon skipping",
          "intron-retention": "Intron retention",
          "opposite-side-intron-retention": "Intron retention",
          "alternative-5'-splice-site": "Alternative 5'SS",
          "alternative-3'-splice-site": "Alternative 3'SS",
          "intronic-alternative-5'-splice-site": "Alternative 5'SS",
          "intronic-alternative-3'-splice-site": "Alternative 3'SS"}


target_key = ["Cancer_Type", "Gene_Symbol", "Sample_Name", "Mutation_Key", "Mutation_Type", "Is_Canonical",
              "Splicing_Key", "Splicing_Class", "Supporting_Read_Num", "Score"]


def read_print_savnet_file(infile, hout):

    with open(infile, 'r') as hin:
        header2ind = {}
        header = hin.readline().rstrip('\n').split('\t')
        for (i, cname) in enumerate(header):
            header2ind[cname] = i
        
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[header2ind["IR_filtered"]] == "TRUE": continue
            
            F[header2ind["Splicing_Class"]] = sptype[F[header2ind["Splicing_Class"]]]
            
            F[header2ind["Mutation_Type"]] = F[header2ind["Mutation_Type"]].replace("splicing ", '')
            F[header2ind["Mutation_Type"]] = F[header2ind["Mutation_Type"]][0].upper() + F[header2ind["Mutation_Type"]][1:]
            
            F[header2ind["Is_Canonical"]] = F[header2ind["Is_Canonical"]][0].upper() + F[header2ind["Is_Canonical"]][1:]
            
            F[header2ind["Mutation_Key"]] = "chr" + F[header2ind["Mutation_Key"]]
            F[header2ind["Splicing_Key"]] = "chr" + F[header2ind["Splicing_Key"]]
            
            print >> hout, '\t'.join([F[header2ind[x]] for x in target_key])



hout = open(output_file, 'w')    
print >> hout, '\t'.join(["Cancer type", "Gene", "Sample name", "Mutation key", "Mutation type",
                          "Canonical or not", "Splicing position", "Splicing outcome", "Supporting read count", "Bayes factor"])
 
read_print_savnet_file(input_file, hout)
read_print_savnet_file(rescue_file, hout)

hout.close()



