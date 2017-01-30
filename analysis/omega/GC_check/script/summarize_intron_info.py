#! /usr/bin/env python

import sys, subprocess
import pysam

from annot_utils.coding import *


gsm_file = sys.argv[1]
output_file = sys.argv[2]
reference = sys.argv[3]

def get_seq(reference, pos):

    seq = ""    
    for item in pysam.faidx(reference, pos):
        seq = seq + item.rstrip('\n')
    seq = seq.replace('>', '')
    seq = seq.replace(pos, '')

    return seq



make_coding_info(output_file + ".refCoding.bed.gz", "refseq", "hg19", True, False)
ref_coding_tb = pysam.TabixFile(output_file + ".refCoding.bed.gz")

mut2info = {}
mut2count = {}
with open(gsm_file, 'r') as hin:

    header2ind = {}
    header = hin.readline().rstrip('\n').split('\t')
    for (i, cname) in enumerate(header):
        header2ind[cname] = i

    for line in hin:
        F = line.rstrip('\t').split('\t')
        # if F[header2ind["Mutation_Type"]] not in ["splicing donor disruption", "splicing acceptor disruption"]: continue
        # mut = F[header2ind["Mutation_Key"]]

        mut = ','.join([F[header2ind[x]] for x in ["Sample_Name", "Chr_Mut", "Start_Mut", "End_Mut", "Ref_Mut", "Alt_Mut"]])
        gsms = F[header2ind["GenomonSplicingMutation"]].split(';')
        if gsms[0] == "---":
            gsm = "no-change"
        elif len(gsms) == 2 and "retention" in gsms[0] and "retention" in gsms[1]:
            gsm = "intron-retention"
        elif len(gsms) >= 2:
            gsm = "complex"
        else:
            gsm = gsms[0]

        # if mut not in mut2count:
            # mut2count[mut] = int(F[header2ind["Supporting_Read_Num"]])
        mut2info[mut] = F[header2ind["Type_Motif"]] + '\t' + F[header2ind["Gene_Symbol"]] + '\t' + gsm
        # elif int(F[header2ind["Supporting_Read_Num"]]) >= mut2count[mut]:
        #     mut2count[mut] = int(F[header2ind["Supporting_Read_Num"]])
        #     mut2info[mut] = F[header2ind["Type_Motif"]] + '\t' + F[header2ind["Gene_Symbol"]] + '\t' + gsm


hout = open(output_file, 'w')
print >> hout, "Mutation_Key" + '\t' + "Type_Motif" + '\t' + "Splice_Class" + '\t' + "Gene_Symbol" + '\t' + "Strand" + '\t' + \
               "Len_intron_5prime" + '\t' + "Len_exon" + '\t' + "Len_intron_3prime" + '\t' + \
               "GC_intron_5prime" + '\t' + "GC_exon" + '\t' + "GC_intron_3prime" 


process_n = 0
for mut in mut2info:
    MM = mut.split(',')
    mut_type, gene_symbol, sp_class = mut2info[mut].split('\t')
    gene_strand = ""

    if mut == "17,7578291,T,G":
        pass

    process_n = process_n + 1
    if process_n % 1000 == 0:
        print >> sys.stderr, str(process_n) + " records have been processed!"

    ####################
    # get the records for control junction data for the current position
    tabixErrorFlag = 0
    try:
        records = ref_coding_tb.fetch(MM[1], int(MM[2]) - 10, int(MM[2]) + 10)
    except Exception as inst:
        # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
        tabixErrorMsg = str(inst.args)
        tabixErrorFlag = 1

    introns_direct = []
    exons = []
    if tabixErrorFlag == 0:
        for record_line in records:
            record = record_line.split('\t')
            if gene_symbol != record[3]: continue
            gene_strand = record[5]
            if record[4] == "intron": introns_direct.append((record[0], record[1], record[2]))
            if record[4] in ["noncoding", "coding"]: exons.append((record[0], record[1], record[2]))
            

    introns_direct = list(set(introns_direct))
    exons = list(set(exons))

    # multiple introns or exons
    if len(introns_direct) != 1 or len(exons) != 1: continue


    introns_opposite = []
    # intron -> exon
    if introns_direct[0][2] == exons[0][1]:

        ####################
        # get the records for control junction data for the current position
        tabixErrorFlag2 = 0
        try:
            records = ref_coding_tb.fetch(MM[1], int(exons[0][2]) - 10, int(exons[0][2]) + 10)
        except Exception as inst:
            # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
            tabixErrorMsg = str(inst.args)
            tabixErrorFlag = 1

        if tabixErrorFlag2 == 0:
            for record_line in records:
                record = record_line.split('\t')
                if gene_symbol != record[3]: continue
                if record[4] == "intron": introns_opposite.append((record[0], record[1], record[2]))

    # exon -> intron
    elif introns_direct[0][1] == exons[0][2]:

        ####################
        # get the records for control junction data for the current position
        tabixErrorFlag2 = 0
        try:
            records = ref_coding_tb.fetch(MM[1], int(exons[0][1]) - 10, int(exons[0][1]) + 10)
        except Exception as inst:
            # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
            tabixErrorMsg = str(inst.args)
            tabixErrorFlag = 1
        
        if tabixErrorFlag2 == 0:
            for record_line in records:
                record = record_line.split('\t')
                if gene_symbol != record[3]: continue
                if record[4] == "intron": introns_opposite.append((record[0], record[1], record[2]))


    introns_opposite = list(set(introns_opposite))

    # multiple introns or exons
    if len(introns_opposite) != 1: continue


    tchr, tstart, tend = exons[0][0], exons[0][1], exons[0][2]
    intron_prev_seq = get_seq(reference, tchr + ':' + str(int(tstart) - 150) + '-' + str(int(tstart) - 20))
    exon_seq = get_seq(reference, tchr + ':' + str(int(tstart) + 6) + '-' + str(int(tend) - 5))
    intron_after_seq = get_seq(reference, tchr + ':' + str(int(tend) + 20 + 1) + '-' + str(int(tend) + 150 + 1))

    if len(exon_seq) < 30: continue

    intron_prev_GC_ratio = float(intron_prev_seq.count("G") + intron_prev_seq.count("C")) / len(intron_prev_seq)
    exon_GC_ratio = float(exon_seq.count("G") + exon_seq.count("C")) / len(exon_seq)
    intron_after_GC_ratio = float(intron_after_seq.count("G") + intron_after_seq.count("C")) / len(intron_after_seq)

    if gene_strand == "+":
        GC_intron_5prime = str(round(intron_prev_GC_ratio, 3))
        GC_exon = str(round(exon_GC_ratio, 3))
        GC_intron_3prime = str(round(intron_after_GC_ratio, 3))
    else:
        GC_intron_5prime = str(round(intron_after_GC_ratio, 3))
        GC_exon = str(round(exon_GC_ratio, 3))
        GC_intron_3prime = str(round(intron_prev_GC_ratio, 3))


    # intron -> exon and "+" or exon -> intron and "-"
    if (introns_direct[0][2] == exons[0][1] and gene_strand == "+") or (introns_direct[0][1] == exons[0][2] and gene_strand == "-"):
        len_intron_5prime = str(int(introns_direct[0][2]) - int(introns_direct[0][1]))
        len_exon = str(int(exons[0][2]) - int(exons[0][1]))
        len_introns_3prime = str(int(introns_opposite[0][2]) - int(introns_opposite[0][1]))
    else:
        len_intron_5prime = str(int(introns_opposite[0][2]) - int(introns_opposite[0][1]))
        len_exon = str(int(exons[0][2]) - int(exons[0][1]))
        len_introns_3prime = str(int(introns_direct[0][2]) - int(introns_direct[0][1]))


    print >> hout, mut + '\t' + mut_type + '\t' + sp_class + '\t' + gene_symbol + '\t'+  gene_strand + '\t' + \
                    len_intron_5prime + '\t' + len_exon + '\t' + len_introns_3prime + '\t' + \
                    GC_intron_5prime + '\t' + GC_exon + '\t' + GC_intron_3prime


subprocess.call(["rm", "-rf", output_file + ".refCoding.bed.gz"])
subprocess.call(["rm", "-rf", output_file + ".refCoding.bed.gz.tbi"])


