#! /usr/bin/env python

import sys, re, glob
from annot_utils.junction import *

margin = 30
input_file = sys.argv[1]
output_dir = sys.argv[2]

hout = open(output_dir + "/indel_creation.bed", 'w')
with open(input_file, 'r') as hin:

    header2ind = {}

    header = hin.readline().rstrip('\n').split('\t')
    for (i, cname) in enumerate(header):
        header2ind[cname] = i
    # print '\t'.join(header)

    for line in hin:
        F = line.rstrip('\n').split('\t')
        
        if F[2] == "TCGA-2H-A9GF-01":
            pass

        FF = F[header2ind["Mutation_Key"]].split(',')
        if F[header2ind["Mutation_Type"]] not in ["splicing donor creation", "splicing acceptor creation"]: continue
        if F[header2ind["Splicing_Class"]] not in ["alternative-3'-splice-site", "alternative-5'-splice-site"]: continue

        if FF[1] == "100773560":
            pass

        if len(FF[2]) > 1:
        #  or len(FF[3]) > 1:
            start = int(FF[1]) + 1 - margin
            end = int(FF[1]) + len(FF[2]) - 1 + margin
            print >> hout, FF[0] + '\t' + str(start) + '\t' + str(end) + '\t' + '\t'.join(F)
        elif len(FF[3]) > 1:
            start = int(FF[1]) - margin
            end = int(FF[1]) + margin
            print >> hout, FF[0] + '\t' + str(start) + '\t' + str(end) + '\t' + '\t'.join(F)


hout.close()


make_junction_info(output_dir + "/refJunc.bed", "hg19", True, "1,0", "0,1")

hout = open(output_dir + "/indel_creation_junc.bed", 'w')
subprocess.call(["bedtools", "intersect", "-a", output_dir + "/indel_creation.bed",
                 "-b", output_dir + "/refJunc.bed", "-wa", "-wb"], stdout = hout)
hout.close()



# generate information on
# 1. normal junction, splicing junction, mutation

key2apparent = {}
hout = open(output_dir + "/apparent_splicing.txt", 'w')
with open(output_dir + "/indel_creation_junc.bed", 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')

        exon_intron_junc_pos = int(F[20])
        motif_type = F[22]
        motif_strand = F[23]
        normal_juncs = F[24].split(',')

        splice_key = F[10]
        splice_type = F[11]


        mut_start = int(F[1]) + margin
        mut_end = int(F[2]) - margin

        mut_keys = F[6].split(',')
        mut_ref, mut_alt = mut_keys[2][1:], mut_keys[3][1:]
        if mut_ref == "": mut_ref = "-"
        if mut_alt == "": mut_alt = "-"

        # check on motif splice-type consistency
        if motif_type == "donor" and splice_type == "alternative-3'-splice-site": continue
        if motif_type == "acceptor" and splice_type == "alternative-5'-splice-site": continue

        # just consider indel within exon
        if motif_type == "donor" and motif_strand == '+' or motif_type == "acceptor" and motif_strand == '-':
            if mut_end > exon_intron_junc_pos: 
                print >> hout, '\t'.join(F[3:18])
                key2apparent['\t'.join(F[3:6])] = 1
                continue
        if motif_type == "donor" and motif_strand == '-' or motif_type == "acceptor" and motif_strand == '+':
            if mut_start < exon_intron_junc_pos: 
                print >> hout, '\t'.join(F[3:18])
                key2apparent['\t'.join(F[3:6])] = 1
                continue

        junc_match = re.match(r'([^ \t\n\r\f\v,]+)\:(\d+)\-(\d+)', splice_key)
        if junc_match is None: 
            print "Something is wrong"
            sys.exit(1)
        junc_chr, junc_start, junc_end = junc_match.group(1), int(junc_match.group(2)), int(junc_match.group(3))
 
        if F[6] == "5,67593206,GCTCAAAAGACAGTTTTTCTTCTCTCCT,G":
            pass

        # check the closeness of indel and junction
        if motif_type == "donor" and motif_strand == '+' or motif_type == "acceptor" and motif_strand == '-':
            if abs(junc_start - exon_intron_junc_pos) > 50: continue
        if motif_type == "donor" and motif_strand == '-' or motif_type == "acceptor" and motif_strand == '+':
            if abs(junc_end - exon_intron_junc_pos) > 50: continue


        # check spliced length
        if motif_type == "donor" and motif_strand == '+' or motif_type == "acceptor" and motif_strand == '-':
            spliced_size = exon_intron_junc_pos - junc_start + 1
        if motif_type == "donor" and motif_strand == '-' or motif_type == "acceptor" and motif_strand == '+':
            spliced_size = junc_end - exon_intron_junc_pos + 1

        if abs(spliced_size - len(mut_ref)) <= 3:
            print >> hout, '\t'.join(F[3:18])
            key2apparent['\t'.join(F[3:6])] = 1

hout.close()


hout = open(output_dir + "/filtered_indel_creation.txt", 'w')
with open(input_file, 'r') as hin:

    header2ind = {}

    header = hin.readline().rstrip('\n').split('\t')
    for (i, cname) in enumerate(header):
        header2ind[cname] = i

    for line in hin:
        F = line.rstrip('\n').split('\t')

        FF = F[header2ind["Mutation_Key"]].split(',')
        if F[header2ind["Mutation_Type"]] not in ["splicing donor creation", "splicing acceptor creation"]: continue
        if F[header2ind["Splicing_Class"]] not in ["alternative-3'-splice-site", "alternative-5'-splice-site"]: continue

        FF = F[header2ind["Mutation_Key"]].split(',')
        if len(FF[2]) == 1 and len(FF[3]) == 1: continue
        key = '\t'.join(F[:3])
        if key in key2apparent: continue
        print >> hout, '\t'.join(F)

hout.close()


"""
# temp mutation file
hout = open(output_dir + "/temp_mut.txt", 'w')
print >> hout, F[0] + '\t' + str(mut_start) + '\t' + str(mut_end) + '\t' + mut_ref + '\t' + mut_alt 
hout.close()


# print '\t'.join(F[:7]) + '\t' + F[10] + '\t' + str(exon_intron_junc_pos) + '\t' + str(spliced_size) + '\t' + str(len(mut_ref))
params = ["/home/w3varann/database/GRCh37/GRCh37.fa", 
          "/home/yshira/mysoftware/intron_retention_utils/resource/refGene.txt.gz",
          "--chr_name_list", "/home/yshira/mysoftware/intron_retention_utils/resource/ucsc2grch.txt",
          "--donor_size", "30,0", "--acceptor_size", "0,30", "--template_size", "20", "--template_score_margin", "1", "--debug"]

bam_file = glob.glob("/home/omega3/omega_rna/star/*/" + F[5] + "*/" + F[5] + "*.Aligned.sortedByCoord.out.bam")[0]

# print '\t'.join(F[:7]) + '\t' + F[10] + '\t' + str(exon_intron_junc_pos) + '\t' + str(spliced_size) + '\t' + str(len(mut_ref))
# print '\t'.join(F)
subprocess.call(["intron_retention_utils", "allele_count", bam_file, output_dir + "/temp_mut.txt", output_dir + "/temp_mut_allele_count.txt"] + params)

with open(output_dir + "/temp_mut_allele_count.txt", 'r') as hin:
    for line in hin:
        FFF = line.rstrip('\n').split('\t')
        if FFF[0] == F[4]:
            print '\t'.join(FFF)
"""


"""
9       138839716       138839760       UCEC    UBAC1   TCGA-B5-A0JU-01 9,138839725,ATCTGGAGCTTTCTGGTCTTGTTTTT,A        9:138839726-138839734,-;9:138839751-138839759,- splicing acceptor creation;splicing acceptor disruption canonical;non-canonical 9:138839727-138845525   alternative-3'-splice-site      ---     13.3117 0.0275  0.0277  ---     ---     FALSE   9       138839750       138839751       UBAC1   acceptor
        -       9:138839751-138845526   NM_016172       6

"""



 
