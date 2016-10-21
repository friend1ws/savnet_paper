#! /usr/bin/env python

import sys, re, numpy
# from subscript_matome import my_seq
import my_seq

input_file = sys.argv[1]
output_file = sys.argv[2]
reference = sys.argv[3]


key2exists = {}
header2ind = {}

seq_margin = 100

splicing_donor_motif = "AGGTRAGT"
splicing_acceptor_motif = "YYYYNCAGG"


nuc2vec = {'A': [1, 0, 0, 0], 'C': [0, 1, 0, 0], 'G': [0, 0, 1, 0], 'T': [0, 0, 0, 1], \
           'W': [1, 0, 0, 1], 'S': [0, 1, 1, 0], 'M': [1, 1, 0, 0], 'K': [0, 0, 1, 1], \
           'R': [1, 0, 1, 0], 'Y': [0, 1, 0, 1], 'B': [0, 1, 1, 1], 'D': [1, 0, 1, 1], \
           'H': [1, 1, 0, 1], 'V': [1, 1, 1, 0], 'N': [1, 1, 1, 1]}


def get_edit_dist(seq1, seq2):

    if len(seq1) != len(seq2):
        print >> sys.stderr, "The lengths of seq1 and seq2 are different"
        sys.exit(1)

    edit_dist = 0
    for i in range(len(seq1)):
        edit_dist = edit_dist + 1 - numpy.dot(nuc2vec[seq1[i]], nuc2vec[seq2[i]])

    return(edit_dist)


hout = open(output_file, 'w')
print >> hout, '\t'.join(["Cancer_Type", "Sample_Name", "Gene_Symbol", "Mutation_Key", 
                          "Motif_Seq", "Rel_Pos", "Ref_Base", "Alt_Base", "Motif_Type", "Is_Canonocal"])
 
with open(input_file, 'r') as hin:

    header = hin.readline().rstrip('\n').split('\t')
    for i in range(len(header)):
        header2ind[header[i]] = i

    for line in hin:
        F = line.rstrip('\n').split('\t')
        if F[header2ind["Mutation_Type"]] not in ["splicing donor creation", "splicing acceptor creation"]: continue
        key = F[header2ind["Cancer_Type"]] + '\t' + F[header2ind["Sample_Name"]] + '\t' + F[header2ind["Gene_Symbol"]] + '\t' + F[header2ind["Mutation_Key"]]
        if key in key2exists: continue


        if F[header2ind["Sample_Name"]] == "TCGA-DU-A7TB-01":
            pass

        mchr, mpos, mref, malt = F[header2ind["Mutation_Key"]].split(',')
        if len(mref) > 1 or len(malt) > 1: continue

        motif_pos, strand = F[header2ind["Motif_Pos"]].split(',')

        pos_match = re.match(r'([\w\d]+)\:(\d+)\-(\d+)', motif_pos)
        schr, sstart, send = pos_match.group(1), pos_match.group(2), pos_match.group(3)

        motif_size = int(send) - int(sstart) + 1
        seq_get_pos = schr + ':' + str(int(sstart) - seq_margin) + '-' + str(int(send) + seq_margin)

        # if strand == "+":
        ref_seq = my_seq.get_seq(reference, motif_pos)
        tmp_seq = my_seq.get_seq(reference, seq_get_pos)

        mut_start = int(mpos) - (int(sstart) - seq_margin)
        mut_end = int(mpos) - (int(sstart) - seq_margin) + len(mref) - 1

        if tmp_seq[mut_start:(mut_end + 1)] != mref:
            print "inconsinstent:" + '\t' + mref + '\t' + tmp_seq[mut_start:(mut_end + 1)]
            sys.exit(1)

        alt_seq = tmp_seq[:mut_start] + malt + tmp_seq[(mut_end + 1):]

        if F[header2ind["Mutation_Type"]] == "splicing donor creation" and strand == '+' or \
           F[header2ind["Mutation_Type"]] == "splicing acceptor creation" and strand == '-':
            alt_seq = alt_seq[(seq_margin + 0):(seq_margin + motif_size + 0)]
        else:
            alt_seq = alt_seq[(- seq_margin - motif_size):(-seq_margin)]


        if strand == "-":
            ref_seq = my_seq.reverse_complement(ref_seq)
            alt_seq = my_seq.reverse_complement(alt_seq)

        if F[header2ind["Mutation_Type"]] == "splicing donor creation":
            edit_dist_ref = get_edit_dist(ref_seq, splicing_donor_motif)
            edit_dist_alt = get_edit_dist(alt_seq, splicing_donor_motif)
        else:
            edit_dist_ref = get_edit_dist(ref_seq, splicing_acceptor_motif) 
            edit_dist_alt = get_edit_dist(alt_seq, splicing_acceptor_motif)


        print key + '\t' + ref_seq + '\t' + alt_seq + '\t' + str(edit_dist_ref) + '\t' + str(edit_dist_alt) + '\t' + F[header2ind["Motif_Pos"]] + '\t' + F[header2ind["Mutation_Type"]] + '\t' + F[header2ind["Is_Canonical"]]

hout.close()


