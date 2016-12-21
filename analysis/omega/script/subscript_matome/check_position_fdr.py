#! /usr/bin/env python

import sys, glob, os, re

gsm_dir = sys.argv[1]
score_thres = sys.argv[2]

all_gsm_perm_file = glob.glob(gsm_dir + "*/*.splicing_mutation.count_summary.anno.perm_all.txt") + \
                    glob.glob(gsm_dir + "*/*.genomon_splicing_mutation.result.txt")

key2count = {}
for gsm_perm_file in sorted(all_gsm_perm_file):

    cancer_type = os.path.basename(os.path.dirname(gsm_perm_file))
    file_type = "permutation" if gsm_perm_file.endswith("splicing_mutation.count_summary.anno.perm_all.txt") else "original"

    for rel_pos in range(1, 21):
        key2count[cancer_type + '\t' + file_type + '\t' + "splicing donor disruption" + '\t' + str(rel_pos)] = 0
        key2count[cancer_type + '\t' + file_type + '\t' + "splicing donor creation" + '\t' + str(rel_pos)] = 0
        key2count[cancer_type + '\t' + file_type + '\t' + "splicing acceptor disruption" + '\t' + str(rel_pos)] = 0
        key2count[cancer_type + '\t' + file_type + '\t' + "splicing acceptor creation" + '\t' + str(rel_pos)] = 0

    with open(gsm_perm_file, 'r') as hin:
        header2ind = {}
        header = hin.readline().rstrip('\n').split('\t')
        for (i, cname) in enumerate(header):
            header2ind[cname] = i

        for line in hin:
            F = line.rstrip('\n').split('\t')
            if float(F[header2ind["Score"]]) < float(score_thres): continue
            FF = F[header2ind["Mutation_Key"]].split(',')
            if len(FF[2]) > 1 or len(FF[3]) > 1: continue

            pos_match = re.match(r'([\w\d]+)\:(\d+)\-(\d+),([\+\-])', F[header2ind["Motif_Pos"]])
            motif_chr, motif_start, motif_end, motif_strand = pos_match.group(1), pos_match.group(2), pos_match.group(3), pos_match.group(4)

            if motif_strand == '+':
                rel_pos = int(FF[1]) - int(motif_start) + 1
            elif motif_strand == '-':
                rel_pos = int(motif_end) - int(FF[1]) + 1
            else:
                print >> sys.stderr, "something is wrong"
                sys.exit(1)

            key = cancer_type + '\t' + file_type + '\t' + F[header2ind["Mutation_Type"]] + '\t' + str(rel_pos)

            if "donor" in F[header2ind["Mutation_Type"]]:
                if F[header2ind["Is_Canonical"]] == "canonical" and rel_pos not in [6, 7]:
                    print >> sys.stderr, "donor is wrong!"
                    print >> sys.stderr, '\t'.join(F)
                    print >> sys.stderr, key
                    sys.exit(1)
                elif F[header2ind["Is_Canonical"]] == "non-canonical" and rel_pos in [6, 7]:
                    print >> sys.stderr, "donor is wrong!"
                    print >> sys.stderr, '\t'.join(F)
                    print >> sys.stderr, key
                    sys.exit(1)
                else:
                    pass

            elif "acceptor" in F[header2ind["Mutation_Type"]]:
                if F[header2ind["Is_Canonical"]] == "canonical" and rel_pos not in [14, 15]:
                    print >> sys.stderr, "acceptor is wrong!"
                    print >> sys.stderr, '\t'.join(F)
                    print >> sys.stderr, key
                    sys.exit(1)
                elif F[header2ind["Is_Canonical"]] == "non-canonical" and rel_pos in [14, 15]:
                    print >> sys.stderr, "acceptor is wrong!"
                    print >> sys.stderr, '\t'.join(F)
                    print >> sys.stderr, key
                    sys.exit(1)
                else:
                    pass
            else:
                print >> sys.stderr, "something is wrong"
                sys.exit(1)


            key2count[key] = key2count[key] + 1


print "Cancer_Type" + '\t' + "Is_Original" + '\t' + "Mutation_Type" + '\t' + "Rel_Pos" + '\t' + "Count"
for key in sorted(key2count):
    print key + '\t' + str(key2count[key])



