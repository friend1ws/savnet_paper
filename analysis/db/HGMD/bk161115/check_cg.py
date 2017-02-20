#! /usr/bin/env python

CG = {}
with open("/home/yshira/mysoftware/sv_utils/cancer_gene/cancer_gene.txt", 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        if F[4] == "CG": CG[F[0]] = 1

with open("HGMD.CS.bed", 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        FF = F[5].split(';')

        if FF[4] in CG:
            print '\t'.join(F)


