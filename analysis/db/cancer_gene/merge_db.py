#! /usr/bin/env python

gene2lawrence = {}
with open("LawrenceEtAl/LawrenceEtAl_sig_gene.txt", 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        gene2lawrence[F[0]] = F[1]

gene2cgc = {}
with open("CancerGeneCensus/CancerGeneCensus.txt", 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        gene2cgc[F[0]] = F[1].strip('"')

# gene2kinase = {}
# with open("Kinase_com_manual.proc.txt", 'r') as hin:
#     for line in hin:
#         F = line.rstrip('\n').split('\t')
#         gene2kinase[F[0]] = F[1]


gene2vogelstein = {}
with open("Vogelstein_Sicence_2013/Vogelstein_Sicence_2013.txt", 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        gene2vogelstein[F[0]] = F[1] 


gene2ye = {}
with open("YeEtAl/YeEtAl_gene.txt", 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        gene2ye[F[0]] = "CG"

 
genes = list(set(gene2lawrence.keys() + gene2cgc.keys() + gene2vogelstein.keys() + gene2ye.keys()))


hout = open("cancer_gene.txt", 'w')
for gene in sorted(genes):
    lawrence = gene2lawrence[gene] if gene in gene2lawrence else "---"
    cgc = gene2cgc[gene] if gene in gene2cgc else "---"
    vogelstein = gene2vogelstein[gene] if gene in gene2vogelstein else "---"
    ye = gene2ye[gene] if gene in gene2ye else "---"
    print >> hout, gene + '\t' + lawrence + '\t' + cgc + '\t' + vogelstein + '\t' + ye


hout.close()

