#! /usr/bin/env python

import sys, os, glob

output_file = sys.argv[1]
mut_dir = sys.argv[2]
sj_dir1 = sys.argv[3]
sj_dir2 = sys.argv[4]
ir_dir1 = sys.argv[5]
ir_dir2 = sys.argv[6]

mut_files = glob.glob(mut_dir + "/*.genomon.genomon_mutation.result.filt.blacklist_filtered.txt")
sj_files = glob.glob(sj_dir1 + "/*/*.SJ.out.tab") + glob.glob(sj_dir2 + "/*/*.SJ.out.tab")
ir_files = glob.glob(ir_dir1 + "/*/*.ir_simple_count.txt") + glob.glob(ir_dir2 + "/*/*.ir_simple_count.txt")


def get_weight(log_final_file):

    with open(log_final_file) as hin:
        for line in hin:
            if "Uniquely mapped reads number |" in line:
                read_num = int(line.replace("Uniquely mapped reads number |", '').strip(' '))
                weight = float(read_num) / float(10000000)
                break
    
    return(weight) 


sample2mut_file = {}
for mut_file in sorted(mut_files):
    # sample = os.path.basename(os.path.dirname(mut_file))
    sample = os.path.basename(mut_file).replace(".genomon.genomon_mutation.result.filt.blacklist_filtered.txt", '')
    sample = sample[:15]
    sample2mut_file[sample] = mut_file

sample2sj_file = {}
sample2weight = {}
for sj_file in sorted(sj_files):
    sample = os.path.basename(os.path.dirname(sj_file))
    sample = sample[:15]
    sample2sj_file[sample] = sj_file
    log_final_file = sj_file.replace("SJ.out.tab", "Log.final.out")
    weight = get_weight(log_final_file)
    sample2weight[sample] = weight


sample2ir_file = {}
for ir_file in sorted(ir_files):
    sample = os.path.basename(ir_file).replace(".genomon_intron_retention.result.txt", '')
    sample = sample[:15]
    sample2ir_file[sample] = ir_file

hout = open(output_file, 'w')
print >> hout, '\t'.join(["Sample_Name", "Weight", "Mutation_File", "SJ_File", "IR_File"])
for sample in sorted(sample2mut_file):
    if sample not in sample2sj_file: continue
    if sample not in sample2ir_file: continue
    print >> hout, sample + '\t' + str(round(sample2weight[sample], 4)) + '\t' + sample2mut_file[sample] + '\t' + sample2sj_file[sample] + '\t' + sample2ir_file[sample]

hout.close()



