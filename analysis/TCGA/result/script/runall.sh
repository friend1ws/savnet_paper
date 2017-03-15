#! /bin/bash

if [ ! -d ../temporary ]
then
    mkdir ../temporary
fi

if [ ! -d ../figure ]
then
    mkdir ../figure
fi

if [ ! -d ../table ]
then
    mkdir ../table
fi

# Splicing Factor related table
cp ../../output/splicing_factor/TCGA.SJ_IR.sig_summary.txt ../table/
cp ../../output/splicing_factor/TCGA.savnet.sf_mut.summary.txt ../table/

# spling summary table
cut -f1-14 ../../output/rescue/TCGA.savnet.with_rescued.result.txt | grep -v "Rescued" > ../table/TCGA.savnet.result.proc.txt 

# rescued table
cut -f1-14 ../../output/rescue/TCGA.savnet.with_rescued.result.txt | grep "Sample_Name\|Rescued" > ../table/TCGA.savnet.rescued.proc.txt

# SNV data 
python subscript_matome/summarize_snv_info.py \
    ../../output/savnet_out/d3.6_a6.1_8_ka/TCGA.savnet.result.txt \
    ../temporary/TCGA.savnet.motif_summary.txt \
    ../../../db/GRCh37/GRCh37.fa \
    ../../../db/spidex/hg19_spidex.bed.gz

Rscript subscript_matome/add_mes_score.R ../temporary/TCGA.savnet.motif_summary.txt ../temporary/TCGA.savnet.motif_summary.mes.txt 


# mutation data
python subscript_matome/summarize_mut_count.py \
    ../../output/savnet_out/d3.6_a6.1_8_ka/TCGA.savnet.result.txt \
    ../../data/input_list \
    ../temporary/omega.mut_count.txt

cp ../temporary/omega.mut_count.txt ../table/omega.mut_count.txt


# position wise false positive check
python subscript_matome/check_position_fdr.py ../../output/savnet_out/d5.15_a15.5_8_ka/ 3.0 > ../temporary/position_fdr.txt
Rscript subscript_matome/position_fdr.R

# cancer type wise false positive check
python subscript_matome/check_cancer_type_fdr.py ../../output/savnet_out/d3.6_a6.1_8_ka/ 3.0 > ../temporary/cancer_type_fdr.txt
Rscript subscript_matome/cancer_type_fdr.R 


# intron retention VAF check
python subscript_matome/add_allele_count.py \
    ../../output/savnet_out/d3.6_a6.1_8_ka/TCGA.savnet.result.txt \
    ../../data/allele_count \
    ../temporary/TCGA.savnet.allele_count.txt

Rscript subscript_matome/IR_VAF.R


# comparison with Jung et al.
Rscript subscript_matome/JungEtAl_venn.R
rm -rf ../figure/JungEtAl*.log


# category count summary
Rscript subscript_matome/category_count.R

# cancer type count summary
Rscript subscript_matome/cancer_type_count_summary.R 



# SAV generating multi splicing changes
Rscript subscript_matome/multi_splice_sav_count.R
Rscript subscript_matome/multi_splice_sav_example.R


# allele count summary
python subscript_matome/summarize_allele_count.py \
    ../../output/savnet_out/d3.6_a6.1_8_ka/TCGA.savnet.result.txt \
    ../../data/allele_count \
    ../../data/expression \
    ../../../db/GRCh37/GRCh37.fa \
    ../temporary/TCGA.savnet.allele_count.summary.txt

Rscript subscript_matome/splicing_mut_count.R
Rscript subscript_matome/mes_hb_diff_score.R
Rscript subscript_matome/seqlogo_summary.R
Rscript subscript_matome/snv_pos_total_gsm_ratio.R


# splicing motif creating snv  distribution 
Rscript subscript_matome/snv_motif_dist_creation.R


# signature related figures
python subscript_matome/summarize_all_mut_membership.py ../../data/pmsignature/ ../temporary/TCGA.mut_membership.all.result.txt

python subscript_matome/summarize_all_mut_membership_savnet.py \
    ../../output/savnet_out/d3.6_a6.1_8_ka/TCGA.savnet.result.txt \
    ../temporary/TCGA.mut_membership.savnet.result.txt \
    ../../data/pmsignature/

Rscript subscript_matome/signature_membership.R


# alternative junction
python subscript_matome/check_alt_junc.py \
    ../../output/savnet_out/d3.6_a6.1_8_ka/TCGA.savnet.result.txt \
    ../temporary/TCGA.savnet.alt_junc.txt

Rscript subscript_matome/alt_junc_pos.R


# GC contents
python subscript_matome/summarize_intron_info.py \
    ../temporary/TCGA.savnet.allele_count.summary.txt \
    ../temporary/TCGA.savnet.gc_intron.txt \
    ../../../db/GRCh37/GRCh37.fa

Rscript subscript_matome/GC_contents_diff.R
Rscript subscript_matome/GC_Len_test_comp.R


# splicing pattern wise relative expression 
Rscript subscript_matome/mut_func_exp.R 


# cancer gene ratio
python subscript_matome/merge_mut.py \
    ../../output/savnet_out/d3.6_a6.1_8_ka/TCGA.savnet.result.txt \
    ../../data/input_list > ../temporary/TCGA.mutation.merged.txt


cat ../../output/savnet_out/d3.6_a6.1_8_ka/TCGA.savnet.result.txt \
    ../../output/rescue/TCGA.savnet.rescued.result.txt > \
    ../../output/rescue/TCGA.savnet.with_rescued.result.txt

Rscript subscript_matome/cancer_gene_ratio.R
Rscript subscript_matome/cancer_gene_ratio2.R 
Rscript subscript_matome/cancer_gene_summary.R
Rscript subscript_matome/cancer_gene_ratio_table.R

# profile
Rscript subscript_matome/sav_profile.R


# motif read count
python subscript_matome/gather_read_num.py ../../output/rescue/TCGA.savnet.with_rescued.result.txt ../../data/input_list/ ../temporary/TCGA.motif_read_num.txt

Rscript subscript_matome/top_splicing_read_ratio.R 




