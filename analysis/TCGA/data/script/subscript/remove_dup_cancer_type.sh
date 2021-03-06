#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e log/ -o log/

CTYPE=$1

if [ ! -d ../mutation/${CTYPE} ]
then
    mkdir -p ../mutation/${CTYPE}
fi

for mutfile in `ls /home/kchiba/work_directory/work_hotspot/black_list_output_min10/${CTYPE}/*.genomon.genomon_mutation.result.filt.blacklist_filtered.txt`
do
    bmutfile=`basename $mutfile`
    SAMPLE=${bmutfile%.genomon.genomon_mutation.result.filt.blacklist_filtered.txt}
    mSAMPLE=${SAMPLE%-0?}

    mut_file=/home/kchiba/work_directory/work_hotspot/black_list_output_min10/${CTYPE}/${SAMPLE}.genomon.genomon_mutation.result.filt.blacklist_filtered.txt
    oxog_file=/home/kchiba/work_directory/work_hotspot/work_oxog/script_black_list_output_min10/result/${CTYPE}/${mSAMPLE}_oxog.result_test.result.txt

    if [ -f ${oxog_file} ]
    then
        echo "python subscript/remove_dup_oxog.py ${mut_file} ${oxog_file} ../mutation/canonical.bed.gz ../mutation/noncanonical.bed.gz > ../mutation/${CTYPE}/${SAMPLE}.mutation.filt.txt"
        python subscript/remove_dup_oxog.py ${mut_file} ${oxog_file} ../mutation/canonical.bed.gz ../mutation/noncanonical.bed.gz > ../mutation/${CTYPE}/${SAMPLE}.mutation.filt.txt
    fi

done


