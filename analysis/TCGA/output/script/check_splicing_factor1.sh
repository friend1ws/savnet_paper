#! /bin/bash

if [ ! -d ../splicing_factor ]
then
    mkdir -p ../splicing_factor
fi

# python subscript_splicing_factor/get_sample_with_sf.py ../../data/mutation > ../splicing_factor/sample_list_with_sfm.txt
python subscript_splicing_factor/get_sample_with_sf.py ../../data/input_list/ 1 0.01 ../splicing_factor/sample_list_with_sfm.txt ../splicing_factor/ctype_sfgene_list.txt

while read line
do
    ctype=`echo ${line} | cut -d ' ' -f 1`
    sfgene=`echo ${line} | cut -d ' ' -f 2`
    echo "qsub subscript_splicing_factor/perform_t_test.sh ${ctype} ${sfgene}"
    qsub subscript_splicing_factor/perform_t_test.sh ${ctype} ${sfgene}

done < ../splicing_factor/ctype_sfgene_list.txt

<<_COM_

qsub subscript_splicing_factor/perform_t_test.sh BLCA SF3B1

qsub subscript_splicing_factor/perform_t_test.sh BLCA ZRSR2

qsub subscript_splicing_factor/perform_t_test.sh BRCA SF3B1

qsub subscript_splicing_factor/perform_t_test.sh CESC U2AF1

qsub subscript_splicing_factor/perform_t_test.sh COAD ZRSR2

qsub subscript_splicing_factor/perform_t_test.sh GBM ZRSR2

qsub subscript_splicing_factor/perform_t_test.sh HNSC ZRSR2

qsub subscript_splicing_factor/perform_t_test.sh LIHC SF3B1

qsub subscript_splicing_factor/perform_t_test.sh LIHC ZRSR2 

qsub subscript_splicing_factor/perform_t_test.sh LUAD U2AF1

qsub subscript_splicing_factor/perform_t_test.sh LUAD SF3B1 

qsub subscript_splicing_factor/perform_t_test.sh LUAD ZRSR2

qsub subscript_splicing_factor/perform_t_test.sh PAAD U2AF1

qsub subscript_splicing_factor/perform_t_test.sh PRAD SF3B1

qsub subscript_splicing_factor/perform_t_test.sh SKCM SF3B1 

qsub subscript_splicing_factor/perform_t_test.sh STAD ZRSR2

qsub subscript_splicing_factor/perform_t_test.sh UCEC U2AF1

qsub subscript_splicing_factor/perform_t_test.sh UCEC SF3B1

qsub subscript_splicing_factor/perform_t_test.sh UCEC ZRSR2

qsub subscript_splicing_factor/perform_t_test.sh UCS U2AF1

qsub subscript_splicing_factor/perform_t_test.sh UVM SF3B1 

_COM_
