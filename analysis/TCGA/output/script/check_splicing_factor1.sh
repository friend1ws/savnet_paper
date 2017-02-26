#! /bin/bash

if [ ! -d ../splicing_factor ]
then
    mkdir -p ../splicing_factor
fi

python subscript_splicing_factor/get_sample_with_sf.py ../../data/mutation > ../splicing_factor/sample_list_with_sfm.txt

qsub subscript_splicing_factor/perform_t_test.sh BLCA SF3B1

qsub subscript_splicing_factor/perform_t_test.sh BRCA SF3B1

qsub subscript_splicing_factor/perform_t_test.sh CESC U2AF1

qsub subscript_splicing_factor/perform_t_test.sh LUAD U2AF1

qsub subscript_splicing_factor/perform_t_test.sh PAAD U2AF1

qsub subscript_splicing_factor/perform_t_test.sh SKCM SKCM

qsub subscript_splicing_factor/perform_t_test.sh UCEC U2AF1

qsub subscript_splicing_factor/perform_t_test.sh UCEC SF3B1

qsub subscript_splicing_factor/perform_t_test.sh UCS U2AF1

qsub subscript_splicing_factor/perform_t_test.sh UVM SF3B1 


