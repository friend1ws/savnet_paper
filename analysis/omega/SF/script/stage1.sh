#! /bin/bash

python get_sample_with_sf.py > ../output/sample_list_with_sfm.txt

qsub perform_t_test.sh BLCA SF3B1

qsub perform_t_test.sh BRCA SF3B1

qsub perform_t_test.sh CESC U2AF1

qsub perform_t_test.sh LUAD U2AF1

qsub perform_t_test.sh PAAD U2AF1

qsub perform_t_test.sh SKCM SKCM

qsub perform_t_test.sh UCEC U2AF1

qsub perform_t_test.sh UCEC SF3B1

qsub perform_t_test.sh UCS U2AF1

qsub perform_t_test.sh UVM SF3B1 


