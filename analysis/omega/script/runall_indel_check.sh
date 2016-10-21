#! /bin/bash

if [ ! -d ../matome/indel_check ]
then
    mkdir -p ../matome/indel_check
fi

python subscript_matome/get_indel_creation.py ../matome/omega.genomon_splicing_mutation.result.txt ../matome/indel_check
 
 
