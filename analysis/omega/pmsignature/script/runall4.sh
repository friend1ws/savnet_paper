#! /bin/bash

python summarize_all_mut_membership.py > ../output/omega.mut_membership.all.result.txt

python summarize_all_mut_membership_gsm.py ../output/omega.mut_membership.result.txt

Rscript signature_membership.R

