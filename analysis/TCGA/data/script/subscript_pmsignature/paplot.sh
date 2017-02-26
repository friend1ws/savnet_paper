# !/bin/bash
#$ -S /bin/bash         # set shell in UGE
#$ -cwd                 # execute at the submitted dir

# set python environment
# export PYTHONHOME=/usr/local/package/python2.7/current
# export PATH=$PYTHONHOME/bin:$PATH
export PYTHONPATH=/home/w3varann/.genomon_local/genomon_pipeline-2.5.0/python2.7-packages/lib/python

CTYPE=$1

for i in `seq 2 6`
do
    /home/w3varann/.genomon_local/genomon_pipeline-2.5.0/python2.7-packages/bin/paplot pmsignature ../pmsignature/${CTYPE}/pmsignature/pmsignature.ind.result.${i}.json ../pmsignature/${CTYPE}/paplot SAVNet_paper --config_file /home/w3varann/.genomon_local/genomon_pipeline-2.5.0/genomon_conf/paplot/paplot_dna.cfg --title 'Pmsignature' --overview 'Pmsignature type=ind.' --ellipsis ind${i}

    /home/w3varann/.genomon_local/genomon_pipeline-2.5.0/python2.7-packages/bin/paplot signature ../pmsignature/${CTYPE}/pmsignature/pmsignature.full.result.${i}.json ../pmsignature/${CTYPE}/paplot SAVNet_paper --config_file /home/w3varann/.genomon_local/genomon_pipeline-2.5.0/genomon_conf/paplot/paplot_dna.cfg --title 'Signature' --overview 'Pmsignature type=full.' --ellipsis full${i}
done


/home/w3varann/.genomon_local/genomon_pipeline-2.5.0/python2.7-packages/bin/paplot index ../pmsignature/${CTYPE}/paplot --config_file /home/w3varann/.genomon_local/genomon_pipeline-2.5.0/genomon_conf/paplot/paplot_dna.cfg

