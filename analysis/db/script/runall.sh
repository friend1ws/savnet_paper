#! /bin/bash
#$ -S /bin/bash
#$ -cwd

# cancer gene
# see the memo.txt in the cancer_gene directory

# HGMD
python make_HGMD_CS_bed.py > ../HGMD/HGMD.CS.bed

bgzip -f ../HGMD/HGMD.CS.bed

tabix -p bed ../HGMD/HGMD.CS.bed.gz

# spidex

unzip  ../spidex/hg19_spidex.zip -d ../spidex/

python spidex_to_bed.py ../spidex/hg19_spidex.txt > ../spidex/hg19_spidex.bed

bgzip -f ../spidex/hg19_spidex.bed

tabix -p bed ../spidex/hg19_spidex.bed.gz

rm -rf ../spidex/hg19_spidex.txt*
rm -rf ../spidex/DeepGenomics_ANNOVAR_EULA_2.2.docx


