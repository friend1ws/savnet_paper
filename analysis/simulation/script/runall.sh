#! /bin/bash

if [ ! -d ../output ]
then
    mkdir ../output
fi

echo "python simulation1.py"
python simulation1.py

echo "python simulation2.py"
python simulation2.py

echo "Rscript get_boxplot.R"
Rscript get_boxplot.R

