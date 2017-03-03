#! /bin/sh

python subscript_rescue/merge_hotspot_result.py ../rescue/sav.anno.txt ../rescue ../rescue/TCGA.savnet.hotspot.result.txt

python subscript_rescue/rescue_hotspot.py ../../data/input_list/ ../savnet_out/d3.6_a6.1_8_ka/TCGA.savnet.result.txt ../rescue/TCGA.savnet.hotspot.result.txt > ../rescue/TCGA.savnet.rescued.result.txt

cat ../savnet_out/d3.6_a6.1_8_ka/TCGA.savnet.result.txt ../rescue/TCGA.savnet.rescued.result.txt > ../rescue/TCGA.savnet.with_rescued.result.txt

