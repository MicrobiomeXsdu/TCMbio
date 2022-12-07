
###### kraken2宏基因组无参分析流程 #############
# 1. 多样本并行物种注释
cd /data_200t/shengdahsuang2/TCGA_RNA/CESC/
mkdir -p temp/kraken2_fungi
mkdir -p temp/kraken2_fungi/aligned
conda activate humann
time parallel -j 3 'kraken2 --db /data_200t/shengdahsuang2/db/kraken2_fungi/ \
  --paired --classified-out temp/kraken2_fungi/aligned/{1}#.fq \
  temp/trimqc/{1}_gdc_realn_rehead.bamR1.fq_paired.fq.gz \
  temp/trimqc/{1}_gdc_realn_rehead.bamR2.fq_paired.fq.gz \
  --threads 3 --use-names --report-zero-counts \
  --report temp/kraken2_fungi/{1}_report \
  --output temp/kraken2_fungi/{1}_output --memory-mapping' \
  ::: `cat temp/ls.txt` &> temp/kraken2_fungi/kraken220221117.log


# 2. 利用bracken将kraken结果转成物种丰度
mkdir -p temp/bracken_report_fungi
# conda deactivate
time parallel -j 3 'bracken -d /data_200t/shengdahsuang2/db/kraken2_fungi/ \
  -i temp/kraken2_fungi/{1}_report -o temp/kraken2_fungi/{1}.bracken -w temp/bracken_report_fungi/{1}.bracken.report \
  -r 100 -l S' ::: `cat temp/ls.txt` > temp/bracken_report_fungi/bracken.log
  
# 3. 将Braken的report格式转换成--use-mpa-style格式
time parallel -j 3 \
     'kreport2mpa.py -r temp/bracken_report_fungi/{1}.bracken.report \
      -o temp/bracken_report_fungi/{1}.new.report' \
     ::: `cat temp/ls.txt`

