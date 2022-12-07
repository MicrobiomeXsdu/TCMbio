####### kraken2进行注释 ######################


## 分析前准备##############
# 进入对应癌症的结果存储文件夹
wd=/mnt/ba76940c-bfd0-42f0-a133-112be8e744c9/TCGA_RESULT/ECAD
mkdir -p $wd && cd $wd
chmod 777 $wd
data=/mnt/ba76940c-bfd0-42f0-a133-112be8e744c9/TCGA_RNA/ECAD/nohumanfastq
# 创建临时和结果目录，放实验设计于结果目录
mkdir -p temp result
chmod 777 temp
chmod 777 result
# 把/mnt/ba76940c-bfd0-42f0-a133-112be8e744c9/TCGA_RNA/ECAD/下的ls.txt 拷贝到/mnt/ba76940c-bfd0-42f0-a133-112be8e744c9/TCGA_RESULT/ECAD
cp /mnt/ba76940c-bfd0-42f0-a133-112be8e744c9/TCGA_RNA/ECAD/ls.txt /mnt/ba76940c-bfd0-42f0-a133-112be8e744c9/TCGA_RESULT/ECAD
# 激活工作环境
conda activate

######## 一. FastQC质量评估 #################
# 1.fastqc每个文件一个线程，6个双端样本12个文件，设置6线程
# 结果见result/fastqc目录，解读见[数据的质量控制软件——fastQC](https://mp.weixin.qq.com/s/MMfierO-8H2MVRkJKGCOVQ)
mkdir result/fastqc
time fastqc ${data}/*.fq -t 20 -o result/fastqc 

# 2.生成多样品报告比较
# 查看右侧result/qc目录中multiqc_report.html，可交互式报告
mkdir  result/qc
multiqc -d result/fastqc/ -o result/qc

# 3.首先判断测序类型，随便选一个下载的测序文件即可
## ***记得把fddad662-6df4-44d9-91a4-b99139c8460a_gdc_realn_rehead.bamR2.fq更换成自己的测序文件
qs=/mnt/ba76940c-bfd0-42f0-a133-112be8e744c9/TCGA_RNA/ECAD/nohumanfastq/fddad662-6df4-44d9-91a4-b99139c8460a_gdc_realn_rehead.bamR2.fq
less $qs | head -n 1000 | awk '{if(NR%4==0) printf("%s",$0);}' \
| od -A n -t u1 -v \
| awk 'BEGIN{min=100;max=0;} \
{for(i=1;i<=NF;i++) {if($i>max) max=$i; if($i<min) min=$i;}}END \
{if(max<=126 && min<59) print "Phred33"; \
else if(max>73 && min>=64) print "Phred64"; \
else if(min>=59 && min<64 && max>73) print "Solexa64"; \
else print "Unknown score encoding"; \
print "( " min ", " max, ")";}'
# Phred33

# 4. 利用trimmomatic对测序read进行质控修剪
mkdir -p temp/trimqc
cd ${data}
time parallel -j 4 --xapply \
  'trimmomatic PE -threads 8 -phred33 \
  -summary /mnt/ba76940c-bfd0-42f0-a133-112be8e744c9/TCGA_RESULT/ECAD/temp/trimqc/{1}sum.txt \
  {1} {2} \
  /mnt/ba76940c-bfd0-42f0-a133-112be8e744c9/TCGA_RESULT/ECAD/temp/trimqc/{1}_paired.fq.gz \
  /mnt/ba76940c-bfd0-42f0-a133-112be8e744c9/TCGA_RESULT/ECAD/temp/trimqc/{1}_unpaired.fq.gz \
  /mnt/ba76940c-bfd0-42f0-a133-112be8e744c9/TCGA_RESULT/ECAD/temp/trimqc/{2}_paired.fq.gz \
  /mnt/ba76940c-bfd0-42f0-a133-112be8e744c9/TCGA_RESULT/ECAD/temp/trimqc/{2}_unpaired.fq.gz \
  ILLUMINACLIP:/mnt/ba76940c-bfd0-42f0-a133-112be8e744c9/miniconda3/share/trimmomatic/adapters/TruSeq2-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:35' \
  ::: *R1.fq ::: *R2.fq

  
# 5. 合并质控后的双端文件
mkdir temp/concats
cat ls.txt | while read line
do
cat temp/trimqc/${line}_gdc_realn_rehead.bamR1.fq_paired.fq.gz temp/trimqc/${line}_gdc_realn_rehead.bamR2.fq_paired.fq.gz > temp/concats/${line}.fastq.gz
done

# 6. 质控后质量再评估 (可选)
fastqc temp/trimqc/*_gdc_realn_rehead.bamR?.fq_paired.fq.gz -t 6 -o temp/trimqc/
multiqc -d temp/trimqc/ -o result/qc/trimqc/


###### 二、kraken2宏基因组无参分析流程 #############
# 1. 多样本并行物种注释
cd /mnt/ba76940c-bfd0-42f0-a133-112be8e744c9/TCGA_RESULT/ECAD/
mkdir -p temp/kraken2
time parallel -j 3 'kraken2 --db /data_200t/shareAll/zhaolanlan/db/kraken2/210122/ \
  --paired --classified-out {1}#.fq temp/trimqc/{1}_gdc_realn_rehead.bamR1.fq_paired.fq.gz \
  temp/trimqc/{1}_gdc_realn_rehead.bamR2.fq_paired.fq.gz --threads 3 --use-names --report-zero-counts \
  --report temp/kraken2/{1}_report \
  --output temp/kraken2/{1}_output --memory-mapping' \
  ::: `cat ls.txt` &> temp/kraken2/kraken2.log

# 2. 利用bracken将kraken结果转成物种丰度
mkdir temp/bracken_report
time parallel -j 3 'bracken -d /data_200t/shareAll/zhaolanlan/db/kraken2/210122/ \
  -i temp/kraken2/{1}_report -o temp/kraken2/{1}.bracken -w temp/bracken_report/{1}.bracken.report \
  -r 100 -l S' ::: `cat ls.txt` > temp/bracken.log
  
# 3. 将Braken的report格式转换成--use-mpa-style格式
time parallel -j 3 'kreport2mpa.py -r temp/bracken_report/{1}.bracken.report   -o temp/bracken_report/{1}.new.report' ::: `cat ls.txt`

