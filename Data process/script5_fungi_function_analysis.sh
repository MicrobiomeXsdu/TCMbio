cd /data_200t/shengdahsuang2/TCGA_RNA/CESC
# 一、组装分析流程 Assemble-based
## 1.1 拼接 Assembly
# 提取微生物序列，输入的是kraken下游文件
# ls.txt文件是TCGA文件编号
head temp/ls.txt

# 从fastq文件中剔除比对至人的序列
mkdir -p temp/kraken2_fungi/classified_micro
cat temp/ls.txt | while read line
do
cat temp/kraken2_fungi/aligned/${line}_1.fq | sed '/kraken:taxid|131567/,+3d' | sed '/kraken:taxid|9606/,+3d' > temp/kraken2_fungi/classified_micro/${line}_1_micro.fq
cat temp/kraken2_fungi/aligned/${line}_2.fq | sed '/kraken:taxid|131567/,+3d' | sed '/kraken:taxid|9606/,+3d' > temp/kraken2_fungi/classified_micro/${line}_2_micro.fq
done

#注意文件路径，去除其他的保留微生物

### 1.1.1 MEGAHIT拼接  ###py2环境
# 有旧文件夹时megahit无法运行
# 组装，62m54.445s，TB级数据需几天至几周
mkdir -p function_analysis_fungi && cd function_analysis_fungi
mkdir -p result
mkdir -p script
mkdir -p temp
conda activate humann
# ls2.txt文件与ls.txt是同一个
# 这一步是把编号补充成完整的路径+文件名，直接进行拼接
time megahit -t 3 \
-1 `cat /data_200t/shengdahsuang2/TCGA_RNA/CESC/temp/ls.txt | sed 's/^/..\/temp\/kraken2_fungi\/classified_micro\//;s/$/_1_micro.fq/'| tr '\n' ','|sed 's/,$//'` \
-2 `cat /data_200t/shengdahsuang2/TCGA_RNA/CESC/temp/ls.txt | sed 's/^/..\/temp\/kraken2_fungi\/classified_micro\//;s/$/_2_micro.fq/'| tr '\n' ','|sed 's/,$//'` \
-o /data_200t/shengdahsuang2/TCGA_RNA/CESC/function_analysis_fungi/temp/megahit \
&> /data_200t/shengdahsuang2/TCGA_RNA/CESC/function_analysis_fungi/temp/megahit20220828.log

# 检查点：查看拼接结果大小，通常300M~5G
ls -sh /data_200t/shengdahsuang2/TCGA_RNA/CESC/function_analysis_fungi/temp/megahit/final.contigs.fa

cd /data_200t/shengdahsuang2/TCGA_RNA/CESC/function_analysis_fungi/

# 统计
seqkit stat temp/megahit/final.contigs.fa
#file                           format  type  num_seqs    sum_len  min_len  avg_len  max_len
#temp/megahit/final.contigs.fa  FASTA   DNA      6,477  6,660,634      200  1,028.4   38,823

# 预览重叠群最前6行，前60列字符
head -n6 temp/megahit/final.contigs.fa | cut -c1-60
# 如果contigs太多，可以按长度筛选，降低数据量，提高基因完整度，详见附录megahit
# 备份重要结果
mkdir -p result/megahit/
ln -f temp/megahit/final.contigs.fa result/megahit/
# 删除临时文件
rm -rf temp/megahit/intermediate_contigs/
  
### 1.1.3 QUAST评估 py2环境
quast.py temp/megahit/final.contigs.fa -o result/megahit/quast -t 2
# 生成report文本tsv/txt、网页html、PDF等格式报告

## 1.2 基因预测、去冗余和定量 
# Gene prediction, cluster & quantitfy
### 1.2.1 metaProdigal基因预测
# 输入文件：拼装好的序列文件 result/megahit/final.contigs.fa
# 输出文件：prodigal预测的基因序列 temp/prodigal/gene.fa
# 基因文件大，可参考附录prodigal拆分基因文件，并行计算
# conda activate py2
mkdir -p temp/prodigal
# prodigal的meta模式预测基因，35s，>和2>&1记录分析过程至gene.log
time prodigal -i result/megahit/final.contigs.fa \
-d temp/prodigal/gene.fa \
-o temp/prodigal/gene.gff \
-p meta -f gff > temp/prodigal/gene.log 2>&1 

# 查看日志是否运行完成，有无错误
tail temp/prodigal/gene.log
# 统计基因数量（10702)
grep -c '>' temp/prodigal/gene.fa 
# 统计完整基因数量，数据量大可只用完整基因部分(2913)
grep -c 'partial=00' temp/prodigal/gene.fa 
# # 提取完整基因(完整片段获得的基因全为完整，如成环的细菌基因组)
# grep 'partial=00' temp/prodigal/gene.fa | cut -f1 -d ' '| sed 's/>//' > temp/prodigal/full_length.id
# seqkit grep -f temp/prodigal/full_length.id temp/prodigal/gene.fa > temp/prodigal/full_length.fa
# seqkit stat temp/prodigal/full_length.fa

### 1.2.2 基因聚类/去冗余cd-hit
# 输入文件：prodigal预测的基因序列 temp/prodigal/gene.fa
# 输出文件：去冗余后的基因和蛋白序列：result/NR/nucleotide.fa
#                                     result/NR/protein.fa

mkdir -p result/NR
# aS覆盖度，c相似度，G局部比对，g最优解，T多线程，M内存0不限制
# 2万基因2m，2千万需要2000h，多线程可加速
time cd-hit-est -i temp/prodigal/gene.fa \
-o result/NR/nucleotide.fa \
-aS 0.9 -c 0.95 -G 0 -g 0 -T 0 -M 0
# 统计非冗余基因数量，单次拼接结果数量下降不大，多批拼接冗余度高(10695)
grep -c '>' result/NR/nucleotide.fa
# 翻译核酸为对应蛋白序列，emboss
transeq -sequence result/NR/nucleotide.fa \
-outseq result/NR/protein.fa -trim Y 
# 序列名自动添加了_1，为与核酸对应要去除
sed -i 's/_1 / /' result/NR/protein.fa
# 两批数据去冗余使用cd-hit-est-2d加速，见附录

### 1.2.3 基因定量salmon
# 输入文件：去冗余后的基因和蛋白序列：result/NR/nucleotide.fa
# 输出文件：Salmon定量后的结果：result/salmon/gene.count
#                               result/salmon/gene.TPM
conda activate salmon
mkdir -p temp/salmon
salmon -v # 1.8.0

# 建索引, -t序列, -i 索引，10s
time salmon index \
-t result/NR/nucleotide.fa \
-p 9 \
-i temp/salmon/index 

# 定量，l文库类型自动选择，p线程，--meta宏基因组模式, 2个任务并行2个样
# 注意parallel中待并行的命令必须是双引号，内部变量需要使用原始绝对路径 

time parallel -j 2 \
"salmon quant \
    -i /data_200t/shengdahsuang2/TCGA_RNA/CESC/function_analysis_fungi/temp/salmon/index -l A -p 8 --meta \
    -1 /data_200t/shengdahsuang2/TCGA_RNA/CESC/temp/kraken2_fungi/classified_micro/{1}_1_micro.fq \
    -2 /data_200t/shengdahsuang2/TCGA_RNA/CESC/temp/kraken2_fungi/classified_micro/{1}_2_micro.fq \
    -o temp/salmon/{1}.quant" \
::: `cat ../temp/ls.txt`

# 合并
mkdir -p result/salmon
salmon quantmerge \
--quants /data_200t/shengdahsuang2/TCGA_RNA/CESC/function_analysis_fungi/temp/salmon/*.quant \
-o /data_200t/shengdahsuang2/TCGA_RNA/CESC/function_analysis_fungi/result/salmon/gene.TPM

salmon quantmerge \
--quants temp/salmon/*.quant \
--column NumReads -o result/salmon/gene.count
sed -i '1 s/.quant//g' result/salmon/gene.*
  
  # 预览结果表格
head -n3 result/salmon/gene.*
  
## 3.3 功能基因注释
  # 输入数据：上一步预测的蛋白序列 result/NR/protein.fa
  # 中间结果：temp/eggnog/protein.emapper.seed_orthologs
  #           temp/eggnog/output.emapper.annotations
  #           temp/eggnog/output
  
  # COG定量表：result/eggnog/cogtab.count
  #            result/eggnog/cogtab.count.spf (用于STAMP)
  
  # KO定量表：result/eggnog/kotab.count
#           result/eggnog/kotab.count.spf  (用于STAMP)

# 这部分可以拓展到其它数据库
### 1.3.1 基因注释eggNOG(COG/KEGG/CAZy)
# https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2
# 记录软件版本
cd /data_200t/shengdahsuang2/TCGA_RNA/CESC/function_analysis_fungi
conda humann
emapper.py --version # 2.0.5-6
# diamond比对基因至eggNOG 5.0数据库, 1~9h，默认diamond 1e-3
mkdir -p temp/eggnog

time emapper.py --override \
-i result/NR/protein.fa --tax_scope all \
--cpu 12 -m diamond \
-o temp/eggnog/protein

# 添表头, 1列为ID，9列KO，16列CAZy，21列COG，22列描述
sed '1 i query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs' \
temp/eggnog/protein.emapper.annotations \
> temp/eggnog/output
grep -v '#' temp/eggnog/output > temp/eggnog/output1
head temp/eggnog/output1

# summarizeAbundance生成COG/KO/CAZy丰度汇总表

mkdir -p result/eggnog
# 显示帮助，需要Python3环境，可修改软件第一行指定python位置
# pip install pandas #下边脚本需要的module
# /home/worker019/TCGA_analysis/summarizeAbundance.py -h
# 汇总，7列COG按字母分隔，12列KO，13列KEGG通路按逗号分隔，原始值累加
#gene.TPM 标准化后的基因丰度表
/data_200t/shareAll/zhaolanlan/code/summarizeAbundance.py \
-i result/salmon/gene.TPM \
-m temp/eggnog/output \
-c '7,12,13' -s '*+,+,' -n raw \
-o result/eggnog/eggnog
sed -i 's/^ko://' result/eggnog/eggnog.KEGG_ko.raw.txt
# eggnog.CAZy.raw.txt  eggnog.COG.raw.txt  eggnog.KO.raw.txt