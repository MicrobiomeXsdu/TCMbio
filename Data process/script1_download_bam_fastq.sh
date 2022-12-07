################### Note ###########################
# 整体分四步：下载数据、过滤non-human reads、
# 在进行每一步时一定要使用pwd看清楚自己的位置、创建好相应的文件夹
# 创建文件夹请以TCGA项目名称，例如食管癌为ECAD
# 以下脚本以食管癌为例,并在28的服务器上运行。
# 用到了两个文件夹：
## 第一个是存储数据和过滤的文件夹，/mnt/ba76940c-bfd0-42f0-a133-112be8e744c9/TCGA_RNA
## 第二个是存储结果的文件夹， /mnt/ba76940c-bfd0-42f0-a133-112be8e744c9/TCGA_RESULT

####### 一. 下载前准备###############
#1.进入数据存储的位置。
cd /mnt/ba76940c-bfd0-42f0-a133-112be8e744c9/TCGA_RNA
#2.创建各自负责下载的癌症的文件夹,并进入。
## ***ECAD是食管癌的，到时候需要改。
mkdir -p ECAD
cd ECAD
#3. 利用filezilla上传manefest文件和token（token是可以下载数据的权限的文件）到所在文件夹下
## ***因为我们可以并行下载数据，所以我们可以写多个manefest，同时下载
### 即将原始的manefest文件拆成n个manefest（manefest1，2，3...）
### 保证每个manifest有50-100个需要下载的样本的名称
#4. 创建存储下载测序数据的文件夹并进入。不需要改文件名
mkdir -p rna_seq
cd rna_seq


############## 二.利用gdc-client从TCGA中下载数据##########
## 切记要进入各自存储下载的测序数据的文件夹下
#1. 创建下载的数据存储的文件夹：rna_seq/
cd /mnt/ba76940c-bfd0-42f0-a133-112be8e744c9/TCGA_RNA/ECAD/rna_seq

#2. 利用gdc-client下载rna测序数据到/rna_seq、文件夹下
## *** 其中-m后面接的是manefest存储的位置和文件名
## *** -t后面接的是token存储的位置和文件名
## *** 切记修改文件名
/mnt/ba76940c-bfd0-42f0-a133-112be8e744c9/apps/gdc/gdc-client download -m /mnt/ba76940c-bfd0-42f0-a133-112be8e744c9/TCGA_RNA/manifest.txt -t /mnt/ba76940c-bfd0-42f0-a133-112be8e744c9/TCGA_RNA/token.txt


############## 三. 从原始的BAM文件中提取未比对上的序列 ############
####备注:用到了samtools 和 bedtools，这两个插件。。
####samtools:https://www.cnblogs.com/emanlee/p/4316581.html
# 进入各自负责的癌症的文件夹下.
pw=/mnt/ba76940c-bfd0-42f0-a133-112be8e744c9/TCGA_RNA/ESCA/
cd $pw 
# 1. 创建存储 过滤掉人的测序文件 的文件夹nohumanBAM。文件夹名称不需要改。
mkdir nohumanBAM
# 2. 写出测序文件夹下 以bam结尾的 文件名，并将其存储到BAMfile.txt里面。
ls rna_seq/*/*.bam > BAMfile.txt
# 3. 依次提取 未比对上human的 reads并写到 nohumanBAM 文件夹下。不需要改动。
cat BAMfile.txt | while read line
do
samtools view -b -h -f 4 ${line} > nohumanBAM/${line#*/*/}
samtools view nohumanBAM/${line#*/*/} | awk '$3=="*"' > nohumanBAM/${line:45:53}.sam
samtools view -bS nohumanBAM/${line:45:53}.sam > nohumanBAM/${line#*/*/}
rm nohumanBAM/${line:45:53}.sam #删除中间产生的sam文件（未压缩格式）节省空间
done

######## 四. 对BAM文件进行排序并转成fastq文件 #################
# 进入各自负责的癌症文件夹下
cd /mnt/ba76940c-bfd0-42f0-a133-112be8e744c9/TCGA_RNA/ESCA
# 1. 写出nohumanBAM文件夹下的所有文件名，并将其存储到nohumanBAM.txt里面。
ls nohumanBAM > nohumanBAM.txt

# 2. 删除windows处理文件时带来的^M，并对bam文件进行排序
cat nohumanBAM.txt| tr -d "\r" > nohumanBAM.1.txt 

# 3. 创建文件夹nohumanfastq，以存储下一步bam转fastq的文件，但不需要进入该文件夹。
## bedtools将将nohumanBAM下的bam测序文件 转为 R1和R2双端的fastq测序文件，并存储到nohumanfastq文件夹下
mkdir nohumanfastq
cat nohumanBAM.1.txt | while read line
do
samtools sort -l 9 -n nohumanBAM/${line} -o nohumanBAM/sort_${line}
bedtools bamtofastq -i nohumanBAM/sort_${line} \
-fq nohumanfastq/${line}R1.fq \
-fq2 nohumanfastq/${line}R2.fq
done

# 4. 删除未排序的bam文件，保留排序后的bam文件以节省空间
cat nohumanBAM.txt | while read line
do
rm nohumanBAM/${line}
done

# 5. 统计未必对上的序列数，并将结果输出到另一个存储输出结果的文件夹下。
## ***先进入存储result的文件夹下，创捷对应癌症的文件夹
cd /mnt/ba76940c-bfd0-42f0-a133-112be8e744c9/TCGA_RESULT/
mkdir ESCA
## ***切记创建完成后，回到存储下载数据的文件夹下
cd /mnt/ba76940c-bfd0-42f0-a133-112be8e744c9/TCGA_RNA/ESCA
cat nohumanBAM.txt | while read line
do
echo ${line} >> /mnt/ba76940c-bfd0-42f0-a133-112be8e744c9/TCGA_RESULT/ESCA/readscount.txt
samtools view -c nohumanBAM/sort_${line} >> /mnt/ba76940c-bfd0-42f0-a133-112be8e744c9/TCGA_RESULT/ESCA/readscount.txt
fpath=`find . -name ${line}`
samtools view -c ${fpath} >> /mnt/ba76940c-bfd0-42f0-a133-112be8e744c9/TCGA_RESULT/ESCA/readscount.txt
done



