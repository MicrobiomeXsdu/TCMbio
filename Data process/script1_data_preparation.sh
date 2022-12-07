cd /data_200t/shengdahsuang2/TCGA_RNA/
mkdir -p CESC &&cd CESC
mkdir temp
mkdir result

#######################################
# 从原始的BAM文件中提取未比对上的序列 #
#######################################

pw=/data_200t/micro_biome/TCGA_RNA/CESC/
cd $pw 
mkdir nohumanBAM

ls raw/*/*.bam > BAMfile.txt

# Extract unaligned sequences and output them as BAM files
cat BAMfile.txt | while read line
do
samtools view -b -h -f 4 ${line} > nohumanBAM/${line#*/*/}
samtools view nohumanBAM/${line#*/*/} | awk '$3=="*"' > nohumanBAM/${line:42:77}.sam
samtools view -bS nohumanBAM/${line:42:77}.sam > nohumanBAM/${line#*/*/}
rm nohumanBAM/${line:42:77}.sam #删除中间产生的sam文件（未压缩格式）节省空间
done

####################################
#       BAM文件转成fastq文件       #
####################################

cd /data_200t/micro_biome/TCGA_RNA/CESC
mkdir -p nohumanfastq
ls nohumanBAM > nohumanBAM.txt
# convert
cat nohumanBAM.txt | while read line
do
samtools sort -l 9 -n nohumanBAM/${line} -o nohumanBAM/sort_${line}
#convert
bedtools bamtofastq -i nohumanBAM/sort_${line} \
-fq nohumanfastq/${line}R1.fq \
-fq2 nohumanfastq/${line}R2.fq
done
# 不幸中断了
ls nohumanfastq > nohumanfastq.txt
# convert
cat nohumanBAM1.txt | while read line
do
samtools sort -l 9 -n nohumanBAM/${line} -o nohumanBAM/sort_${line}
#convert
bedtools bamtofastq -i nohumanBAM/sort_${line} \
-fq nohumanfastq/${line}R1.fq \
-fq2 nohumanfastq/${line}R2.fq
done

# 统计未必对上的序列数
echo "TCGA_CESC"  > readscount.txt
cat nohumanBAM.txt | while read line
do
echo ${line} >> readscount.txt
samtools view -c nohumanBAM/${line} >> readscount.txt
fpath=`find ./raw -name ${line}`
samtools view -c ${fpath} >> readscount.txt
done

#########################
# FastQC质量评估(可选)  #
#########################
cd /data_200t/shengdahsuang2/TCGA_RNA/CESC/
# 进入子目录，使输入文件没有目录更简洁
# fastqc每个文件一个线程，6个双端样本12个文件，设置6线程
data=/data_200t/micro_biome/TCGA_RNA/CESC/nohumanfastq
mkdir -p result/fastqc
time fastqc ${data}/*.fq -t 20 -o result/fastqc 

# 结果见result/fastqc目录，解读见[数据的质量控制软件——fastQC](https://mp.weixin.qq.com/s/MMfierO-8H2MVRkJKGCOVQ)

# 生成多样品报告比较
mkdir  result/qc
multiqc -d result/fastqc/ -o result/qc
# 查看右侧result/qc目录中multiqc_report.html，可交互式报告

# 利用trimomatic进行序列质控
mkdir -p temp/trimqc
cd ${data}
time parallel -j 4 --xapply \
  'trimmomatic PE -threads 5 -phred33 \
  -summary /data_200t/shengdahsuang2/TCGA_RNA/CESC/temp/trimqc/{1}sum.txt \
  {1} {2} \
  /data_200t/shengdahsuang2/TCGA_RNA/CESC/temp/trimqc/{1}_paired.fq.gz \
  /data_200t/shengdahsuang2/TCGA_RNA/CESC/temp/trimqc/{1}_unpaired.fq.gz \
  /data_200t/shengdahsuang2/TCGA_RNA/CESC/temp/trimqc/{2}_paired.fq.gz \
  /data_200t/shengdahsuang2/TCGA_RNA/CESC/temp/trimqc/{2}_unpaired.fq.gz \
  ILLUMINACLIP:/home/shengdashuang/miniconda2/share/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:35' \
  ::: *R1.fq ::: *R2.fq &> /data_200t/shengdahsuang2/TCGA_RNA/CESC/temp/trimqc/trimmomatic220324.log
# real    38m59.212s
# user    216m59.415s
# sys     12m12.587s


# 合并质控后的双端文件
cd /data_200t/micro_biome/TCGA_RNA/CESC/nohumanfastq/
ls *_gdc_realn_rehead.bamR1.fq > /data_200t/shengdahsuang2/TCGA_RNA/CESC/temp/ls.txt

cd /data_200t/shengdahsuang2/TCGA_RNA/CESC/temp
sed -i "s/_gdc_realn_rehead.bamR1.fq//g" ls.txt
mkdir -p concats
cat ls.txt | while read line
do
cat trimqc/${line}_gdc_realn_rehead.bamR1.fq_paired.fq.gz trimqc/${line}_gdc_realn_rehead.bamR2.fq_paired.fq.gz \
  > concats/${line}.fastq.gz
done

# 1.4 质控后质量再评估 (可选)
mkdir -p trimqc/fastqc
fastqc trimqc/*_gdc_realn_rehead.bamR?.fq_paired.fq.gz -t 6 -o trimqc/fastqc
multiqc -d trimqc/fastqc -o ../result/qc/trimqc/
# 整理bowtie2, trimmomatic, fastqc报告

##########
## 2.1 kraken2基于NCBI数据库注释reads层面
# 还可以进行contig、gene、bin层面的序列物种注释

### 2.1.1 多样本并行物种注释
cd /data_200t/shengdahsuang2/TCGA_RNA/CESC/
mkdir -p temp/kraken2
mkdir -p temp/kraken2/aligned
time parallel -j 3 'kraken2 --db /data_200t/shareAll/zhaolanlan/db/kraken2/210122/ \
  --paired --classified-out temp/kraken2/aligned/{1}#.fq temp/trimqc/{1}_gdc_realn_rehead.bamR1.fq_paired.fq.gz \
  temp/trimqc/{1}_gdc_realn_rehead.bamR2.fq_paired.fq.gz --threads 5 --use-names --report-zero-counts \
  --report temp/kraken2/{1}_report \
  --output temp/kraken2/{1}_output \
  --memory-mapping' \
  ::: `cat temp/ls.txt` &> temp/kraken2/kraken2.log

# 利用bracken将kraken结果转成物种丰度
mkdir temp/bracken_report
time parallel -j 3 'bracken -d /data_200t/shareAll/zhaolanlan/db/kraken2/210122/ \
  -i temp/kraken2/{1}_report -o temp/kraken2/{1}.bracken -w temp/bracken_report/{1}.bracken.report \
  -r 100 -l S' ::: `cat temp/ls.txt` > temp/bracken.log
# real    4m0.073s
# user    4m43.229s
# sys     0m15.505s
# 将Braken的report格式转换成--use-mpa-style格式
# 用Bracken的子程序
# conda install krakentools
time parallel -j 3 'kreport2mpa.py -r temp/bracken_report/{1}.bracken.report \
  -o temp/bracken_report/{1}.new.report' \
  ::: `cat temp/ls.txt`
# real    0m32.716s
# user    0m12.418s
# sys     0m4.778s

# # 统计微生物序列文件大小
# cd /data_200t/shengdahsuang2/TCGA_RNA/CESC/temp/
# ls -lsh micro_seqbind/ > micro_size.txt

### 2.1.2 汇总样品物种组成表
# Rscript kraken结果合并
mkdir result/composition
cd result/composition
sed "s/,/\t/g" kresult_redun.csv > kresult_redun.tsv

humann_renorm_table --input kresult_redun.tsv --units relab \
  --output kresult_redun_relab.tsv
  
# 提取比对至种水平的微生物
cd /data_200t/shengdahsuang2/TCGA_RNA/CESC/result/composition
echo "species,s2001,s2002,s2003,s2004,s2005,s2006,s2007,s2008,s2009,s2010,s2011,s2012,s2013,s2014,s2015,s2016,s2017,s2018,s2019,s2020,s2021,s2022,s2023,s2024,s2025,s2026,s2027,s2028,s2029,s2030,s2031,s2032,s2033,s2034,s2035,s2036,s2037,s2038,s2039,s2040,s2041,s2042,s2043,s2044,s2045,s2046,s2047,s2048,s2049,s2050,s2051,s2052,s2053,s2054,s2055,s2056,s2057,s2058,s2059,s2060,s2061,s2062,s2063,s2064,s2065,s2066,s2067,s2068,s2069,s2070,s2071,s2072,s2073,s2074,s2075,s2076,s2077,s2078,s2079,s2080,s2081,s2082,s2083,s2084,s2085,s2086,s2087,s2088,s2089,s2090,s2091,s2092,s2093,s2094,s2095,s2096,s2097,s2098,s2099,s2100,s2101,s2102,s2103,s2104,s2105,s2106,s2107,s2108,s2109,s2110,s2111,s2112,s2113,s2114,s2115,s2116,s2117,s2118,s2119,s2120,s2121,s2122,s2123,s2124,s2125,s2126,s2127,s2128,s2129,s2130,s2131,s2132,s2133,s2134,s2135,s2136,s2137,s2138,s2139,s2140,s2141,s2142,s2143,s2144,s2145,s2146,s2147,s2148,s2149,s2150,s2151,s2152,s2153,s2154,s2155,s2156,s2157,s2158,s2159,s2160,s2161,s2162,s2163,s2164,s2165,s2166,s2167,s2168,s2169,s2170,s2171,s2172,s2173,s2174,s2175,s2176,s2177,s2178,s2179,s2180,s2181,s2182,s2183,s2184,s2185,s2186,s2187,s2188,s2189,s2190,s2191,s2192,s2193,s2194,s2195,s2196,s2197,s2198,s2199,s2200,s2201,s2202,s2203,s2204,s2205,s2206,s2207,s2208,s2209,s2210,s2211,s2212,s2213,s2214,s2215,s2216,s2217,s2218,s2219,s2220,s2221,s2222,s2223,s2224,s2225,s2226,s2227,s2228,s2229,s2230,s2231,s2232,s2233,s2234,s2235,s2236,s2237,s2238,s2239,s2240,s2241,s2242,s2243,s2244,s2245,s2246,s2247,s2248,s2249,s2250,s2251,s2252,s2253,s2254,s2255,s2256,s2257,s2258,s2259,s2260,s2261,s2262,s2263,s2264,s2265,s2266,s2267,s2268,s2269,s2270,s2271,s2272,s2273,s2274,s2275,s2276,s2277,s2278,s2279,s2280,s2281,s2282,s2283,s2284,s2285,s2286,s2287,s2288,s2289,s2290,s2291,s2292,s2293,s2294,s2295,s2296,s2297,s2298,s2299,s2300,s2301,s2302,s2303,s2304,s2305,s2306,s2307,s2308,s2309" > kresult_s.csv
cat kresult.csv |grep -E 's__' >> kresult_s.csv

# 提取比对至属水平的微生物
cat kresult.csv |sed -e '/s__/d'   > kresult_m.csv
echo "genus,s2001,s2002,s2003,s2004,s2005,s2006,s2007,s2008,s2009,s2010,s2011,s2012,s2013,s2014,s2015,s2016,s2017,s2018,s2019,s2020,s2021,s2022,s2023,s2024,s2025,s2026,s2027,s2028,s2029,s2030,s2031,s2032,s2033,s2034,s2035,s2036,s2037,s2038,s2039,s2040,s2041,s2042,s2043,s2044,s2045,s2046,s2047,s2048,s2049,s2050,s2051,s2052,s2053,s2054,s2055,s2056,s2057,s2058,s2059,s2060,s2061,s2062,s2063,s2064,s2065,s2066,s2067,s2068,s2069,s2070,s2071,s2072,s2073,s2074,s2075,s2076,s2077,s2078,s2079,s2080,s2081,s2082,s2083,s2084,s2085,s2086,s2087,s2088,s2089,s2090,s2091,s2092,s2093,s2094,s2095,s2096,s2097,s2098,s2099,s2100,s2101,s2102,s2103,s2104,s2105,s2106,s2107,s2108,s2109,s2110,s2111,s2112,s2113,s2114,s2115,s2116,s2117,s2118,s2119,s2120,s2121,s2122,s2123,s2124,s2125,s2126,s2127,s2128,s2129,s2130,s2131,s2132,s2133,s2134,s2135,s2136,s2137,s2138,s2139,s2140,s2141,s2142,s2143,s2144,s2145,s2146,s2147,s2148,s2149,s2150,s2151,s2152,s2153,s2154,s2155,s2156,s2157,s2158,s2159,s2160,s2161,s2162,s2163,s2164,s2165,s2166,s2167,s2168,s2169,s2170,s2171,s2172,s2173,s2174,s2175,s2176,s2177,s2178,s2179,s2180,s2181,s2182,s2183,s2184,s2185,s2186,s2187,s2188,s2189,s2190,s2191,s2192,s2193,s2194,s2195,s2196,s2197,s2198,s2199,s2200,s2201,s2202,s2203,s2204,s2205,s2206,s2207,s2208,s2209,s2210,s2211,s2212,s2213,s2214,s2215,s2216,s2217,s2218,s2219,s2220,s2221,s2222,s2223,s2224,s2225,s2226,s2227,s2228,s2229,s2230,s2231,s2232,s2233,s2234,s2235,s2236,s2237,s2238,s2239,s2240,s2241,s2242,s2243,s2244,s2245,s2246,s2247,s2248,s2249,s2250,s2251,s2252,s2253,s2254,s2255,s2256,s2257,s2258,s2259,s2260,s2261,s2262,s2263,s2264,s2265,s2266,s2267,s2268,s2269,s2270,s2271,s2272,s2273,s2274,s2275,s2276,s2277,s2278,s2279,s2280,s2281,s2282,s2283,s2284,s2285,s2286,s2287,s2288,s2289,s2290,s2291,s2292,s2293,s2294,s2295,s2296,s2297,s2298,s2299,s2300,s2301,s2302,s2303,s2304,s2305,s2306,s2307,s2308,s2309" > kresult_g.csv
cat kresult_m.csv |grep -E 'g__' >> kresult_g.csv
# 提取比对至科水平的微生物
cat kresult.csv |sed -e '/g__/d'   > kresult_m.csv
echo "family,s2001,s2002,s2003,s2004,s2005,s2006,s2007,s2008,s2009,s2010,s2011,s2012,s2013,s2014,s2015,s2016,s2017,s2018,s2019,s2020,s2021,s2022,s2023,s2024,s2025,s2026,s2027,s2028,s2029,s2030,s2031,s2032,s2033,s2034,s2035,s2036,s2037,s2038,s2039,s2040,s2041,s2042,s2043,s2044,s2045,s2046,s2047,s2048,s2049,s2050,s2051,s2052,s2053,s2054,s2055,s2056,s2057,s2058,s2059,s2060,s2061,s2062,s2063,s2064,s2065,s2066,s2067,s2068,s2069,s2070,s2071,s2072,s2073,s2074,s2075,s2076,s2077,s2078,s2079,s2080,s2081,s2082,s2083,s2084,s2085,s2086,s2087,s2088,s2089,s2090,s2091,s2092,s2093,s2094,s2095,s2096,s2097,s2098,s2099,s2100,s2101,s2102,s2103,s2104,s2105,s2106,s2107,s2108,s2109,s2110,s2111,s2112,s2113,s2114,s2115,s2116,s2117,s2118,s2119,s2120,s2121,s2122,s2123,s2124,s2125,s2126,s2127,s2128,s2129,s2130,s2131,s2132,s2133,s2134,s2135,s2136,s2137,s2138,s2139,s2140,s2141,s2142,s2143,s2144,s2145,s2146,s2147,s2148,s2149,s2150,s2151,s2152,s2153,s2154,s2155,s2156,s2157,s2158,s2159,s2160,s2161,s2162,s2163,s2164,s2165,s2166,s2167,s2168,s2169,s2170,s2171,s2172,s2173,s2174,s2175,s2176,s2177,s2178,s2179,s2180,s2181,s2182,s2183,s2184,s2185,s2186,s2187,s2188,s2189,s2190,s2191,s2192,s2193,s2194,s2195,s2196,s2197,s2198,s2199,s2200,s2201,s2202,s2203,s2204,s2205,s2206,s2207,s2208,s2209,s2210,s2211,s2212,s2213,s2214,s2215,s2216,s2217,s2218,s2219,s2220,s2221,s2222,s2223,s2224,s2225,s2226,s2227,s2228,s2229,s2230,s2231,s2232,s2233,s2234,s2235,s2236,s2237,s2238,s2239,s2240,s2241,s2242,s2243,s2244,s2245,s2246,s2247,s2248,s2249,s2250,s2251,s2252,s2253,s2254,s2255,s2256,s2257,s2258,s2259,s2260,s2261,s2262,s2263,s2264,s2265,s2266,s2267,s2268,s2269,s2270,s2271,s2272,s2273,s2274,s2275,s2276,s2277,s2278,s2279,s2280,s2281,s2282,s2283,s2284,s2285,s2286,s2287,s2288,s2289,s2290,s2291,s2292,s2293,s2294,s2295,s2296,s2297,s2298,s2299,s2300,s2301,s2302,s2303,s2304,s2305,s2306,s2307,s2308,s2309" > kresult_f.csv
cat kresult_m.csv |grep -E 'f__' >> kresult_f.csv
# 提取比对至目水平的微生物
cat kresult.csv |sed -e '/f__/d'   > kresult_m.csv
echo "order,s2001,s2002,s2003,s2004,s2005,s2006,s2007,s2008,s2009,s2010,s2011,s2012,s2013,s2014,s2015,s2016,s2017,s2018,s2019,s2020,s2021,s2022,s2023,s2024,s2025,s2026,s2027,s2028,s2029,s2030,s2031,s2032,s2033,s2034,s2035,s2036,s2037,s2038,s2039,s2040,s2041,s2042,s2043,s2044,s2045,s2046,s2047,s2048,s2049,s2050,s2051,s2052,s2053,s2054,s2055,s2056,s2057,s2058,s2059,s2060,s2061,s2062,s2063,s2064,s2065,s2066,s2067,s2068,s2069,s2070,s2071,s2072,s2073,s2074,s2075,s2076,s2077,s2078,s2079,s2080,s2081,s2082,s2083,s2084,s2085,s2086,s2087,s2088,s2089,s2090,s2091,s2092,s2093,s2094,s2095,s2096,s2097,s2098,s2099,s2100,s2101,s2102,s2103,s2104,s2105,s2106,s2107,s2108,s2109,s2110,s2111,s2112,s2113,s2114,s2115,s2116,s2117,s2118,s2119,s2120,s2121,s2122,s2123,s2124,s2125,s2126,s2127,s2128,s2129,s2130,s2131,s2132,s2133,s2134,s2135,s2136,s2137,s2138,s2139,s2140,s2141,s2142,s2143,s2144,s2145,s2146,s2147,s2148,s2149,s2150,s2151,s2152,s2153,s2154,s2155,s2156,s2157,s2158,s2159,s2160,s2161,s2162,s2163,s2164,s2165,s2166,s2167,s2168,s2169,s2170,s2171,s2172,s2173,s2174,s2175,s2176,s2177,s2178,s2179,s2180,s2181,s2182,s2183,s2184,s2185,s2186,s2187,s2188,s2189,s2190,s2191,s2192,s2193,s2194,s2195,s2196,s2197,s2198,s2199,s2200,s2201,s2202,s2203,s2204,s2205,s2206,s2207,s2208,s2209,s2210,s2211,s2212,s2213,s2214,s2215,s2216,s2217,s2218,s2219,s2220,s2221,s2222,s2223,s2224,s2225,s2226,s2227,s2228,s2229,s2230,s2231,s2232,s2233,s2234,s2235,s2236,s2237,s2238,s2239,s2240,s2241,s2242,s2243,s2244,s2245,s2246,s2247,s2248,s2249,s2250,s2251,s2252,s2253,s2254,s2255,s2256,s2257,s2258,s2259,s2260,s2261,s2262,s2263,s2264,s2265,s2266,s2267,s2268,s2269,s2270,s2271,s2272,s2273,s2274,s2275,s2276,s2277,s2278,s2279,s2280,s2281,s2282,s2283,s2284,s2285,s2286,s2287,s2288,s2289,s2290,s2291,s2292,s2293,s2294,s2295,s2296,s2297,s2298,s2299,s2300,s2301,s2302,s2303,s2304,s2305,s2306,s2307,s2308,s2309" > kresult_o.csv
cat kresult_m.csv |grep -E 'o__' >> kresult_o.csv
# 提取比对至纲水平的微生物
cat kresult.csv |sed -e '/o__/d'   > kresult_m.csv
echo "family,s2001,s2002,s2003,s2004,s2005,s2006,s2007,s2008,s2009,s2010,s2011,s2012,s2013,s2014,s2015,s2016,s2017,s2018,s2019,s2020,s2021,s2022,s2023,s2024,s2025,s2026,s2027,s2028,s2029,s2030,s2031,s2032,s2033,s2034,s2035,s2036,s2037,s2038,s2039,s2040,s2041,s2042,s2043,s2044,s2045,s2046,s2047,s2048,s2049,s2050,s2051,s2052,s2053,s2054,s2055,s2056,s2057,s2058,s2059,s2060,s2061,s2062,s2063,s2064,s2065,s2066,s2067,s2068,s2069,s2070,s2071,s2072,s2073,s2074,s2075,s2076,s2077,s2078,s2079,s2080,s2081,s2082,s2083,s2084,s2085,s2086,s2087,s2088,s2089,s2090,s2091,s2092,s2093,s2094,s2095,s2096,s2097,s2098,s2099,s2100,s2101,s2102,s2103,s2104,s2105,s2106,s2107,s2108,s2109,s2110,s2111,s2112,s2113,s2114,s2115,s2116,s2117,s2118,s2119,s2120,s2121,s2122,s2123,s2124,s2125,s2126,s2127,s2128,s2129,s2130,s2131,s2132,s2133,s2134,s2135,s2136,s2137,s2138,s2139,s2140,s2141,s2142,s2143,s2144,s2145,s2146,s2147,s2148,s2149,s2150,s2151,s2152,s2153,s2154,s2155,s2156,s2157,s2158,s2159,s2160,s2161,s2162,s2163,s2164,s2165,s2166,s2167,s2168,s2169,s2170,s2171,s2172,s2173,s2174,s2175,s2176,s2177,s2178,s2179,s2180,s2181,s2182,s2183,s2184,s2185,s2186,s2187,s2188,s2189,s2190,s2191,s2192,s2193,s2194,s2195,s2196,s2197,s2198,s2199,s2200,s2201,s2202,s2203,s2204,s2205,s2206,s2207,s2208,s2209,s2210,s2211,s2212,s2213,s2214,s2215,s2216,s2217,s2218,s2219,s2220,s2221,s2222,s2223,s2224,s2225,s2226,s2227,s2228,s2229,s2230,s2231,s2232,s2233,s2234,s2235,s2236,s2237,s2238,s2239,s2240,s2241,s2242,s2243,s2244,s2245,s2246,s2247,s2248,s2249,s2250,s2251,s2252,s2253,s2254,s2255,s2256,s2257,s2258,s2259,s2260,s2261,s2262,s2263,s2264,s2265,s2266,s2267,s2268,s2269,s2270,s2271,s2272,s2273,s2274,s2275,s2276,s2277,s2278,s2279,s2280,s2281,s2282,s2283,s2284,s2285,s2286,s2287,s2288,s2289,s2290,s2291,s2292,s2293,s2294,s2295,s2296,s2297,s2298,s2299,s2300,s2301,s2302,s2303,s2304,s2305,s2306,s2307,s2308,s2309" > kresult_c.csv
cat kresult_m.csv |grep -E 'c__' >> kresult_c.csv
# 提取比对至门水平的微生物
cat kresult.csv |sed -e '/c__/d'   > kresult_m.csv
echo "family,s2001,s2002,s2003,s2004,s2005,s2006,s2007,s2008,s2009,s2010,s2011,s2012,s2013,s2014,s2015,s2016,s2017,s2018,s2019,s2020,s2021,s2022,s2023,s2024,s2025,s2026,s2027,s2028,s2029,s2030,s2031,s2032,s2033,s2034,s2035,s2036,s2037,s2038,s2039,s2040,s2041,s2042,s2043,s2044,s2045,s2046,s2047,s2048,s2049,s2050,s2051,s2052,s2053,s2054,s2055,s2056,s2057,s2058,s2059,s2060,s2061,s2062,s2063,s2064,s2065,s2066,s2067,s2068,s2069,s2070,s2071,s2072,s2073,s2074,s2075,s2076,s2077,s2078,s2079,s2080,s2081,s2082,s2083,s2084,s2085,s2086,s2087,s2088,s2089,s2090,s2091,s2092,s2093,s2094,s2095,s2096,s2097,s2098,s2099,s2100,s2101,s2102,s2103,s2104,s2105,s2106,s2107,s2108,s2109,s2110,s2111,s2112,s2113,s2114,s2115,s2116,s2117,s2118,s2119,s2120,s2121,s2122,s2123,s2124,s2125,s2126,s2127,s2128,s2129,s2130,s2131,s2132,s2133,s2134,s2135,s2136,s2137,s2138,s2139,s2140,s2141,s2142,s2143,s2144,s2145,s2146,s2147,s2148,s2149,s2150,s2151,s2152,s2153,s2154,s2155,s2156,s2157,s2158,s2159,s2160,s2161,s2162,s2163,s2164,s2165,s2166,s2167,s2168,s2169,s2170,s2171,s2172,s2173,s2174,s2175,s2176,s2177,s2178,s2179,s2180,s2181,s2182,s2183,s2184,s2185,s2186,s2187,s2188,s2189,s2190,s2191,s2192,s2193,s2194,s2195,s2196,s2197,s2198,s2199,s2200,s2201,s2202,s2203,s2204,s2205,s2206,s2207,s2208,s2209,s2210,s2211,s2212,s2213,s2214,s2215,s2216,s2217,s2218,s2219,s2220,s2221,s2222,s2223,s2224,s2225,s2226,s2227,s2228,s2229,s2230,s2231,s2232,s2233,s2234,s2235,s2236,s2237,s2238,s2239,s2240,s2241,s2242,s2243,s2244,s2245,s2246,s2247,s2248,s2249,s2250,s2251,s2252,s2253,s2254,s2255,s2256,s2257,s2258,s2259,s2260,s2261,s2262,s2263,s2264,s2265,s2266,s2267,s2268,s2269,s2270,s2271,s2272,s2273,s2274,s2275,s2276,s2277,s2278,s2279,s2280,s2281,s2282,s2283,s2284,s2285,s2286,s2287,s2288,s2289,s2290,s2291,s2292,s2293,s2294,s2295,s2296,s2297,s2298,s2299,s2300,s2301,s2302,s2303,s2304,s2305,s2306,s2307,s2308,s2309" > kresult_p.csv
cat kresult_m.csv |grep -E 'p__' >> kresult_p.csv


# 3.1.4 提取比对至微生物的序列
cd /data_200t/shengdahsuang2/TCGA_RNA/CESC/
mkdir -p temp/kraken2/classified_micro
cat temp/ls.txt | while read line
do
cat temp/kraken2/aligned/${line}_1.fq | sed '/kraken:taxid|131567/,+3d' | sed '/kraken:taxid|9606/,+3d' > temp/kraken2/classified_micro/${line}_micro_1.fq
cat temp/kraken2/aligned/${line}_2.fq | sed '/kraken:taxid|131567/,+3d' | sed '/kraken:taxid|9606/,+3d' > temp/kraken2/classified_micro/${line}_micro_2.fq
done
