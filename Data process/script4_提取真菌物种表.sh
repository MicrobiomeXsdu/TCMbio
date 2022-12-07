### 2.1.2 汇总样品物种组成表
# Rscript kraken结果合并
cd /data_200t/shengdahsuang2/TCGA_RNA/CESC/result/composition_fungi

# 提取比对至种水平的微生物
cut -d  ","  -f 2 CESC_kresult_name.csv > CESC_kresult_s.csv
sed -i ":a;N;s/\n/,/g;$!ba" CESC_kresult_s.csv
cat CESC_kresult.csv |grep -E 's__' >> CESC_kresult_s.csv

# 提取比对至属水平的微生物
cat CESC_kresult.csv |sed -e '/s__/d'   > kresult_m.csv

cut -d  ","  -f 2 CESC_kresult_name.csv > CESC_kresult_g.csv
sed -i ":a;N;s/\n/,/g;$!ba" CESC_kresult_g.csv
cat kresult_m.csv |grep -E 'g__' >> CESC_kresult_g.csv


# 提取比对至科水平的微生物
cat CESC_kresult.csv |sed -e '/s__/d'   > kresult_m.csv
sed -i -e '/g__/d' kresult_m.csv

cut -d  ","  -f 2 CESC_kresult_name.csv > CESC_kresult_f.csv
sed -i ":a;N;s/\n/,/g;$!ba" CESC_kresult_f.csv
cat kresult_m.csv |grep -E 'f__' >> CESC_kresult_f.csv

# 提取比对至目水平的微生物
cat CESC_kresult.csv |sed -e '/s__/d'   > kresult_m.csv
sed -i -e '/g__/d' kresult_m.csv
sed -i -e '/f__/d' kresult_m.csv

cut -d  ","  -f 2 CESC_kresult_name.csv > CESC_kresult_o.csv
sed -i ":a;N;s/\n/,/g;$!ba" CESC_kresult_o.csv
cat kresult_m.csv |grep -E 'o__' >> CESC_kresult_o.csv

# 提取比对至纲水平的微生物
cat CESC_kresult.csv |sed -e '/s__/d'   > kresult_m.csv
sed -i -e '/g__/d' kresult_m.csv
sed -i -e '/f__/d' kresult_m.csv
sed -i -e '/o__/d' kresult_m.csv

cut -d  ","  -f 2 CESC_kresult_name.csv > CESC_kresult_c.csv
sed -i ":a;N;s/\n/,/g;$!ba" CESC_kresult_c.csv
cat kresult_m.csv |grep -E 'c__' >> CESC_kresult_c.csv

# 提取比对至门水平的微生物
cat CESC_kresult.csv |sed -e '/s__/d'   > kresult_m.csv
sed -i -e '/g__/d' kresult_m.csv
sed -i -e '/f__/d' kresult_m.csv
sed -i -e '/o__/d' kresult_m.csv
sed -i -e '/c__/d' kresult_m.csv

cut -d  ","  -f 2 CESC_kresult_name.csv > CESC_kresult_p.csv
sed -i ":a;N;s/\n/,/g;$!ba" CESC_kresult_p.csv
cat kresult_m.csv |grep -E 'p__' >> CESC_kresult_p.csv

