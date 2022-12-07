# 下载历史版本和现在版本的对应信息表

cd /data_200t/shengdahsuang2/TCGA_RNA/CESC/temp/match_id
cat -v gdc_manifest_20220405_061224.txt
sed -i "s/\r//g" gdc_manifest_20220405_061224.txt
cat -v gdc_manifest_20220405_061224.txt
a1=https://portal.gdc.cancer.gov/auth/api/history/
a2='?size=10000&attachment=true&format=JSON&fields=&filters=%7B%7D&pretty=true&downloadCookieKey=13173eb06&downloadCookiePath=%2F'

cat gdc_manifest_20220405_061224.txt | while read line
do
wget ${a1}${line}${a2}
done

# # R脚本
# setwd("/data_200t/shengdahsuang2/TCGA_RNA/CESC/temp/match_id")
# name <- read.table("gdc_manifest_20220405_061224.txt")
# match_id <- data.frame(current = rep(NA,309), old = NA)
# library(rjson)
# for (i in 1:309) {
#   a <- fromJSON(file = paste0(name[i,1], "?size=10000&attachment=true&format=JSON&fields=&filters={}&pretty=true&downloadCookieKey=13173eb06&downloadCookiePath=%2F"))
#   if(a[[1]][["version"]]=="1"){
#     match_id[i,1] <- a[[2]][["uuid"]]
#     match_id[i,2] <- a[[1]][["uuid"]]
#   }else{
#     match_id[i,1] <- a[[1]][["uuid"]]
#     match_id[i,2] <- a[[2]][["uuid"]]
#     print("false")
#   }
# }
# 
# write.csv(match_id, "match_id.csv", row.names = F)