#!/usr/bin/Rscript

# 1. maf ≤ 0.0004 AND avg_mapq ≥ 20
#    Note: The maf filter is dependent on using the maf data provided in one of the *data.zip files
# 2. (pct_double_split > 0.1 AND blast_hit != "no_hit") OR pct_double_split ≤ 0.1
# 3. ddg2p != "NA" AND sr_total ≥ 5 AND exonic = "True"
# 4. mom_sr < 2 AND dad_sr < 2

args = commandArgs(trailingOnly=TRUE)

input = args[1]
output = args[2]

annotated = read.table(input, header=T, sep="\t")
filtered = subset(annotated, 
                  (is.na(annotated$maf) | annotated$maf <= 0.0004) & annotated$avg_mapq >= 20 &
                  ((annotated$pct_double_split > 0.1 & annotated$blast_hit != "no_hit") | annotated$pct_double_split <= 1) &
                  !is.na(annotated$ddg2p) & annotated$sr_total >= 5 & annotated$exonic == "True" &
                  annotated$mom_sr < 2 & annotated$dad_sr < 2)

write.table(filtered, output, col.names=T, row.names=F, quote=F, sep="\t")

