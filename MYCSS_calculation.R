rm(list=ls())
load('MYCsignature.RData')
load('Reddy_exp.RData')

myoutf = "MYCSS.txt" # optput

mywt = SignatureMYC
source("base5.R")
data = Reddy_pmid28985567_exp

xx = base5(data, mywt, perm=1000, myoutf, median.norm=T)


ES = read.table(myoutf, sep="\t", header=T, row.names=1)
cnum = ncol(ES)/2
ES = ES[, 1:cnum]
tmp = colnames(ES)
tmp = gsub("\\.ES", "", tmp)
colnames(ES) = tmp
ES$MYCSS = ES$MYC_up - ES$MYC_dn  # MYCSS