# MYC signature calculation
rm(list=ls())
load('Reddy.RData')
load('Reddy_exp.RData')

myoutf1 = "BASE_signature.txt" # optput

#----------------------------
tmp1 = Reddy_pmid28985567_mutation # somatic mutation
tmp2 = Reddy_pmid28985567_CopyNumber # CNV
soma = cbind(tmp1, tmp2)
dim(soma)
soma = soma[,grepl('MYC',colnames(soma)),drop=F]
#---------------------------
data = Reddy_pmid28985567_exp # expression

info = Reddy_pmid28985567_clinical # clinical
se = c("age.at.diagnosis", "Gender", "IPI", "ABC.GCB.RNAseq")
info = info[,se]
colnames(info) = c("age", "gender", "stage", "type")
info$gender = as.character(info$gender)
info$age = as.numeric(info$age)
info$type = as.character(info$type)
xx = as.character(info$stage)
info$stage = xx

comSam = intersect(colnames(data), row.names(info))
data = data[, comSam]
info = info[comSam,]
raw.info = info
raw.data = data
#~~~~~~~~~~~~~~~~~~~~
tmp = matrix(0, nrow(data), ncol(soma))
row.names(tmp) = row.names(data)
colnames(tmp) = colnames(soma)
p1 = p2 = p3 = p4  =  tmp
## uni.noj --> p1:   univariate using each genomic feature without adjust by clinical variables
## uni.adj --> p2:   univariate using each genomic feature adjusted by clinical variables

for(s in 1:ncol(soma))
{
  mut = row.names(soma)[soma[,s]>0]
  wt =  row.names(soma)[soma[,s]==0]	
  mut = mut[mut%in%colnames(data)]
  wt = wt[wt%in%colnames(data)]
  xx = c(mut, wt)
  info = raw.info[xx, ]
  data = raw.data[, xx]
  mut = ifelse(row.names(info)%in%wt, 0, 1)
  info = cbind(mut, info)
  
  pval1 = beta1 = rep(0, nrow(data))
  pval2 = beta2 = rep(0, nrow(data))
  
  for(k in 1:nrow(data))
  {
    cat("\r\r\r", s, "-->", k)
    mytf = as.numeric(data[k,])
    xx = cbind(mytf, info)
    mylog <- lm(mytf ~ mut, data = xx)
    mylog = summary(mylog)
    pval1[k] =  mylog[["coefficients"]][2, 4]
    beta1[k] =  mylog[["coefficients"]][2, 1]
    
    mylog <- lm(mytf ~ mut + age + gender + stage + type, data = xx, family = "binomial")
    mylog = summary(mylog)
    pval2[k] =  mylog[["coefficients"]][2, 4]
    beta2[k] =  mylog[["coefficients"]][2, 1]
  }
  xx = -log10(pval1)
  xx = ifelse(beta1>0, xx, -xx)
  p1[,s] = xx
  xx = -log10(pval2)
  xx = ifelse(beta2>0, xx, -xx)
  p2[,s] = xx		
}


#--------------------------
res = p1
res1 = res2 = res
for(k in 1:ncol(res))
{
  tmp = res[,k]
  tmp1 = ifelse(tmp>0, tmp, 0)
  tmp1[tmp1>10] = 10
  res1[,k] = tmp1
  tmp2 = ifelse(tmp<0, -tmp, 0)
  tmp2[tmp2>10] = 10
  res2[,k] = tmp2
}
colnames(res1)= paste(colnames(res1), ".up", sep="")
colnames(res2)= paste(colnames(res2), ".dn", sep="")
myres = cbind(res1, res2)
row.names(myres) = row.names(res)
colnames(myres) = paste("uni.noj", colnames(myres), sep="__")
minv = min(myres)
maxv = max(myres)
myres = (myres-minv)/(maxv-minv)
prof1 = myres


res = p2
res1 = res2 = res
for(k in 1:ncol(res))
{
  tmp = res[,k]
  tmp1 = ifelse(tmp>0, tmp, 0)
  tmp1[tmp1>10] = 10
  res1[,k] = tmp1
  tmp2 = ifelse(tmp<0, -tmp, 0)
  tmp2[tmp2>10] = 10
  res2[,k] = tmp2
}
colnames(res1)= paste(colnames(res1), ".up", sep="")
colnames(res2)= paste(colnames(res2), ".dn", sep="")
myres = cbind(res1, res2)
row.names(myres) = row.names(res)
colnames(myres) = paste("uni.adj", colnames(myres), sep="__")
minv = min(myres)
maxv = max(myres)
myres = (myres-minv)/(maxv-minv)
prof2 = myres


profile = cbind(prof1, prof2)
se1 = grep("\\.up", colnames(profile))
se2 = grep("\\.dn", colnames(profile))
profile = profile[, c(se1, se2)]
dim(profile)	## 

write.table(profile, myoutf1, sep="\t", quote=F)






