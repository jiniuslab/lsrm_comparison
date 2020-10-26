rm(list = ls())

tryCatch( # set wording directory to current source file folder
  {setwd(getSrcDirectory()[1])}
  ,error=function(e){setwd(dirname(rstudioapi::getActiveDocumentContext()$path))}
)
# load R library for files
library(MCMCpack)
library(coda)
library(xtable)
library(kknn)
library(tictoc)

library('gplots')
library(scales)
#install.packages('fpc')
library(fpc)
#install.packages('pheatmap')
library(pheatmap)
library(gridExtra)
library(ggplot2)
library(Rfast)
library(ks)
getwd()
# function define ---------------------------------------------
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
Mode <- function(x,inp_round=0) {
  if(inp_round==0){
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  else{
    x=round(x,inp_round)
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
}
colModes <-function(mat, inp_round = 1) {
  if ((max(mat) < 1)&(inp_round=1)) {
    inp_round = 2
  }
  mat = round(mat, inp_round)
  return(apply(X = mat, MARGIN = 2, FUN = function(x) Mode(x)))
}

# Calculate numerical integral of duplicated pdf mass for detecting different distribution
num_integral=function(inp_samp1,inp_samp2,inp_round=1,inp_lb=0.1){
  inp_eval=seq(min(min(inp_samp1),min(inp_samp2)),max(max(inp_samp1),max(inp_samp2)),length.out = 400)
  samp1_kde=kde(inp_samp1,eval.points = inp_eval)
  samp2_kde=kde(inp_samp2,eval.points = inp_eval)
  cdf1=cumsum(samp1_kde$estimate)/sum(samp1_kde$estimate)
  cdf2=cumsum(samp2_kde$estimate)/sum(samp2_kde$estimate)    
  
  min1=min(samp1_kde$eval.points[which((0.025<cdf1)&(cdf1<0.975))]);max1=max(samp1_kde$eval.points[which((0.025<cdf1)&(cdf1<0.975))])
  min2=min(samp2_kde$eval.points[which((0.025<cdf2)&(cdf2<0.975))]);max2=max(samp2_kde$eval.points[which((0.025<cdf2)&(cdf2<0.975))])
  mmin=max(min1,min2);mmax=min(max1,max2)
  
  pdf1=samp1_kde$estimate/sum(samp1_kde$estimate)
  pdf2=samp2_kde$estimate/sum(samp2_kde$estimate)
  
  tmp_data=rbind(pdf1[(mmin<samp1_kde$eval.points)&(samp1_kde$eval.points<mmax)],pdf2[(mmin<samp2_kde$eval.points)&(samp2_kde$eval.points<mmax)])
  return(sum(colMins(tmp_data,value = T)))
}

num_integral_v2=function(inp_samp1,inp_samp2,inp_round=1,inp_lb=0.1){
  inp_eval=seq(min(min(inp_samp1),min(inp_samp2)),max(max(inp_samp1),max(inp_samp2)),length.out = 400)
  samp1_kde=kde(inp_samp1,eval.points = inp_eval)
  samp2_kde=kde(inp_samp2,eval.points = inp_eval)
  
  mmin=min(inp_eval);mmax=max(inp_eval)
  
  pdf1=samp1_kde$estimate/sum(samp1_kde$estimate)
  pdf2=samp2_kde$estimate/sum(samp2_kde$estimate)
  
  tmp_data=rbind(pdf1[(mmin<samp1_kde$eval.points)&(samp1_kde$eval.points<mmax)],pdf2[(mmin<samp2_kde$eval.points)&(samp2_kde$eval.points<mmax)])
  return(sum(colMins(tmp_data,value = T)))
}

#calculate KL-divergence
KL_divergence_v2=function(inp_samp1,inp_samp2){
  mmin=min(inp_samp1,inp_samp2)
  mmax=max(inp_samp1,inp_samp2)
  
  kde1=kde(inp_samp1,eval.points = seq(mmin,mmax,length.out = 400))
  kde2=kde(inp_samp2,eval.points = seq(mmin,mmax,length.out = 400))
  P=kde1$estimate/sum(kde1$estimate)+1e-6
  Q=kde2$estimate/sum(kde2$estimate)+1e-6
  
  return(sum(P*log(P/Q)))
}

KL_divergence=function(inp_samp1,inp_samp2,inp_round=1,inp_lb=0.1){
  stopifnot(length(inp_samp1)==length(inp_samp2))
  inp_samp1=round(inp_samp1,inp_round)
  inp_samp2=round(inp_samp2,inp_round)
  
  dist2=as.data.frame(table(inp_samp2))
  colnames(dist2)=c('axis','Freq')
  
  dist1=as.data.frame(table(inp_samp1))
  colnames(dist1)=c('axis','Freq')
  
  merged_dist=merge(x = dist1,  y = dist2,  by = 'axis', all = TRUE)
  merged_dist$axis=as.numeric(as.character(merged_dist$axis))
  merged_dist=merged_dist[order(merged_dist$axis),]
  merged_dist[is.na(merged_dist)]=0
  
  dist1=merged_dist$Freq.x
  dist2=merged_dist$Freq.y
  
  P=dist1+inp_lb
  Q=dist2+inp_lb
  P=P/sum(P)# We need to change P&Q sum 1 to be pdf
  Q=Q/sum(Q)# We need to change P&Q sum 1 to be pdf
  
  return(sum(P*log(P/Q)))
}

Jaccard_sim=function(inp.mat){
  stopifnot(max(inp.mat)<=1)
  x.y=as.matrix(t(inp.mat))%*%as.matrix(inp.mat)
  x.yn=as.matrix(t(inp.mat))%*%as.matrix(1-inp.mat)
  xn.y=as.matrix(t(1-inp.mat))%*%as.matrix(inp.mat)
  
  Jaccard_sim.mat=(x.y/(x.y+x.yn+xn.y))
  return(Jaccard_sim.mat)
}

# Load fitted results------------------------------------------------------------
# setup parameter values
ndim = 2; nmcmc = 5000

raw_data = read.csv('../Item_Meta.txt',header = T,skip = 1,sep = '\t')
raw_data$Variable=sapply(raw_data$Variable,FUN = function(x) substr(x,6,nchar(as.character(x))))

item_name=as.character(as.matrix(read.table('../item_name.txt')))
item_name=sapply(item_name,FUN = function(x) substr(x,6,nchar(as.character(x))))


YSR_data = as.matrix(read.table("../YSR/data/item.txt"))
YSR_data[YSR_data==999]=NA

CBCL_data = as.matrix(read.table("../CBCL/data/item.txt"))
CBCL_data[CBCL_data==999]=NA

# defining group index---------------------------------------------------
AB = c('03', '16', '19', '20', '21', '22', '23', '37', '57', '68', '86', '87', '89', '94', '95', '97', '104')
AD = c('14', '29', '30', '31', '32', '33', '35', '45', '50', '52', '71', '91', '112')
AP = c('01', '04', '08', '10', '13', '17', '41', '61', '78')
RBB = c('02', '26', '28', '39', '43', '63', '67', '72', '81', '82', '90', '96', '99', '101', '105')
SP = c('11', '12', '25', '27', '34', '36', '38', '48', '62', '64', '79')
TP = c('09', '18', '40', '46', '58', '66', '70', '76', '83', '84', '85', '100')
WD = c('05', '42', '65', '69', '75', '102', '103', '111')
SC = c('47', '51', '54', '56A', '56B', '56C', '56D', '56E', '56F', '56G')
# AB=as.character(sapply(X=AB, FUN = function(x) paste('CBCL_',x,sep = '')))
# AD=as.character(sapply(X=AD, FUN = function(x) paste('CBCL_',x,sep = '')))
# AP=as.character(sapply(X=AP, FUN = function(x) paste('CBCL_',x,sep = '')))
# RBB=as.character(sapply(X=RBB, FUN = function(x) paste('CBCL_',x,sep = '')))
# SP=as.character(sapply(X=SP, FUN = function(x) paste('CBCL_',x,sep = '')))
# TP=as.character(sapply(X=TP, FUN = function(x) paste('CBCL_',x,sep = '')))            
# WD=as.character(sapply(X=WD, FUN = function(x) paste('CBCL_',x,sep = '')))
# SC=as.character(sapply(X=SC, FUN = function(x) paste('CBCL_',x,sep = '')))
data_list=list(AB,AD,AP,RBB,SP,TP,WD,SC)
name_list=c('AB','AD','AP','RBB','SP','TP','WD','SC')

group_name=as.vector(item_name)
group_sorted_idx=c(which(group_name=='RBB'),which(group_name=='AB'),which(group_name=='AP'),which(group_name=='SP'),which(group_name=='TP'),which(group_name=='WD'),which(group_name=='SC'),which(group_name=='AD'))

#load CBCL result-----------------------------------
tic()
# load result files
{CBCL_impute = as.matrix(read.table("../CBCL/result/impute.txt"))
  CBCL_beta = as.matrix(read.table("../CBCL/result/beta.txt"))
  CBCL_theta = as.matrix(read.table("../CBCL/result/theta.txt"))
  nsample = ncol(CBCL_theta)
  nitem = ncol(CBCL_beta)
  CBCL_z = array(NA,dim=c(nmcmc,nsample,ndim))
  CBCL_w = array(NA,dim=c(nmcmc,nitem,ndim))
  for(j in 1:ndim){
    fopen = paste("../CBCL/result/z",j,".txt",sep="")
    CBCL_z[,,j] = as.matrix(read.table(fopen))
  }
  for(j in 1:ndim){
    fopen = paste("../CBCL/result/w",j,".txt",sep="")
    CBCL_w[,,j] = as.matrix(read.table(fopen))
  }
  CBCL_var_beta = scan("../CBCL/result/var_beta.txt")
  CBCL_var_theta = scan("../CBCL/result/var_theta.txt")
  CBCL_mle = scan("../CBCL/result/mle.txt")}
toc()

tic()
# preparation of Procrustes matching
{CBCL_max.address = which.max(CBCL_mle)
  CBCL_w.star = CBCL_w[CBCL_max.address,,]
  CBCL_z.star = CBCL_z[CBCL_max.address,,]
  CBCL_w.proc = array(0,dim=c(nmcmc,nitem,ndim))
  CBCL_z.proc = array(0,dim=c(nmcmc,nsample,ndim))}

# Procrustes matching for latent spaces Z (respondent) and W (item)
for(iter in 1:nmcmc){
  CBCL_z.iter = CBCL_z[iter,,]
  CBCL_w.iter = CBCL_w[iter,,]
  if(iter != CBCL_max.address){
    CBCL_w.proc[iter,,] = procrustes(CBCL_w.iter,CBCL_w.star)$X.new
    CBCL_z.proc[iter,,] = procrustes(CBCL_z.iter,CBCL_z.star)$X.new
  }
  else{
    CBCL_w.proc[iter,,] = CBCL_w.iter
    CBCL_z.proc[iter,,] = CBCL_z.iter
  }
}
toc()
#Load YSR result----------------------------------
tic()
# load result files
{YSR_impute = as.matrix(read.table("../YSR/result/impute.txt"))
  YSR_beta = as.matrix(read.table("../YSR/result/beta.txt"))
  YSR_theta = as.matrix(read.table("../YSR/result/theta.txt"))
  nsample = ncol(YSR_theta)
  nitem = ncol(YSR_beta)
  YSR_z = array(NA,dim=c(nmcmc,nsample,ndim))
  YSR_w = array(NA,dim=c(nmcmc,nitem,ndim))
  for(j in 1:ndim){
    fopen = paste("../YSR/result/z",j,".txt",sep="")
    YSR_z[,,j] = as.matrix(read.table(fopen))
  }
  for(j in 1:ndim){
    fopen = paste("../YSR/result/w",j,".txt",sep="")
    YSR_w[,,j] = as.matrix(read.table(fopen))
  }
  YSR_var_beta = scan("../YSR/result/var_beta.txt")
  YSR_var_theta = scan("../YSR/result/var_theta.txt")
  YSR_mle = scan("../YSR/result/mle.txt")}
toc()

tic()
# preparation of Procrustes matching
{YSR_max.address = which.max(YSR_mle)
  YSR_w.star = YSR_w[YSR_max.address,,]
  YSR_z.star = YSR_z[YSR_max.address,,]
  YSR_w.proc = array(0,dim=c(nmcmc,nitem,ndim))
  YSR_z.proc = array(0,dim=c(nmcmc,nsample,ndim))}

# Procrustes matching for latent spaces Z (respondent) and W (item)
for(iter in 1:nmcmc){
  YSR_z.iter = YSR_z[iter,,]
  YSR_w.iter = YSR_w[iter,,]
  if(iter != YSR_max.address){
    YSR_w.proc[iter,,] = procrustes(YSR_w.iter,YSR_w.star)$X.new
    YSR_z.proc[iter,,] = procrustes(YSR_z.iter,YSR_z.star)$X.new
  }
  else{
    YSR_w.proc[iter,,] = YSR_w.iter
    YSR_z.proc[iter,,] = YSR_z.iter
  }
}
toc()

#Calculate point estimate----------------------------------------------------

tic()
YSR_w.est = matrix(NA,nitem,ndim)
YSR_z.est = matrix(NA,nsample,ndim)
for(i in 1:nitem){
  for(j in 1:ndim){
    YSR_w.est[i,j] = mean(YSR_w.proc[,i,j])
  }
}
for(k in 1:nsample){
  for(j in 1:ndim){
    YSR_z.est[k,j] = mean(YSR_z.proc[,k,j])
  }
}
toc()

YSR_beta.est=colMeans(YSR_beta)
YSR_theta.est=colMeans(YSR_theta)


tic()
CBCL_w.est = matrix(NA,nitem,ndim)
CBCL_z.est = matrix(NA,nsample,ndim)
for(i in 1:nitem){
  for(j in 1:ndim){
    CBCL_w.est[i,j] = mean(CBCL_w.proc[,i,j])
  }
}
for(k in 1:nsample){
  for(j in 1:ndim){
    CBCL_z.est[k,j] = mean(CBCL_z.proc[,k,j])
  }
}
toc()

CBCL_beta.est=colMeans(CBCL_beta)
CBCL_theta.est=colMeans(CBCL_theta)


CBCL_w.est_pr=procrustes(CBCL_w.est,YSR_w.est)$X.new #procruste, just for visulize
CBCL_z.est_pr=procrustes(CBCL_z.est,YSR_z.est)$X.new #procruste, just for visulize

# Caculate basic matrix---------------------------------------------------------------------------------------

tic()
integ_mat=matrix(NA,nrow=nitem,ncol=nitem)
for(i in 1:nitem){
  for(j in 1:nitem){
    if(i>=j){next}
    else{
      CBCL_dist=(rowSums((CBCL_w.proc[,i,]-CBCL_w.proc[,j,])**2))**0.5
      YSR_dist=(rowSums((YSR_w.proc[,i,]-YSR_w.proc[,j,])**2))**0.5
      integ_mat[i,j]=num_integral(YSR_dist,CBCL_dist)
    }
  }
}
toc()
integ_mat=as.matrix(integ_mat)
colnames(integ_mat)=item_name
rownames(integ_mat)=item_name
sym_integ_mat=t(replace(integ_mat,is.na(integ_mat),0))+replace(integ_mat,is.na(integ_mat),0)
# write.table(sym_integ_mat,'../sym_integ_mat.txt',row.names = F,col.names = F)

tic()
integ_mat=matrix(NA,nrow=nitem,ncol=nitem)
for(i in 1:nitem){
    for(j in 1:nitem){
        if(i==j){next}
        CBCL_dist=(rowSums((CBCL_w.proc[,i,]-CBCL_w.proc[,j,])**2))**0.5
        YSR_dist=(rowSums((YSR_w.proc[,i,]-YSR_w.proc[,j,])**2))**0.5
        integ_mat[i,j]=num_integral_v2(YSR_dist,CBCL_dist)
    }
}
toc()
asym_integ_mat=integ_mat
colnames(asym_integ_mat)=(group_name)
rownames(asym_integ_mat)=(group_name)
# write.table(asym_integ_mat,'../asym_integ_mat.txt',row.names = F,col.names = F)
# asym_integ_mat=as.matrix(read.table('../asym_integ_mat.txt'))


tic()
CBCL_dist_mat=matrix(NA,nrow=nitem,ncol=nitem)
YSR_dist_mat=matrix(NA,nrow=nitem,ncol=nitem)
i=1;j=2
for(i in 1:nitem){
    for(j in 1:nitem){
        if(i==j){next}
        CBCL_dist=(rowSums((CBCL_w.proc[,i,]-CBCL_w.proc[,j,])**2))**0.5
        tmp1=kde(CBCL_dist)
        CBCL_dist_mat[i,j]=tmp1$eval.points[which.max(tmp1$estimate)]

        YSR_dist=(rowSums((YSR_w.proc[,i,]-YSR_w.proc[,j,])**2))**0.5
        tmp2=kde(YSR_dist)
        YSR_dist_mat[i,j]=tmp2$eval.points[which.max(tmp2$estimate)]
    }
}
toc()
# write.table(CBCL_dist_mat,'../CBCL_dist_mat.txt',row.names = F,col.names = F)
# write.table(YSR_dist_mat,'../YSR_dist_mat.txt',row.names = F,col.names = F)
# YSR_dist_mat=as.matrix(read.table('../YSR_dist_mat.txt'))
# CBCL_dist_mat=as.matrix(read.table('../CBCL_dist_mat.txt'))
colnames(CBCL_dist_mat)=(group_name)
rownames(CBCL_dist_mat)=(group_name)
colnames(YSR_dist_mat)=(group_name)
rownames(YSR_dist_mat)=(group_name)

# Calculate KL_v2
tic()
asy_dist_mat=matrix(NA,nrow=nitem,ncol=nitem)
for(i in 1:nitem){
  for(j in 1:nitem){
    if(i==j){next}
    CBCL_dist=(rowSums((CBCL_w.proc[,i,]-CBCL_w.proc[,j,])**2))**0.5
    YSR_dist=(rowSums((YSR_w.proc[,i,]-YSR_w.proc[,j,])**2))**0.5

    if(i>=j){
      asy_dist_mat[i,j]=KL_divergence_v2(CBCL_dist,YSR_dist)
    }
    else{
      asy_dist_mat[i,j]=KL_divergence_v2(YSR_dist,CBCL_dist)
    }
  }
}
toc()
asy_dist_mat=as.matrix(asy_dist_mat)
colnames(asy_dist_mat)=item_name
rownames(asy_dist_mat)=item_name
# write.table(asy_dist_mat,'../scaled_KL_asymm_dist_v2.txt',row.names = F,col.names = F)

tic()
asy_dist_mat=matrix(NA,nrow=nitem,ncol=nitem)
for(i in 1:nitem){
  for(j in 1:nitem){
    if(i==j){next}
    CBCL_dist=(rowSums((CBCL_w.proc[,i,]-CBCL_w.proc[,j,])**2))**0.5
    YSR_dist=(rowSums((YSR_w.proc[,i,]-YSR_w.proc[,j,])**2))**0.5

    if(i>=j){
      asy_dist_mat[i,j]=mean(CBCL_dist)-mean(YSR_dist) 
    }
    else{
      asy_dist_mat[i,j]=mean(YSR_dist)-mean(CBCL_dist) 
    }
  }
}
toc()
asy_dist_mat=as.matrix(asy_dist_mat)
colnames(asy_dist_mat)=item_name
rownames(asy_dist_mat)=item_name
# write.table(asy_dist_mat,'../map_comparison_for_wi.txt',row.names = F,col.names = F)

nitem=dim(CBCL_w.proc)[2]
dim(CBCL_w.proc)
# i=1;j=2
tic()
dist_mat=matrix(NA,nrow=nitem,ncol=nitem)
for(i in 1:nitem){
  for(j in 1:nitem){
    if(i>=j){next}
    else{
      CBCL_dist=(rowSums((CBCL_w.proc[,i,]-CBCL_w.proc[,j,])**2))**0.5
      YSR_dist=(rowSums((YSR_w.proc[,i,]-YSR_w.proc[,j,])**2))**0.5
      dist_mat[i,j]=KL_divergence(YSR_dist,CBCL_dist) 
    }
  }
}
toc()
dist_mat=as.matrix(dist_mat)
colnames(dist_mat)=item_name
rownames(dist_mat)=item_name

hist(as.numeric(dist_mat),nclass=50)
sym_dist=t(replace(dist_mat,is.na(dist_mat),0))+replace(dist_mat,is.na(dist_mat),0)
# write.table(sym_dist,'../scaled_KL_symm_dist.txt',row.names = F,col.names = F)
# write.table(sym_dist,'../KL_symm_dist.txt',row.names = F,col.names = F)


# dup_dist_itemwise.png-------------------------------------
# dup_dist_groupwise.png-------------------------------------

# Item pair with smallest duplicated pdf_mass ------------------------------------------------------------
asym_integ_mat=as.matrix(read.table('../asym_integ_mat.txt'))
row_idx=which((asym_integ_mat<sort(asym_integ_mat,decreasing = F)[20]))%/%nitem+1
col_idx=which((asym_integ_mat<sort(asym_integ_mat,decreasing = F)[20]))%%nitem
chk_idx2=cbind(row_idx,col_idx)
cbind(as.numeric(item_name[row_idx]),as.numeric(item_name[col_idx]))


# Item pair with largest KL-divergence ------------------------------------------------------------
sym_dist=as.matrix(read.table('../scaled_KL_symm_dist.txt'))
row_idx=which((sym_dist>sort(sym_dist,decreasing = T)[10]))%/%nitem+1
col_idx=which((sym_dist>sort(sym_dist,decreasing = T)[10]))%%nitem
chk_idx=cbind(row_idx,col_idx)
cbind(as.numeric(item_name[row_idx]),as.numeric(item_name[col_idx]))

# Compare KL and num_int ------------------------------------------------------------
row_idx=which((asym_integ_mat>sort(asym_integ_mat,decreasing = T)[20]))%/%nitem+1
col_idx=which((asym_integ_mat>sort(asym_integ_mat,decreasing = T)[20]))%%nitem
chk_idx3=cbind(row_idx,col_idx)

tmp1=as.numeric(asym_integ_mat)
tmp2=as.numeric(sym_dist)
cor(tmp1[!is.na(tmp1)],tmp2[!is.na(tmp1)])
par(mfrow=c(2,1))
#chk_ii=3
for(chk_ii in 1:dim(chk_idx3)[1]){
  i=chk_idx3[chk_ii,1]
  j=chk_idx3[chk_ii,2]
  CBCL_dist=(rowSums((CBCL_w.proc[,i,]-CBCL_w.proc[,j,])**2))**0.5
  YSR_dist=(rowSums((YSR_w.proc[,i,]-YSR_w.proc[,j,])**2))**0.5
  tmp_KL=KL_divergence_v2(YSR_dist,CBCL_dist)
  tmp_int=num_integral_v2(YSR_dist,CBCL_dist)
  hist(CBCL_dist,main=sprintf('KL is %s',tmp_KL),nclass=100)
  hist(YSR_dist,main=sprintf('num_int  is %s',tmp_int),nclass=100)
}
par(mfrow=c(1,1))



# Calculate t.test p.value for detecting different distribution-----------------------------------------
tic()
ttest_mat=matrix(NA,nrow=nitem,ncol=nitem)
for(i in 1:nitem){
    for(j in 1:nitem){
        if(i==j){next}
        CBCL_dist=(rowSums((CBCL_w.proc[,i,]-CBCL_w.proc[,j,])**2))**0.5
        YSR_dist=(rowSums((YSR_w.proc[,i,]-YSR_w.proc[,j,])**2))**0.5
        ttest_mat[i,j]=t.test(CBCL_dist, YSR_dist, alternative = "two.sided", var.equal = FALSE)$p.value
    }
}
toc()
asym_ttest_mat=ttest_mat
hist(asym_ttest_mat,nclass=100)
pheatmap(asym_ttest_mat,cluster_rows = F,cluster_cols = F)
# write.table(asym_ttest_mat,'../asym_ttest_mat.txt',row.names = F,col.names = F)

# CBCL
tic()
CBCL_DA_jaccard_sum=matrix(0,ncol=118,nrow=118)
niter=dim(CBCL_impute)[1]
CBCL_DA_org=CBCL_data
impute_idx = CBCL_DA_org==999
for(t in 1:niter){
CBCL_DA=CBCL_DA_org
CBCL_DA[impute_idx]=CBCL_impute[t,] # i,k
  CBCL_DA_jaccard_sum=CBCL_DA_jaccard_sum+Jaccard_sim(CBCL_DA)
}
CBCL_DA_jaccard_mean=CBCL_DA_jaccard_sum/niter
toc()


# YSR
tic()
YSR_DA_jaccard_sum=matrix(0,ncol=118,nrow=118)
niter=dim(YSR_impute)[1]
YSR_DA_org=as.matrix(read.table('../YSR/data/item.txt'))
impute_idx = YSR_DA_org==999
for(t in 1:niter){
  YSR_DA=YSR_DA_org
  YSR_DA[impute_idx]=YSR_impute[t,] 
  YSR_DA_jaccard_sum=YSR_DA_jaccard_sum+Jaccard_sim(YSR_DA)
}
YSR_DA_jaccard_mean=YSR_DA_jaccard_sum/niter
toc()

# write.table(CBCL_DA_jaccard_mean,'../CBCL_DA_jaccard_mean.txt',row.names = F,col.names = F)
# write.table(YSR_DA_jaccard_mean,'../YSR_DA_jaccard_mean.txt',row.names = F,col.names = F)

# text_cex=1.2
up_shift=0.25
axis_size=1.5

par(mfrow=c(1,4))
for(idx in 1:10){
  par(mar = c(5,5,3,3))
  hist(CBCL_theta[,4],nclass=50,
       xlab=expression(hat(theta)[CBCL4]),main='',
       cex.lab=1.5,cex.axis=axis_size)
}

par(mfrow=c(1,1))
plot(CBCL_w.proc[1,1,],type='n',ylim=c(-4,4),xlim=c(-4,4))
for(i_idx in 1:10){
  for(idx in 1:500){
    j_idx=idx*10
    points(CBCL_w.proc[j_idx,i_idx,1],CBCL_w.proc[j_idx,i_idx,2],col=alpha(i_idx,0.3),pch=19)
  }
}



# tmporal function
gr_cent=function(inp_group='CBCL_w.proc_w_YSR 2'){
  return(colMedians(as.matrix(dfff1[dfff1$idx==inp_group,1:2])))
}

# visualize ggplot--------------------------------------------------------------------------------



#Visualize latent coordinate of Most simliar W,Z & Most different W,Z for CBCL,YSR-----------------------------------------------------

w.est_diff=rowSums(abs(YSR_w.est-CBCL_w.est_pr))
targ_idx_w=order(w.est_diff)[c(1,2,length(w.est_diff)-1,length(w.est_diff))]

z.est_diff=rowSums(abs(YSR_z.est-CBCL_z.est_pr))
targ_idx_z=order(z.est_diff)[c(1,2,length(z.est_diff)-1,length(z.est_diff))]

z.est_diff[order(z.est_diff)]
w.est_diff[order(w.est_diff)]
plot(YSR_z.est[],type='n',ylim=c(-4,4),xlim=c(-4,4))
tmp_idx=4
points(YSR_z.est[targ_idx_z[tmp_idx],1],YSR_z.est[targ_idx_z[tmp_idx],2],col=2)
points(CBCL_z.est_pr[targ_idx_z[tmp_idx],1],CBCL_z.est_pr[targ_idx_z[tmp_idx],2],col=1)


#CBCL_coord_gg_diff.png-----------------------------------------------------  
item_name[targ_idx_w]
inp_size=1.5
coord_plot=function(inp_w_mat,inp_z_mat,targ_idx_w,targ_idx_z){
  targ_idx=(1:500)*10#1:5000#(1:100)*50
  
  
  i_idx=targ_idx_w[1]
  df1=as.data.frame(inp_w_mat[targ_idx,i_idx,]);df1[,3]=paste('w',as.character(item_name[i_idx]),sep = '');colnames(df1)=c('w1','w2','idx')
  for(i_idx in targ_idx_w[2:4]){
    df0=as.data.frame(inp_w_mat[targ_idx,i_idx,]);df0[,3]=paste('w',as.character(item_name[i_idx]),sep = '');colnames(df0)=c('w1','w2','idx')
    df1=rbind(df1,df0)
  }
  
  i_idx=targ_idx_z[1]
  dff1=as.data.frame(inp_z_mat[targ_idx,i_idx,]);dff1[,3]=paste('z',as.character(i_idx),sep = '');colnames(dff1)=c('w1','w2','idx')
  for(i_idx in targ_idx_z[2:4]){
    dff0=as.data.frame(inp_z_mat[targ_idx,i_idx,]);dff0[,3]=paste('z',as.character(i_idx),sep = '');colnames(dff0)=c('w1','w2','idx')
    dff1=rbind(dff1,dff0)
  }
  
  
  df1[,4]='item'
  dff1[,4]='respondents'
  dfff1=rbind(df1,dff1)
  colnames(dfff1)=c('w1','w2','idx','item_resp')
  # tmporal function
  gr_cent=function(inp_group='inp_w_mat 2'){
    return(colMedians(as.matrix(dfff1[dfff1$idx==inp_group,1:2])))
  }
  
  
  w_char_list=paste0('w',item_name[targ_idx_w])
  z_char_list=paste('z',as.character(targ_idx_z),sep='')  
  # visualize ggplot
  ### CBCL_coord_gg.png
  ggplot(dfff1[,],aes(w1,w2,col=item_resp,shape=item_resp,group=idx,label=idx))+ theme_bw()+###
    geom_point(size=1.7/inp_size,alpha=0.8)+
    geom_density_2d(size=0.7/inp_size)+
    # xlim(-4,4)+ylim(-4,4)+
    xlim(-3,3)+ylim(-4,4)+
    coord_cartesian(expand = FALSE)+
    annotate(geom='text',x=gr_cent(w_char_list[1])[1], y=gr_cent(w_char_list[1])[2], label=w_char_list[1], color='black',size=7/inp_size)+
    annotate(geom='text',x=gr_cent(w_char_list[2])[1], y=gr_cent(w_char_list[2])[2], label=w_char_list[2], color='black',size=7/inp_size)+
    annotate(geom='text',x=gr_cent(w_char_list[3])[1], y=gr_cent(w_char_list[3])[2], label=w_char_list[3], color='black',size=7/inp_size)+
    annotate(geom='text',x=gr_cent(w_char_list[4])[1], y=gr_cent(w_char_list[4])[2], label=w_char_list[4], color='black',size=7/inp_size)+
    annotate(geom='text',x=gr_cent(z_char_list[1])[1], y=gr_cent(z_char_list[1])[2], label=z_char_list[1], color='black',size=7/inp_size)+
    annotate(geom='text',x=gr_cent(z_char_list[2])[1], y=gr_cent(z_char_list[2])[2], label=z_char_list[2], color='black',size=7/inp_size)+
    annotate(geom='text',x=gr_cent(z_char_list[3])[1], y=gr_cent(z_char_list[3])[2], label=z_char_list[3], color='black',size=7/inp_size)+
    annotate(geom='text',x=gr_cent(z_char_list[4])[1], y=gr_cent(z_char_list[4])[2], label=z_char_list[4], color='black',size=7/inp_size)+
    theme(axis.text=element_text(size=16),axis.title=element_text(size=14,face="bold"))+
    theme(plot.margin = unit(c(1,1,1,1), "cm")) +
    theme(legend.position = "none") + 
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    scale_color_manual(values=gg_color_hue(4)[c(2,4)])
  
}

coord_plot(YSR_w.proc,YSR_z.proc,targ_idx_w,targ_idx_z)
# ggsave(file='YSR_coord_gg_diff_highres_v5.png', width=20/inp_size, height=18/inp_size, units=c("cm"),dpi=600)

{YSR_max.address = which.max(YSR_mle)
  YSR_w.star = YSR_w[YSR_max.address,,]
  YSR_z.star = YSR_z[YSR_max.address,,]
  CBCL_w.proc_w_YSR = array(0,dim=c(nmcmc,nitem,ndim))
  CBCL_z.proc_w_YSR = array(0,dim=c(nmcmc,nsample,ndim))}

# Procrustes matching for latent spaces Z (respondent) and W (item)
for(iter in 1:nmcmc){
  CBCL_z.iter = CBCL_z[iter,,]
  CBCL_w.iter = CBCL_w[iter,,]
  CBCL_w.proc_w_YSR[iter,,] = procrustes(CBCL_w.iter,YSR_w.star)$X.new
  CBCL_z.proc_w_YSR[iter,,] = procrustes(CBCL_z.iter,YSR_z.star)$X.new
}

coord_plot(CBCL_w.proc_w_YSR,CBCL_z.proc_w_YSR,targ_idx_w,targ_idx_z)
# ggsave(file='CBCL_coord_gg_diff_highres_v5.png', width=20/inp_size, height=18/inp_size, units=c("cm"),dpi=600)


#beta plot & theta plot-------------------------------------------
outlier_list=(CBCL_beta.est<(0.8))&(YSR_beta.est>2.5)
df=as.data.frame(cbind(CBCL_beta.est,YSR_beta.est))
df[,3]=outlier_list
df[,4]=item_name
colnames(df)=c('CBCL_beta','YSR_beta','outlier','Name')
lm1=lm(YSR_beta.est[!(outlier_list)]~CBCL_beta.est[!(outlier_list)])
lm2=lm(YSR_beta.est~CBCL_beta.est)
cor(YSR_beta.est[!(outlier_list)],CBCL_beta.est[!(outlier_list)]) # correlation : 0.8349648
cor(YSR_beta.est,CBCL_beta.est) # correlation : 0.01974946

ggplot(df,aes(x=CBCL_beta,y=YSR_beta,label=Name))+
  ylim(-2,6.5)+
  geom_point(size=1.5,alpha=1,shape=19)+
  geom_point(size=5,col='black',alpha=1,shape=20,data=df[df[,3]==T,1:4])+
  geom_text(aes(label=ifelse((CBCL_beta<(0.8))&(YSR_beta>2.5),as.character(item_name),'')),hjust=0,vjust=-1,col='blue',size=5)+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=16)) + 
  labs(x = expression(hat(beta)[CBCL]),y = expression(hat(beta)[YSR]))+
  geom_abline(intercept= 0, slope=1, color='red', size = 1)+
  geom_abline(intercept= lm1$coefficients[1], slope=lm1$coefficients[2], color='black', size = 1)

# Automatic outlier detectoin using Robust linear regression
rlm2 <- rlm(YSR_beta.est~CBCL_beta.est)
tmp=rlm2$residuals
tmp_names=names(head(tmp[order(tmp,decreasing = T)],14))
outlier_list22=as.numeric(sapply(tmp_names,FUN = function(x)substring(text = x,first = 2,last = 10)))
lm3=lm(YSR_beta.est[-(outlier_list22)]~CBCL_beta.est[-(outlier_list22)])

df22=df
df22[,3]=F
df22[outlier_list22,3]=T
ggplot(df22,aes(x=CBCL_beta,y=YSR_beta,label=Name))+
  ylim(-2.4,7.4)+xlim(-4,3.9)+coord_cartesian(expand = FALSE)+theme_bw()+
  geom_point(size=2/inp_size,alpha=1,shape=19)+
  geom_point(size=5/inp_size,col='black',alpha=1,shape=20,data=df22[df22[,3]==T,1:4])+
  geom_text(aes(label=ifelse((CBCL_beta<(0.8))&(YSR_beta>2.5),as.character(item_name),'')),hjust=0,vjust=-1,col='blue',size=8/inp_size)+
  
  labs(x = expression(hat(beta)[CBCL]),y = expression(hat(beta)[YSR]))+
  geom_abline(intercept= 0, slope=1, color='black', size = 1/inp_size,linetype = "dashed")+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=16)) + 
  geom_abline(intercept= lm3$coefficients[1], slope=lm3$coefficients[2], color='red', size = 1/inp_size)
# ggsave(file='beta_plot_gg_highres_v4.png', width=20/inp_size, height=18/inp_size, units=c("cm"),dpi=600)

plot(diff(tmp[order(tmp,decreasing = T)[1:20]]))
which.min(diff(tmp[order(tmp,decreasing = T)[1:20]])) # The number of outlier is selected by scree plot!

tmp=rbind(df22[df22[,3],1:3],df22[!df22[,3],1:3])
write.csv(tmp,file="Figure4_data.csv",row.names=F)


# Theata plot
df=as.data.frame(cbind(CBCL_theta.est,YSR_theta.est))
colnames(df)=c('CBCL_theta','YSR_theta')
lm3=lm(YSR_theta.est~CBCL_theta.est)
ggplot(df,aes(x=CBCL_theta,y=YSR_theta))+
  geom_point(size=1.5,alpha=1,shape=19)+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=16)) + 
  labs(x = expression(hat(theta)[CBCL]),y = expression(hat(theta)[YSR]))+
  geom_abline(intercept= 0, slope=1, color='red', size = 1)+theme_bw()+
  geom_abline(intercept= lm3$coefficients[1], slope=lm3$coefficients[2], color='black', size = 1)
cor(YSR_theta.est,CBCL_theta.est) # 0.3503182
ggsave(file='theta_plot_gg_highres_v2.png', width=20, height=18, units=c("cm"),dpi=600)

#rlm
rlm2 <- rlm(YSR_theta.est~CBCL_theta.est)
tmp=rlm2$residuals
plot(tmp[order(tmp,decreasing = T)[1:12]])
plot(diff(tmp[order(tmp,decreasing = T)[1:100]]))
which.min(diff(tmp[order(tmp,decreasing = T)[1:11]]))


tmp_names=names(head(tmp[order(tmp,decreasing = T)],12))
outlier_list22=as.numeric(sapply(tmp_names,FUN = function(x)substring(text = x,first = 2,last = 10)))
lm3=lm(YSR_theta.est[-(outlier_list22)]~CBCL_theta.est[-(outlier_list22)])


df22=df
df22[,3]=F
df22[outlier_list22,3]=T
df22[,4]=seq(1,662)
colnames(df22)=c('CBCL_theta','YSR_theta','outlier','Name')
ggplot(df22,aes(x=CBCL_theta,y=YSR_theta,label=Name))+
  #ylim(-2,6.5)+
  geom_point(size=1.5,alpha=1,shape=19)+
  geom_point(size=5,col='black',alpha=1,shape=20,data=df22[df22[,3]==T,1:4])+
  geom_text(aes(label=ifelse(outlier,as.character(seq(1,662)),'')),hjust=0,vjust=-1,col='blue',size=5)+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=16)) + 
  labs(x = expression(hat(theta)[CBCL]),y = expression(hat(theta)[YSR]))+
  geom_abline(intercept= 0, slope=1, color='red', size = 1)+
  geom_abline(intercept= lm3$coefficients[1], slope=lm3$coefficients[2], color='black', size = 1)


# Coordinate that differs in beta coeff------------------------------------
df1=as.data.frame(CBCL_w.est_pr)
df1[,3]='2'
colnames(df1)=c('w1','w2','source')
df2=as.data.frame(YSR_w.est)
df2[,3]='1'
colnames(df2)=c('w1','w2','source')
df=rbind(df1,df2)
#item_name[j_list]

target_idx=c('06','15','49','59','60','73','80','88','92','98','106','107','108','109')
c(item_name %in%target_idx,item_name %in%target_idx)
ggplot(df,aes(w1,w2,col=source))+
  stat_ellipse(data=df2[item_name %in%target_idx,],type='t',linetype=2,size=1.5)+theme_bw()+
  geom_point(size=4/inp_size,alpha=1,shape=20)+
  
  # geom_point(size=7,alpha=0.9,shape=21,stroke=2.5,col='black',data=df1[250,])+
  # geom_point(aes(fill=source), colour="black",pch=21, size=5,data=df1[item_name %in%target_idx,])+
  geom_point(aes(fill=source), colour="black",pch=21, size=6/inp_size,alpha=1,data=df1[item_name %in%target_idx,])+
  geom_point(aes(fill=source), colour="black",pch=21, size=6/inp_size,alpha=1,data=df2[item_name %in%target_idx,])+
  xlim(-3,3)+ylim(-3,3)+coord_cartesian(expand = FALSE) +
  
  geom_text(aes(label=ifelse(c(item_name %in%target_idx,item_name %in%target_idx),c(as.character(item_name),as.character(item_name)),'')),hjust=0.5,vjust=-1,size=6/inp_size,show.legend = F)+
  
  
  theme(axis.text=element_text(size=16),axis.title=element_text(size=14,face="bold"))+
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),legend.title=element_blank())+ 
  
  theme(legend.position = c(0.9,0.9))+
  scale_fill_manual(values=gg_color_hue(2),labels=c("YSR", "CBCL"))+
  scale_color_manual(values=gg_color_hue(2),labels=c("YSR", "CBCL"))

# ggsave(file='beta_differ_inboth_gg_highres_v7.png', width=20/inp_size, height=18/inp_size, units=c("cm"),dpi=600)


target_idx=high_KL
ggplot(df,aes(w1,w2,col=source))+
  stat_ellipse(data=df2[item_name %in%target_idx,],type='t',linetype=2,size=1.5)+theme_bw()+
  geom_point(size=4/inp_size,alpha=1,shape=20)+
  geom_point(aes(fill=source), colour="black",pch=21, size=6/inp_size,alpha=1,data=df1[item_name %in%target_idx,])+
  geom_point(aes(fill=source), colour="black",pch=21, size=6/inp_size,alpha=1,data=df2[item_name %in%target_idx,])+  
  # geom_point(size=7/inp_size,alpha=1,shape=20,data=df1[item_name %in%target_idx,])+
  # geom_point(size=7/inp_size,alpha=1,shape=20,data=df2[item_name %in%target_idx,])+
  xlim(-3,3)+ylim(-3,3)+coord_cartesian(expand = FALSE) +
  
  geom_text(aes(label=ifelse(c(item_name %in%target_idx,item_name %in%target_idx),c(as.character(item_name),as.character(item_name)),'')),hjust=0.5,vjust=-1,size=6/inp_size,show.legend = F)+
  
  theme(axis.text=element_text(size=16),axis.title=element_text(size=14,face="bold"))+
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),legend.title=element_blank())+ 
  
  theme(legend.position = c(0.9,0.9))+
  scale_fill_manual(values=gg_color_hue(2),labels=c("YSR", "CBCL"))+
  scale_color_manual(values=gg_color_hue(2),labels=c("YSR", "CBCL"))


#item subgroup centroid plot--------------------------

#non centered version 
CBCL_centroidmat=matrix(NA,ncol=2,nrow=length(data_list))
YSR_centroidmat=matrix(NA,ncol=2,nrow=length(data_list))
for (idx in 1:length(data_list)){
  #set target_idx  
  target_idx=data_list[[idx]]
  
  #calculate centroid or group distance, for each target_idx
  CBCL_centroid=c(mean(CBCL_w.est_pr[item_name %in%target_idx,1]),mean(CBCL_w.est_pr[item_name %in%target_idx,2]))
  YSR_centroid=c(mean(YSR_w.est[item_name %in%target_idx,1]),mean(YSR_w.est[item_name %in%target_idx,2]))
  CBCL_centroidmat[idx,]=CBCL_centroid
  YSR_centroidmat[idx,]=YSR_centroid
}



#Respondent specific plot--------------------------------
inp_type='YSR'
{if(inp_type=='CBCL'){
  inp_w.est=CBCL_w.est_pr
  inp_z.est=CBCL_z.est_pr
  inp_centroid_mat=CBCL_centroidmat
}
  else{
    inp_w.est=YSR_w.est
    inp_z.est=YSR_z.est
    inp_centroid_mat=YSR_centroidmat
  }
}


centroidmat=matrix(NA,ncol=2,nrow=length(data_list))
for (idx in 1:length(data_list)){
  #set target_idx  
  target_idx=data_list[[idx]]
  
  #calculate centroid or group distance, for each target_idx
  tmp_centroid=c(mean(inp_w.est[item_name %in%target_idx,1]),mean(inp_w.est[item_name %in%target_idx,2]))
  centroidmat[idx,]=tmp_centroid
}

# KS statistics, not KL-divergence for detecting distribution discrepancy--------------------------------------------------------------
tic()
KS_mat=matrix(NA,nrow=nitem,ncol=nitem)
for(i in 1:nitem){
  for(j in 1:nitem){
    CBCL_dist=(rowSums((CBCL_w.proc[,i,]-CBCL_w.proc[,j,])**2))**0.5
    YSR_dist=(rowSums((YSR_w.proc[,i,]-YSR_w.proc[,j,])**2))**0.5
    KS_mat[i,j]=ks.test(CBCL_dist,YSR_dist)$p.value
  }
}
toc()
KS_mat=as.matrix(KS_mat)
colnames(KS_mat)=item_name
rownames(KS_mat)=item_name
#write.table(KS_mat,'../KS_mat.txt',row.names = F,col.names = F)


# sym_dist=as.matrix(read.table('../KL_symm_dist.txt'))

### indi_plot2_CBCL1_v2 , indi_plot2_CBCL3_v2
# visulization
par(mfrow=c(1,1),mar=rep(3,4))
#for(inp_idx in c(indi_in_out,indi_in_hub)){
for(inp_idx in c(1,3,16,102,250)){
  cat(sprintf('%s, %s',inp_type,inp_idx))
  plot(inp_w.est,type='n',xlim=c(-2.5,2.5),ylim=c(-3,3),xlab='',ylab='',cex.axis=axis_size)
  #points(inp_z.est,col=2)
  #points(inp_w.est)
  
  points(centroidmat,cex=2.5,pch=8,lwd=2,col=ifelse(inp_type=='YSR',2,1))
  text(centroidmat[,1],centroidmat[,2]+up_shift,name_list,lwd=2,col=ifelse(inp_type=='YSR',2,1),cex=text_cex)
  #inp_idx=3
  points(inp_z.est[inp_idx,1],inp_z.est[inp_idx,2],col='blue',cex=3,lwd=3)
  legend(x=1.25,y=3,c(inp_type,'Respondent'),col=c(ifelse(inp_type=='YSR',2,1),'blue'),pch=c(8,21),cex=1.4)
}


### indi_plot2_CBCLall_v2.png
plot(inp_w.est,type='n',xlim=c(-2.5,2.5),ylim=c(-3,3),xlab='',ylab='')
points(inp_z.est,col=alpha('blue',0.5),cex=1.2)
points(centroidmat,cex=2.5,pch=8,lwd=2,col=ifelse(inp_type=='YSR',2,1))
legend(x=1.25,y=3,c(inp_type,'Respondent'),col=c(ifelse(inp_type=='YSR',2,1),'blue'),pch=c(8,21),cex=1.4)


#EDA-------------------------------------------


stopifnot(dim(sym_dist)[1]==length(item_name))
colnames(sym_dist)=item_name
row.names(sym_dist)=item_name
#heatmap.2(sym_dist,scale='none',col=bluered(100),trace='none',symm = T,Rowv=F,Colv=F)

high_KL=names(colSums(sym_dist)[order(colSums(sym_dist,na.rm = T),decreasing = T)][1:10])
raw_data[raw_data$Variable%in%high_KL,c('Question','Variable')]

# Item pair with largest KL-divergence ------------------------------------------------------------
row_idx=which((sym_dist>sort(sym_dist,decreasing = T)[10]))%/%nitem+1
col_idx=which((sym_dist>sort(sym_dist,decreasing = T)[10]))%%nitem
cbind(as.numeric(item_name[row_idx]),as.numeric(item_name[col_idx]))


# Item pair with largest Jaccard-sim ------------------------------------------------------------
diff_of_jaccard=(CBCL_DA_jaccard_mean_dist-YSR_DA_jaccard_mean_dist)**2
to_see=(60+1)*2
row_idx=which((diff_of_jaccard>sort(diff_of_jaccard,decreasing = T)[to_see]))%/%nitem+1
col_idx=which((diff_of_jaccard>sort(diff_of_jaccard,decreasing = T)[to_see]))%%nitem
cbind(as.numeric(item_name[row_idx]),as.numeric(item_name[col_idx]))

unique(as.numeric(cbind(as.numeric(item_name[row_idx]),as.numeric(item_name[col_idx]))))
sort(as.numeric(item_name[outlier_list22]))
#Scatter plot of estimated values(theta, beta)------------------------


### beta_plot_v2.png
par(mar = c(5,5,3,3))
plot(CBCL_beta.est,YSR_beta.est,xlab=expression(hat(beta)[CBCL]),
     ylim=c(-2,6.5),
     ylab=expression(hat(beta)[YSR]),xaxt='n',yaxt='n',cex.lab=1.5,cex=1.2,lwd=1.2)
axis(2,cex.axis=axis_size);axis(1,cex.axis=axis_size)
outlier_list=(CBCL_beta.est<(0.8))&(YSR_beta.est>2.5)
check_item=item_name[outlier_list] 
text(CBCL_beta.est[outlier_list],YSR_beta.est[outlier_list]+up_shift,labels=check_item,cex=text_cex,font=2,col=4)
lm1=lm(YSR_beta.est[!(outlier_list)]~CBCL_beta.est[!(outlier_list)])
abline(a=0, b=1, col="red",lwd=2)
abline(a=lm1$coefficients[1],b=lm1$coefficients[2],lwd=2)

### theta_plot_v2.png
plot(CBCL_theta.est,YSR_theta.est,xlab=expression(hat(theta)[CBCL])
     ,ylab=expression(hat(theta)[YSR])
     ,xaxt='n',yaxt='n',cex.lab=1.5,cex=1.2)
axis(2,cex.axis=axis_size);axis(1,cex.axis=axis_size)
lm3=lm(YSR_theta.est~CBCL_theta.est)
abline(a=0, b=1, col="red",lwd=2)
abline(a=lm3$coefficients[1],b=lm3$coefficients[2],lwd=2)


#visualize unmatching items-------------------------------------------
### KL_differ_inboth_v2.png
#KL
target_idx=high_KL
#plot(CBCL_w.est_pr[,1],CBCL_w.est_pr[,2],col=1,main='different KL-items',type='n',xlim=c(-2.5,2.5),ylim=c(-3,3),xlab = '',ylab='')
plot(CBCL_w.est_pr[,1],CBCL_w.est_pr[,2],col=1,type='n',xlim=c(-2.5,2.5),ylim=c(-3,3),xlab = '',ylab='',yaxt='n',xaxt='n')
axis(2,cex.axis=axis_size);axis(1,cex.axis=axis_size)
points(CBCL_w.est_pr[,1],CBCL_w.est_pr[,2],col=alpha(1,0.4));points(YSR_w.est[,1],YSR_w.est[,2],col=alpha(2,0.4))
points(CBCL_w.est_pr[item_name %in%target_idx,1],CBCL_w.est_pr[item_name %in%target_idx,2],col=1,cex=2,pch=19)
points(YSR_w.est[item_name %in%target_idx,1],YSR_w.est[item_name %in%target_idx,2],col=2,cex=2,pch=19)
text(YSR_w.est[item_name %in%target_idx,1],YSR_w.est[item_name %in%target_idx,2]+up_shift,labels=item_name[item_name %in%target_idx],cex=text_cex,font=2,col=2)
text(CBCL_w.est_pr[item_name %in%target_idx,1],CBCL_w.est_pr[item_name %in%target_idx,2]+up_shift,labels=item_name[item_name %in%target_idx],cex=text_cex,font=2,col=1)
legend(x=1.7,y=3,c('YSR','CBCL'),col=c(2,1),pch=19,cex=1.4)

### beta_differ_inboth_v2.png 
#beta
#target_idx=c('CBCL_06','CBCL_15','CBCL_49','CBCL_59','CBCL_60','CBCL_73','CBCL_80','CBCL_88','CBCL_92','CBCL_98','CBCL_106','CBCL_107','CBCL_108','CBCL_109')
target_idx=c('06','15','49','59','60','73','80','88','92','98','106','107','108','109')
plot(CBCL_w.est_pr[,1],CBCL_w.est_pr[,2],col=1,type='n',xlim=c(-2.5,2.5),ylim=c(-3,3),xlab = '',ylab='',yaxt='n',xaxt='n')
axis(2,cex.axis=axis_size);axis(1,cex.axis=axis_size)
points(CBCL_w.est_pr[,1],CBCL_w.est_pr[,2],col=alpha(1,0.4));points(YSR_w.est[,1],YSR_w.est[,2],col=alpha(2,0.4))
points(CBCL_w.est_pr[item_name %in%target_idx,1],CBCL_w.est_pr[item_name %in%target_idx,2],col=1,cex=2,pch=19)
points(YSR_w.est[item_name %in%target_idx,1],YSR_w.est[item_name %in%target_idx,2],col=2,cex=2,pch=19)
up_shift=0.25
text(YSR_w.est[item_name %in%target_idx,1],YSR_w.est[item_name %in%target_idx,2]+up_shift,labels=item_name[item_name %in%target_idx],cex=1.2,font=2,col=2)
text(CBCL_w.est_pr[item_name %in%target_idx,1],CBCL_w.est_pr[item_name %in%target_idx,2]+up_shift,labels=item_name[item_name %in%target_idx],cex=1.2,font=2,col=1)
legend(x=1.7,y=3,c('YSR','CBCL'),col=c(2,1),pch=19,cex=1.4)

# Calculating item subgroup Correlation matrix-----------------------------------------
group_name=as.vector(item_name)
head(group_name)
for(idx in 1:length(data_list)){
  #idx=1
  this_group=data_list[[idx]]
  group_name[group_name%in%this_group]=name_list[idx]
}
## Corr with all samples start
getwd()
YSR_DA_org=as.matrix(read.table('../YSR/data/item.txt'))
impute_idx = YSR_DA_org==999
YSR_DA=YSR_DA_org
YSR_DA[impute_idx]=YSR_impute[1,]
t_YSR_DA=as.data.frame(t(YSR_DA))
t_YSR_DA[,(dim(t_YSR_DA)[2]+1)]=group_name
t_YSR_DA[,dim(t_YSR_DA)[2]]

tmp=aggregate(t_YSR_DA[,1], by=list(t_YSR_DA[,dim(t_YSR_DA)[2]]), FUN=sum)
tmp_df=t(as.matrix(tmp[24:31,2]))
colnames(tmp_df)=tmp[24:31,1]
merge_df=tmp_df
for(p_idx in 2:(dim(t_YSR_DA)[2]-1)){
  tmp=aggregate(t_YSR_DA[,p_idx], by=list(t_YSR_DA[,dim(t_YSR_DA)[2]]), FUN=sum)
  tmp_df=t(as.matrix(tmp[24:31,2]))
  colnames(tmp_df)=tmp[24:31,1]
  merge_df=rbind(merge_df,tmp_df)
}
merge_df=merge_df[,c('RBB','AB','AP','SP','TP','WD','SC','AD')]
corr_df=cor(merge_df)
niter=dim(CBCL_beta)[1]
tic()
for(t in 2:niter){
  YSR_DA=YSR_DA_org
  YSR_DA[impute_idx]=YSR_impute[t,] # i,k
  t_YSR_DA=as.data.frame(t(YSR_DA))
  t_YSR_DA[,(dim(t_YSR_DA)[2]+1)]=group_name
  t_YSR_DA[,dim(t_YSR_DA)[2]]

  tmp=aggregate(t_YSR_DA[,1], by=list(t_YSR_DA[,dim(t_YSR_DA)[2]]), FUN=sum)
  tmp_df=t(as.matrix(tmp[24:31,2]))
  colnames(tmp_df)=tmp[24:31,1]
  merge_df=tmp_df
  for(p_idx in 2:(dim(t_YSR_DA)[2]-1)){
    tmp=aggregate(t_YSR_DA[,p_idx], by=list(t_YSR_DA[,dim(t_YSR_DA)[2]]), FUN=sum)
    tmp_df=t(as.matrix(tmp[24:31,2]))
    colnames(tmp_df)=tmp[24:31,1]
    merge_df=rbind(merge_df,tmp_df)
  }
  merge_df=merge_df[,c('RBB','AB','AP','SP','TP','WD','SC','AD')]
  corr_df=corr_df+cor(merge_df)
}
toc()
corr_df=corr_df/(t-1)
dim(t_YSR_DA)
(merge_df)
inp_breaks=seq(-1,1,length.out = 100)

### CBCL_cosine_sim_v2.png
pheatmap(round(CBCL_centroid_cosine,3),main='', cluster_rows = F,cluster_cols = F,display_numbers=T,fontsize = 20,breaks = inp_breaks)
### YSR_cosine_sim_v2.png
pheatmap(round(YSR_centroid_cosine,3),main='', cluster_rows = F,cluster_cols = F,display_numbers=T,fontsize = 20,breaks=inp_breaks)


# calcaulte in groupwise and write dbscan result to toatl_df ---------------------
colnames(CBCL_CC_jaccard)
for (idx in 1:length(data_list)){
  target_idx=data_list[[idx]]
  group_KL=sym_dist[item_name %in%target_idx,item_name %in%target_idx]
  # to_clust=F
  # pheatmap(group_KL,cluster_rows = to_clust,cluster_cols = to_clust,main=name_list[idx])
  to_clust=T
  pheatmap(group_KL,cluster_rows = to_clust,cluster_cols = to_clust,main=name_list[idx])
  
  #plot(hclust(dist(group_KL), method="single"))
  db_result=dbscan(dist(group_KL),eps = median(dist(group_KL)),MinPts = 2)$cluster
  for(cl_idx in 1:max(db_result)){
    cl_items=item_name[item_name %in%target_idx][db_result==cl_idx]
    sbgroup_name=name_list[idx]
    cat(sbgroup_name,'\n')
    tmp_df=cbind(sbgroup_name,cl_idx,raw_data[raw_data$Variable%in%cl_items,c('Question','Variable')])
    if((cl_idx==1)&(idx==1)){
      total_df=tmp_df}
    else{
      total_df=rbind(total_df,tmp_df)}
  }
}
(total_df)
#write.table(total_df,'../../total_df.txt',row.names = F,col.names = F)
# total_df=as.data.frame(read.table('./total_df.txt'))
colnames(total_df)=c('subgroup','internal_cluster','Question','Variable')
head(total_df)

tmp=as.dist(group_KL)
high_value=tmp[order(tmp,decreasing = T)[1:3]]
which(group_KL%in%high_value)


# #centered version
centered_YSR_w.est=scale(YSR_w.est,center = T,scale = F)
# plot(YSR_w.est,main='centering YSR_w.est')
# points(centered_YSR_w.est,col=2)
# 
centered_CBCL_w.est_pr=scale(CBCL_w.est_pr,center = T,scale = F)
# plot(CBCL_w.est_pr,main='centering CBCL_w.est')
# points(centered_CBCL_w.est_pr,col=2)


CBCL_centroidmat=matrix(NA,ncol=2,nrow=length(data_list))
YSR_centroidmat=matrix(NA,ncol=2,nrow=length(data_list))
for (idx in 1:length(data_list)){
  #set target_idx  
  target_idx=data_list[[idx]]
  
  #calculate centroid or group distance, for each target_idx
  CBCL_centroid=c(mean(centered_CBCL_w.est_pr[item_name %in%target_idx,1]),mean(centered_CBCL_w.est_pr[item_name %in%target_idx,2]))
  YSR_centroid=c(mean(centered_YSR_w.est[item_name %in%target_idx,1]),mean(centered_YSR_w.est[item_name %in%target_idx,2]))
  CBCL_centroidmat[idx,]=CBCL_centroid
  YSR_centroidmat[idx,]=YSR_centroid
}
# plot(CBCL_w.est_pr,type='n',main='centroid plot, centering version')
# points(CBCL_centroidmat[,1],CBCL_centroidmat[,2],col=1)
# points(YSR_centroidmat[,1],YSR_centroidmat[,2],col=2)
# text(YSR_centroidmat[,1],YSR_centroidmat[,2]+0.2,labels=name_list,cex=0.75,font=2,col=2)
# text(CBCL_centroidmat[,1],CBCL_centroidmat[,2]+0.2,labels=name_list,cex=0.75,font=2,col=1)



# groupwise cosine similarity----------------------------------------------------------

cosine_dist=function(inp_mat){
  Matrix <- as.matrix(inp_mat)
  sim <- Matrix / sqrt(rowSums(Matrix * Matrix))
  cosine_sim <- sim %*% t(sim)  
  return(cosine_sim)
}
CBCL_centroid_cosine=cosine_dist(CBCL_centroidmat)
row.names(CBCL_centroid_cosine)=name_list
colnames(CBCL_centroid_cosine)=name_list

YSR_centroid_cosine=cosine_dist(YSR_centroidmat)
row.names(YSR_centroid_cosine)=name_list
colnames(YSR_centroid_cosine)=name_list

YSR_centroid_cosine=YSR_centroid_cosine[c('RBB','AB','AP','SP','TP','WD','SC','AD'),c('RBB','AB','AP','SP','TP','WD','SC','AD')]
CBCL_centroid_cosine=CBCL_centroid_cosine[c('RBB','AB','AP','SP','TP','WD','SC','AD'),c('RBB','AB','AP','SP','TP','WD','SC','AD')]

### CBCL_cosine_sim_v2.png
# pheatmap(round(CBCL_centroid_cosine,3),main='', cluster_rows = F,cluster_cols = F,display_numbers=T,fontsize = 20,fontsize_number=23,breaks=inp_breaks,filename='CBCL_cosine_sim_v2_highres_v2.png',res=600,angle_col = "0",width=9, height=9)
### YSR_cosine_sim_v2.png
# pheatmap(round(YSR_centroid_cosine,3),main='', cluster_rows = F,cluster_cols = F,display_numbers=T,fontsize = 20,fontsize_number=23,breaks=inp_breaks,filename='YSR_cosine_sim_v2_highres_v2.png',res=600,angle_col = "0",width=9, height=9)

# write.csv(round(YSR_centroid_dist,3),file="Figure11_data_YSRcosine.csv",row.names=T)
# write.csv(round(CBCL_centroid_dist,3),file="Figure11_data_CBCLcosine.csv",row.names=T)
dbscan(dist(CBCL_centroid_cosine),eps = median(dist(CBCL_centroid_cosine)),MinPts = 2)$cluster
dbscan(dist(YSR_centroid_cosine),eps = median(dist(CBCL_centroid_cosine)),MinPts = 2)$cluster

#################If we see  groupwise distance, not cosine similarity
CBCL_centroid_dist=as.matrix(dist(CBCL_centroidmat))
row.names(CBCL_centroid_dist)=name_list
colnames(CBCL_centroid_dist)=name_list

YSR_centroid_dist=as.matrix(dist(YSR_centroidmat))
row.names(YSR_centroid_dist)=name_list
colnames(YSR_centroid_dist)=name_list

YSR_centroid_dist=YSR_centroid_dist[c('RBB','AB','AP','SP','TP','WD','SC','AD'),c('RBB','AB','AP','SP','TP','WD','SC','AD')]
CBCL_centroid_dist=CBCL_centroid_dist[c('RBB','AB','AP','SP','TP','WD','SC','AD'),c('RBB','AB','AP','SP','TP','WD','SC','AD')]

# write.csv(round(YSR_centroid_dist,3),file="Figure11_data_YSRL2.csv",row.names=T)
# write.csv(round(CBCL_centroid_dist,3),file="Figure11_data_CBCLL2.csv",row.names=T)
#################

#Clustering respondent based on the relationship with item subgroups----------------------------------------------------------
cal_dist_subgroup=function(inp_position,inp_centroid_mat){
  return(sqrt(rowSums((inp_centroid_mat-matrix(rep(inp_position,8),nrow=8,ncol=2,byrow = T))**2)))    
}
dist_subgroup_mat=t(apply(CBCL_z.est_pr,MARGIN = 1,FUN = function(x) cal_dist_subgroup(x,CBCL_centroidmat)))

resp_clust_res=kmeans(dist_subgroup_mat, centers = 4,nstart = 25,iter.max = 1e8)
plot(CBCL_z.est_pr)
points(CBCL_z.est_pr[resp_clust_res$cluster==1,1],CBCL_z.est_pr[resp_clust_res$cluster==1,2],col=1,pch=20,cex=2)
points(CBCL_z.est_pr[resp_clust_res$cluster==2,1],CBCL_z.est_pr[resp_clust_res$cluster==2,2],col=2,pch=20,cex=2)
points(CBCL_z.est_pr[resp_clust_res$cluster==3,1],CBCL_z.est_pr[resp_clust_res$cluster==3,2],col=3,pch=20,cex=2)
points(CBCL_z.est_pr[resp_clust_res$cluster==4,1],CBCL_z.est_pr[resp_clust_res$cluster==4,2],col=4,pch=20,cex=2)
points(CBCL_centroidmat,pch=8,cex=5)


dist_subgroup_mat=t(apply(YSR_z.est,MARGIN = 1,FUN = function(x) cal_dist_subgroup(x,YSR_centroidmat)))
#resp_clust_res=dbscan(dist(dist_subgroup_mat),eps = median(dist(dist_subgroup_mat)),MinPts = 60)$cluster
resp_clust_res=kmeans(dist_subgroup_mat, centers = 4,nstart = 25,iter.max = 1e8)
plot(YSR_z.est)
points(YSR_z.est[resp_clust_res$cluster==1,1],YSR_z.est[resp_clust_res$cluster==1,2],col=1,pch=20,cex=2)
points(YSR_z.est[resp_clust_res$cluster==2,1],YSR_z.est[resp_clust_res$cluster==2,2],col=2,pch=20,cex=2)
points(YSR_z.est[resp_clust_res$cluster==3,1],YSR_z.est[resp_clust_res$cluster==3,2],col=3,pch=20,cex=2)
points(YSR_z.est[resp_clust_res$cluster==4,1],YSR_z.est[resp_clust_res$cluster==4,2],col=4,pch=20,cex=2)
points(YSR_centroidmat,pch=8,cex=5)



#find friendly group for each individual respondents------------------------------------------
tmp=procrustes(CBCL_w.est,YSR_w.est)
CBCL_w.est_pr=tmp$X.new #procruste, just for visulize
CBCL_z.est_pr=tmp$s*CBCL_z.est%*%tmp$R

inp_type='YSR'
{if(inp_type=='CBCL'){
  inp_w.est=CBCL_w.est_pr
  inp_z.est=CBCL_z.est_pr
}
  else{
    inp_w.est=YSR_w.est
    inp_z.est=YSR_z.est
  }
}


centroidmat=matrix(NA,ncol=2,nrow=length(data_list))
for (idx in 1:length(data_list)){
  #set target_idx  
  target_idx=data_list[[idx]]
  
  #calculate centroid or group distance, for each target_idx
  tmp_centroid=c(mean(inp_w.est[item_name %in%target_idx,1]),mean(inp_w.est[item_name %in%target_idx,2]))
  centroidmat[idx,]=tmp_centroid
}

cal_dist=function(inp_z=inp_z.est[1,],inp_cent=centroidmat,top_n=1){
  
  #cosine
  rbind(inp_z,inp_cent)
  tmp_dist=cosine_dist(rbind(inp_z,inp_cent))[2:(length(data_list)+1),1]
  
  #L2
  tmp_dist=sqrt(rowSums((matrix(rep(inp_z,length(data_list)),nrow=length(data_list),byrow = T)-inp_cent)**2))
  
  near_group=name_list[order(tmp_dist,decreasing = F)[1:top_n]]
  near_dist=tmp_dist[order(tmp_dist,decreasing = F)[1:top_n]]
  #return(list('name'=near_group,'dist'=near_dist))
  return(near_group)
}

near_item_mat=apply(X = inp_z.est,MARGIN = 1,FUN = function(x) cal_dist(inp_z=x))

indi_in_hub=order(rowSums(inp_z.est**2),decreasing = F)[1]
indi_in_out=order(rowSums(inp_z.est**2),decreasing = T)[1:2]

### indi_plot2_CBCL1_v2 , indi_plot2_CBCL3_v2
# visulization
par(mfrow=c(1,1),mar=rep(3,4))
#for(inp_idx in c(indi_in_out,indi_in_hub)){
for(inp_idx in c(1,3,16,102,250)){
  cat(sprintf('%s, %s',inp_type,inp_idx))
  plot(inp_w.est,type='n',xlim=c(-2.5,2.5),ylim=c(-3,3),xlab='',ylab='',cex.axis=axis_size)
  #points(inp_z.est,col=2)
  #points(inp_w.est)
  
  points(centroidmat,cex=2.5,pch=8,lwd=2,col=ifelse(inp_type=='YSR',2,1))
  text(centroidmat[,1],centroidmat[,2]+up_shift,name_list,lwd=2,col=ifelse(inp_type=='YSR',2,1),cex=text_cex)
  #inp_idx=3
  points(inp_z.est[inp_idx,1],inp_z.est[inp_idx,2],col='blue',cex=3,lwd=3)
  legend(x=1.25,y=3,c(inp_type,'Respondent'),col=c(ifelse(inp_type=='YSR',2,1),'blue'),pch=c(8,21),cex=1.4)
}

### indi_plot2_CBCLall_v2.png
plot(inp_w.est,type='n',xlim=c(-2.5,2.5),ylim=c(-3,3),xlab='',ylab='')
points(inp_z.est,col=alpha('blue',0.5),cex=1.2)
points(centroidmat,cex=2.5,pch=8,lwd=2,col=ifelse(inp_type=='YSR',2,1))
legend(x=1.25,y=3,c(inp_type,'Respondent'),col=c(ifelse(inp_type=='YSR',2,1),'blue'),pch=c(8,21),cex=1.4)

#YSR_near_item_mat=near_item_mat
#CBCL_near_item_mat=near_item_mat
#mean(CBCL_near_item_mat==YSR_near_item_mat)

