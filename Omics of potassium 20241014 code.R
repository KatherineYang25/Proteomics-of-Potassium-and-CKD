rm(list=ls())

library(nloptr)
library(remotes)
library(readstata13)
library(readxl)
library(data.table)
library("survival")
library(ggplot2)
library(VennDiagram)
library("metafor")
library(MASS)
library(sfsmisc)
library(ggrepel)
library(DescTools)
library(plyr)
library(survminer)
library(lme4)
library(reshape2)
library("lmerTest")
library(dplyr)


library(readstata13)
library(readxl)
library(plyr)
library(dplyr)
library(ggplot2)
library(grid)
library(data.table)
library(glmnet)
library(matrixStats)
library(gdata)
library(stringr)
library(pROC)
library(survival)
library(compareC)
library(Hmisc) 
library(calibrate)
library(ggrepel)
library(scales)
library(stringi)
library(gridExtra)
library(DescTools)
library(VennDiagram)
library(metafor)
library(boot)
library(boot.pval)
library(mediation)
library(MASS)
library(lmtest)
library(lme4)
library(survival)


################################################################################################################
### 1. data prep
########################
# proteomics
annot_soma=read.delim(paste0(soma.dir,"soma_visit_3_annot_ANML_SMP_updated.txt"), header = TRUE, sep = "\t", stringsAsFactors =FALSE)

soma_v3_SMP=read.csv(paste0(in.dir,"soma_v3_SMP.csv"))
soma_v3_SMP=soma_v3_SMP[ ,-which(names(soma_v3_SMP) %in% c("X"))]
soma_v3_SMP_d=read.csv(paste0(in.dir,"soma_v3_SMP_d.csv"))
soma_v3_SMP_d=soma_v3_SMP_d[ ,-which(names(soma_v3_SMP_d) %in% c("X"))]
soma_v3_SMP_v=read.csv(paste0(in.dir,"soma_v3_SMP_v.csv"))
soma_v3_SMP_v=soma_v3_SMP_v[ ,-which(names(soma_v3_SMP_v) %in% c("X"))]

pro.soma<-annot_soma[annot_soma$flag2==0,"seqid_in_sample"]
length(pro.soma)
# 4955

# factor variables
ds=c("soma_v3_SMP","soma_v3_SMP_d","soma_v3_SMP_v")
dm=c(rep("pro.soma",3))

for(j in 1:length(ds)) {
  d=get(ds[j]) 
  met=get(dm[j])
  
  colnames(d)[colnames(d)=="age_v3"] <- "age"
  colnames(d)[colnames(d)=="egfrcr_v3"] <- "egfrcr"
  colnames(d)[colnames(d)=="cigt31"] <- "cigt"
  colnames(d)[colnames(d)=="bmi32"] <- "bmi"
  colnames(d)[colnames(d)=="sprt_i31"] <- "sprt"
  colnames(d)[colnames(d)=="hypert35"] <- "hypertens"
  colnames(d)[colnames(d)=="ethanl32"] <- "ethanl"
  colnames(d)[colnames(d)=="diabts34"] <- "diabetes"
  
  d$black<-as.numeric(d$racegrp=="B")
  d$white<-as.numeric(d$racegrp=="W")
  d$asian<-as.numeric(d$racegrp=="A")
  d$indian<-as.numeric(d$racegrp=="I")
  d$female<-as.numeric(d$gender=="F")
  #d$racecen<-as.factor(paste(d$center,"_",d$racegrp,sep=""))
  d$elevel02=as.factor(d$elevel02)
  d$center=as.factor(d$center)
  d$cigt=as.factor(d$cigt)
  d <- within(d, center <- relevel(center, ref = "J"))
  #d <- within(d, racecen <- relevel(racecen, ref = "J_B"))
  d <- within(d, elevel02 <- relevel(elevel02, ref = "1"))
  d <- within(d, cigt <- relevel(cigt, ref = "3"))
  
  # dummy variable for racecen center 
  d<-data.table(d)
  #d[, `:=` (racecenF_B = (racecen=="F_B")*1, racecenF_W = (racecen=="F_W")*1,racecenJ_B = (racecen=="J_B")*1,racecenW_W = (racecen=="W_W")*1,racecenM_W = (racecen=="M_W")*1)]
  d[, `:=` (racegrp_B = (racegrp=="B")*1, racegrp_A = (racegrp=="A")*1,racegrp_I = (racegrp=="I")*1,racegrp_W = (racegrp=="W")*1)]
  d[, `:=` (centerF = (center=="F")*1, centerJ = (center=="J")*1,centerM = (center=="M")*1,centerW = (center=="W")*1)]
  d[, `:=` (elevel02_1 = (elevel02=="1")*1, elevel02_2 = (elevel02=="2")*1,elevel02_3 = (elevel02=="3")*1)]
  d[, `:=` (cigt_current = (cigt=="1")*1, cigt_former = (cigt=="2")*1,cigt_never = (cigt=="3")*1)]
  d<-data.frame(d)
  
  # standardize pota and density
  d[,"pota_orig"]=d[,"pota"]
  d[,"density_orig"]=d[,"density"]
  sd=sd(d[,"pota"],na.rm=T)
  d[,"pota"]=d[,"pota"]/sd
  sd=sd(d[,"density"],na.rm=T)
  d[,"density"]=d[,"density"]/sd
  
  assign(ds[j],d)    
}

# merge in separate ckd outcomes
f.ckd=paste(ckd.dir,"kidney_outcomes_separated.dta",sep="")
d.ckd=read.dta13(f.ckd)
d.ckd=d.ckd[ , -which(names(d.ckd) %in% c("inc_ckd_def2_v3","fu_ckd_def2_v3"))]
soma_v3_SMP=merge(soma_v3_SMP,d.ckd,by.x="SampleId",by.y="id")
dim(soma_v3_SMP)
soma_v3_SMP_d=merge(soma_v3_SMP_d,d.ckd,by.x="SampleId",by.y="id")
dim(soma_v3_SMP_d)
soma_v3_SMP_v=merge(soma_v3_SMP_v,d.ckd,by.x="SampleId",by.y="id")
dim(soma_v3_SMP_v)
# N=10197

#eGFR<60 with >=25% decline
table(soma_v3_SMP$inc_ckd_def1_v3) #n_event=1698
# ckd follow up median
summary(soma_v3_SMP$fu_ckd_def1_v3/365.25)
# median 20.55
summary(soma_v3_SMP_d$fu_ckd_def1_v3/365.25)
# median 20.77
summary(soma_v3_SMP_v$fu_ckd_def1_v3/365.25)
# median 20.27

#hospitalization or death related to CKD
table(soma_v3_SMP$hosp_dth) #n_event=2561
summary(soma_v3_SMP$fu_hosp_dth_v3/365.25) #median 21.34

#eskd
table(soma_v3_SMP$eskd) #n_event=203
summary(soma_v3_SMP$fu_eskd_v3/365.25) #median 23.08


#######################################################################################
### 2. table 1 by discovery and validation
ds=c("soma_v3_SMP","soma_v3_SMP_d","soma_v3_SMP_v")

h<-data.frame(varname=character(),
              mean_sd_soma_v3_SMP=character(),mean_sd_soma_v3_SMP_d=character(),mean_sd_soma_v3_SMP_v=character(),
              n=integer(),
              stringsAsFactors=FALSE)


for(j in 1:length(ds)) {
  line=1
  d=get(ds[j])
  
  h[line,"varname"]="N"
  h[line,paste0("mean_sd_",ds[j])]=nrow(d)
  line=line+1
  
  covar=c("female","black","white","asian","indian","centerF","centerJ","centerM","centerW",
          "elevel02_1","elevel02_2","elevel02_3",
          "cigt_current","cigt_former","cigt_never","diabetes","hypertens")
  for(i in 1:length(covar)) {
    h[line,"varname"]=covar[i]
    n1=length(which(d[,covar[i]]==1))
    n2=length(which(!is.na(d[,covar[i]])))
    pc=round(n1/n2*100,2)
    h[line,paste0("mean_sd_",ds[j])]=paste0(n1," (",pc,"%)")
    line=line+1  
  }
  
  covar=c("age","bmi","egfrcr","tcal","sprt","ethanl","pota_orig","density_orig")
  for(i in 1:length(covar)) {
    h[line,"varname"]=covar[i]
    n1=mean(d[,covar[i]],na.rm =T)
    n2=sd(d[,covar[i]],na.rm =T)
    h[line,paste0("mean_sd_",ds[j])]=paste0(round(n1,2)," (",round(n2,2),")")
    line=line+1  
  }
  
  covar=c("N","age","female","black","white","asian","indian","centerF","centerJ","centerM","centerW",
          "elevel02_1","elevel02_2","elevel02_3","cigt_current","cigt_former","cigt_never","diabetes","hypertens",
          "bmi","egfrcr","tcal","sprt","ethanl","pota_orig","density_orig")
  for(i in 1:length(covar)) {
    h[h$varname==covar[i],"n"]=i
  }
}

h=h[order(h$n),]
h=h[ , -which(names(h) %in% c("n"))]
write.table(h,file=paste(out.dir,"table 1.txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)


### 2. table 1 by CKD outcome
#cohort with CKD
ckd_1 <- subset(soma_v3_SMP, inc_ckd_def2_v3 == 1)
#cohort w/o CKD
ckd_0 <- subset(soma_v3_SMP, inc_ckd_def2_v3 == 0)


ds=c("soma_v3_SMP","ckd_1","ckd_0")

h<-data.frame(varname=character(),
              mean_sd_soma_v3_SMP=character(),mean_sd_ckd_1=character(),mean_sd_ckd_0=character(),
              n=integer(),
              stringsAsFactors=FALSE)


for(j in 1:length(ds)) {
  line=1
  d=get(ds[j])
  
  h[line,"varname"]="N"
  h[line,paste0("mean_sd_",ds[j])]=nrow(d)
  line=line+1
  
  covar=c("female","black","white","asian","indian","centerF","centerJ","centerM","centerW",
          "elevel02_1","elevel02_2","elevel02_3",
          "cigt_current","cigt_former","cigt_never","diabetes","hypertens")
  for(i in 1:length(covar)) {
    h[line,"varname"]=covar[i]
    n1=length(which(d[,covar[i]]==1))
    n2=length(which(!is.na(d[,covar[i]])))
    pc=round(n1/n2*100,2)
    h[line,paste0("mean_sd_",ds[j])]=paste0(n1," (",pc,"%)")
    line=line+1  
  }
  
  covar=c("age","bmi","egfrcr","tcal","sprt","ethanl","pota_orig","density_orig")
  for(i in 1:length(covar)) {
    h[line,"varname"]=covar[i]
    n1=mean(d[,covar[i]],na.rm =T)
    n2=sd(d[,covar[i]],na.rm =T)
    h[line,paste0("mean_sd_",ds[j])]=paste0(round(n1,2)," (",round(n2,2),")")
    line=line+1  
  }
  
  covar=c("N","age","female","black","white","asian","indian","centerF","centerJ","centerM","centerW",
          "elevel02_1","elevel02_2","elevel02_3","cigt_current","cigt_former","cigt_never","diabetes","hypertens",
          "bmi","egfrcr","tcal","sprt","ethanl","pota_orig","density_orig")
  for(i in 1:length(covar)) {
    h[h$varname==covar[i],"n"]=i
  }
}

h=h[order(h$n),]
h=h[ , -which(names(h) %in% c("n"))]
write.table(h,file=paste(out.dir,"table 1-ckd.txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)



############################
### 2) linear models
###density
# discovery
exposure<-c("density")
Model<-c("+age+female+racegrp+center+elevel02+cigt+bmi+egfrcr+tcal+sprt+ethanl")
model<-c(1)


d=soma_v3_SMP_d

for(k in 1:length(exposure)) {
  for(m in 1:length(Model)) {
    h<-data.frame(mid=character(),beta=double(),se=double(),p=double(),stringsAsFactors=FALSE)
    
    for(i in 1:length(pro.soma)) {
      h[i,"mid"]<-pro.soma[i]
      fmla<-as.formula(paste(pro.soma[i],"~",exposure[k],Model[m],sep=""))
      f<-lm(fmla,data=d, na.action=na.exclude)
      h[i,"beta"]<-summary(f)$coefficient[2,1]
      h[i,"se"]<-summary(f)$coefficient[2,2]
      h[i,"p"]<-summary(f)$coefficient[2,4]
    }
    
    h$p_FDR<-p.adjust(h$p, method = "fdr")
    h<-merge(annot_soma[,c("seqid_in_sample","uniprot_id","entrezgeneid","entrezgenesymbol","targetfullname")],h,by.x="seqid_in_sample",by.y="mid")
    h<-h[order(h$p),]
    write.table(h,file=paste(out.dir,"linear/Linear_",exposure[k],"_Model",model[m],"_discovery2.txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
  }
}

k=1
m=1

g<-read.delim(paste0(out.dir,"linear/Linear_",exposure[k],"_Model",model[m],"_discovery2.txt"), header = TRUE, sep = "\t", stringsAsFactors =FALSE)
length(which(g$p<0.05))
# 643

length(which(g$p_FDR<0.05))
# 147

length(which(g$p<0.05/nrow(g)))
# 46

nrow(g) #n=4955 total number of proteins


# step-wise approach - using significant proteins (FDR<0.05) using density to move forward to validation
exposure<-c("density")
Model<-c("+age+female+racegrp+center+elevel02+cigt+bmi+egfrcr+tcal+sprt+ethanl")
model<-c(1)

d=soma_v3_SMP_v

for(k in 1:length(exposure)) {
  for(m in 1:length(Model)) {
    g<-read.delim(paste0(out.dir,"linear/Linear_",exposure[k],"_Model",model[m],"_discovery2.txt"), header = TRUE, sep = "\t", stringsAsFactors =FALSE)
    seq_ids_sel=g[g$p_FDR<0.05,"seqid_in_sample"]
    print(length(seq_ids_sel))
    
    h<-data.frame(mid=character(),beta=double(),se=double(),p=double(),stringsAsFactors=FALSE)
    
    for(i in 1:length(seq_ids_sel)) {
      h[i,"mid"]<-seq_ids_sel[i]
      fmla<-as.formula(paste(seq_ids_sel[i],"~",exposure[k],Model[m],sep=""))
      f<-lm(fmla,data=d, na.action=na.exclude)
      h[i,"beta"]<-summary(f)$coefficient[2,1]
      h[i,"se"]<-summary(f)$coefficient[2,2]
      h[i,"p"]<-summary(f)$coefficient[2,4]
    }
    
    h$p_FDR<-p.adjust(h$p, method = "fdr")
    h<-merge(annot_soma[,c("seqid_in_sample","uniprot_id","entrezgeneid","entrezgenesymbol","targetfullname")],h,by.x="seqid_in_sample",by.y="mid")
    h<-h[order(h$p),]
    write.table(h,file=paste(out.dir,"linear/Linear_",exposure[k],"_Model",model[m],"_validation2.txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
  }
}


g<-read.delim(paste0(out.dir,"linear/Linear_",exposure[k],"_Model",model[m],"_validation2.txt"), header = TRUE, sep = "\t", stringsAsFactors =FALSE)

length(which(g$p<0.05))
# 96

length(which(g$p_FDR<0.05))
# 85

length(which(g$p<0.05/nrow(g)))
# 38

nrow(g) #n=147


###############
## heatmap
exposure<-c("density")
Model<-c("+age+female+racecen+elevel02+cigt+bmi+egfrcr+tcal+sprt+ethanl")
model<-c(1)

m=1
k=1
h<-read.delim(paste0(out.dir,"linear/Linear_",exposure[k],"_Model",model[m],"_validation.txt"), header = TRUE, sep = "\t", stringsAsFactors =FALSE)
d=soma_v3_SMP
h=h[!is.na(h$p),]
seq_ids_en=h[h$p<0.05/147,"seqid_in_sample"]
length(seq_ids_en)
# N=38
seq_ids_en2=h[h$p_FDR<0.05,"seqid_in_sample"]
length(seq_ids_en2)
#N=85
cormat=cor(d[,seq_ids_en],method="spearman")

get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
reorder_cormat <- function(cormat){
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}
# Reorder the correlation matrix
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
lower_tri <- get_lower_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
melted_cormat$i=1:nrow(melted_cormat)
melted_cormat=merge(melted_cormat,annot[,c("seqid_in_sample","entrezgenesymbol")],by.x="Var1",by.y="seqid_in_sample")
colnames(melted_cormat)[colnames(melted_cormat)=="entrezgenesymbol"] <- "entrezgenesymbol1"
melted_cormat=merge(melted_cormat,annot[,c("seqid_in_sample","entrezgenesymbol")],by.x="Var2",by.y="seqid_in_sample")
colnames(melted_cormat)[colnames(melted_cormat)=="entrezgenesymbol"] <- "entrezgenesymbol2"
melted_cormat=melted_cormat[order(melted_cormat$i),]

# Create a ggheatmap
out.plot=paste(out.dir,"plots/heatmaps/","Heatmap_linear_density_model1_validation.tiff",sep="")
tiff(filename = out.plot,width=4800, height=4800, res=400)

ggheatmap <- ggplot(melted_cormat, aes(x=reorder(entrezgenesymbol2,i), y=reorder(entrezgenesymbol1,i), fill = value))+
  geom_tile()+
  geom_text(aes(x=reorder(entrezgenesymbol2,i), y=reorder(entrezgenesymbol1,i), label = sprintf("%.2f",value)), color = "black", size = 5.5) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="") +
  theme_minimal()+ # minimal theme
  coord_fixed() + 
  theme(axis.text.x=element_text(color = "black", size=18, angle=90, vjust=.8, hjust=0.8)) + 
  theme(axis.text.y=element_text(color = "black", size=18, angle=0, vjust=.8, hjust=0.8)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = c(0.2, 0.8)) +
  ylab("") +
  xlab("") +
  scale_y_discrete(position = "right") 
# Print the heatmap
print(ggheatmap)

dev.off()


# race stratified
ds=c("soma_v3_SMP_d","soma_v3_SMP_v")
da=c(rep("annot_soma",2))
exposure<-c("pota","density")
model<-c(1)
Model=c("+age+female+center+elevel02+cigt+bmi+egfrcr+tcal+sprt+ethanl")
dss=c("black","white")



for(j in 1:length(ds)) {
  d=get(ds[j])
  annot=get(da[j])
  
  for(k in 1:length(exposure)) {
    for(m in 1:length(Model)) {
      t=read.delim(paste(out.dir,"linear/Linear_",exposure[k],"_Model",model[m],"_soma_v3_SMP_d.txt",sep=""), header = TRUE, sep = "\t", stringsAsFactors =FALSE)
      t=t[!is.na(t$p),]
      t=t[t$p<0.05/nrow(t),]
      met_sel=t$seqid_in_sample
      print(length(met_sel))
      
      h<-data.frame(mid=character(),
                    N_white=integer(),beta_white=double(),se_white=double(),p_white=double(),
                    N_black=integer(),beta_black=double(),se_black=double(),p_black=double(),
                    p_inter=double(),
                    stringsAsFactors=FALSE)
      
  for(jj in 1:length(dss)) {
    dx=d[d[,dss[jj]]==1,]
    print(dim(dx))
 
        for(i in 1:length(met_sel)) {
          h[i,"mid"]<-met_sel[i]
          h[i,paste0("N_",dss[jj])]<-nrow(dx)
          fmla<-as.formula(paste(met_sel[i],"~",exposure[k],Model[m],sep=""))
          try.test<-try(f<-lm(fmla,data=dx, na.action=na.exclude))
          if(class(try.test)!="try-error") {
            h[i,paste0("beta_",dss[jj])]<-summary(f)$coefficient[2,1]
            h[i,paste0("se_",dss[jj])]<-summary(f)$coefficient[2,2]
            h[i,paste0("p_",dss[jj])]<-summary(f)$coefficient[2,4]
          }
          fmla<-as.formula(paste(met_sel[i],"~",exposure[k],"*",dss[jj],Model[m],sep=""))
          f<-lm(fmla,data=d, na.action=na.exclude)
          h[i,"p_inter"]<-summary(f)$coefficient[nrow(summary(f)$coefficient),4]
        }
  }
      
      h<-merge(annot[,c("seqid_in_sample","uniprot_id","entrezgeneid","entrezgenesymbol","targetfullname")],h,by.x="seqid_in_sample",by.y="mid")
      write.table(h,file=paste(out.dir,"linear/Linear_stratified_",exposure[k],"_Model",model[m],"_",ds[j],"_",dss[jj],".txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
    }
  }
}



## c stat compare using hits from cox model, see if they (individually and together) can improve the predication of incident CKD
ds=c("soma_v3_SMP")
dss=c("")
model<-c(1)
Model=c("age+female+racegrp+center+elevel02+cigt+bmi+egfrcr+tcal+sprt+ethanl+diabetes+hypertens")
outcome<-c("ckd_def2_v3")
event=c("inc_ckd_def2_v3")
fu=c("fu_ckd_def2_v3")  


for(j in 1:length(ds)) {
  d=get(paste0(ds[j]))
  annot=annot_soma
  t<-read.delim(paste0(out.dir, "cox/Cox_ckd_def2_v3_Model2_fdr1.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  t=t[t$p_FDR<0.05,]
  met_sel=t$seqid_in_sample
  met_sel=grep("SeqId_",met_sel,value=T)
  print(length(met_sel))
  
  for(k in 1:length(outcome)) {
    for(m in 1:length(model)) {
      h<-data.frame(mid=character(),
                    Cstat=double(),Cstat_lci=double(),Cstat_uci=double(),
                    Cstat_noprotein=double(),Cstat_noprotein_lci=double(),Cstat_noprotein_uci=double(),
                    Cstat_diff=double(),Cstat_diff_lci=double(),Cstat_diff_uci=double(),
                    p_Cstat_diff=double(),
                    stringsAsFactors=FALSE)
      
      # covariates only
      fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~",Model[m],sep=""))
      f<-coxph(fmla,data=d, na.action=na.exclude)
      c0<-concordance(f)$concordance
      c0_lci=concordance(f)$concordance-1.96*sqrt(concordance(f)$var)
      c0_uci=concordance(f)$concordance+1.96*sqrt(concordance(f)$var)     
      d[,paste0("new.pred0")] = predict(f, d)
      
      # single protein
      for(i in 1:length(met_sel)) {
        h[i,"mid"]<-met_sel[i]
        fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~",met_sel[i],"+",Model[m],sep=""))
        f<-coxph(fmla,data=d, na.action=na.exclude)
        d[,paste0("new.pred1")] = predict(f, d)
        h[i,"Cstat"]<-concordance(f)$concordance
        h[i,"Cstat_lci"]<-concordance(f)$concordance-1.96*sqrt(concordance(f)$var)
        h[i,"Cstat_uci"]<-concordance(f)$concordance+1.96*sqrt(concordance(f)$var)
        h[i,"Cstat_noprotein"]<-c0
        h[i,"Cstat_noprotein_lci"]<-c0_lci
        h[i,"Cstat_noprotein_uci"]<-c0_uci
        compare<-compareC(d[,fu[k]], d[,event[k]], d[,paste0("new.pred0")],  d[,paste0("new.pred1")])
        h[i,"Cstat_diff"]<-compare$est.diff_c
        h[i,"Cstat_diff_lci"]<-compare$est.diff_c-1.96*sqrt(compare$est.vardiff_c)
        h[i,"Cstat_diff_uci"]<-compare$est.diff_c+1.96*sqrt(compare$est.vardiff_c)
        h[i,"p_Cstat_diff"]<-compare$pval                         
      }
      
      # all proteins together
      i=i+1
      h[i,"mid"]<-"all proteins together"
      fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~",paste(met_sel,collapse="+"),"+",Model[m],sep=""))
      f<-coxph(fmla,data=d, na.action=na.exclude)
      d[,paste0("new.pred1")] = predict(f, d)
      h[i,"Cstat"]<-concordance(f)$concordance
      h[i,"Cstat_lci"]<-concordance(f)$concordance-1.96*sqrt(concordance(f)$var)
      h[i,"Cstat_uci"]<-concordance(f)$concordance+1.96*sqrt(concordance(f)$var)
      h[i,"Cstat_noprotein"]<-c0
      h[i,"Cstat_noprotein_lci"]<-c0_lci
      h[i,"Cstat_noprotein_uci"]<-c0_uci
      compare<-compareC(d[,fu[k]], d[,event[k]], d[,paste0("new.pred0")],  d[,paste0("new.pred1")])
      h[i,"Cstat_diff"]<-compare$est.diff_c
      h[i,"Cstat_diff_lci"]<-compare$est.diff_c-1.96*sqrt(compare$est.vardiff_c)
      h[i,"Cstat_diff_uci"]<-compare$est.diff_c+1.96*sqrt(compare$est.vardiff_c)
      h[i,"p_Cstat_diff"]<-compare$pval  
      
      h<-merge(annot[,c("seqid_in_sample","uniprot_id","entrezgeneid","entrezgenesymbol","targetfullname")],h,by.x="seqid_in_sample",by.y="mid",all.y=T)
      write.table(h,file=paste(out.dir,"cox/Cstats to predict ",outcome[k],"_model",model[m],dss[j],".txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
    }
  }
}



#######################################################################################
### 5) plots for linear
## volcano plots

exposure<-c("pota")
Model<-c("+age+female+racecen+elevel02+cigt+bmi+egfrcr+tcal+sprt+ethanl")
model<-c(1)

k=1
m=1

for(k in 1:length(exposure)) {
    for(m in 1:length(Model)) {
      h=read.delim(paste(out.dir,"linear/Linear_",exposure[k],"_Model",model[m],"_discovery.txt",sep=""), header = TRUE, sep = "\t", stringsAsFactors =FALSE)
      sig=h[h$p_FDR<0.05,"seqid_in_sample"]
      h$logp=-log10(h$p)
      h=h[!is.na(h$p),]
      h[h$seqid_in_sample %in% sig,"c"]=1
      h[is.na(h$c),"c"]=2
      h$c=as.factor(h$c)
      cols <- c("1" = "red", "2"= "gray50")
      h[h$c=="1","lab"]=h[h$c=="1","entrezgenesymbol"]

      
      out.plot=paste(out.dir,"plots/linear/","Volcano_linear_",exposure[k],"_Model",model[m],".tiff",sep="")
      tiff(filename = out.plot,width = 1000, height = 1000)
      print(
        ggplot(h, aes(beta, logp, label = lab,color=c)) +
          scale_color_manual(values = cols) +
          guides(colour=FALSE) +
          geom_point(size=4) +
          geom_text_repel(color="black",size=4.5,max.overlaps=35) +
          theme_classic(base_size = 15) +
          geom_vline(xintercept=0,col="black",lty="dashed") +
          ylab("-log10(p)") +
          xlab("Beta discovery")
      )
      dev.off()
  }
}



############################################################################
#scatter plot - linear - compare discovery vs validation results - density
## 5.1) linear
exposure<-c("density")
Model<-c("+age+female+racegrp+center+elevel02+cigt+bmi+egfrcr+tcal+sprt+ethanl")
model<-c(1)

k=1
m=1

for(k in 1:length(exposure)) {
  for(m in 1:length(model)) {
    
    # beta scatter plots
    d<-read.delim(paste0(out.dir,"linear/Linear_",exposure[k],"_Model",model[m],"_discovery2.txt"), header = TRUE, sep = "\t", stringsAsFactors =FALSE)
    v<-read.delim(paste0(out.dir,"linear/Linear_",exposure[k],"_Model",model[m],"_validation2.txt"), header = TRUE, sep = "\t", stringsAsFactors =FALSE)
    
    h=merge(d,v,by="seqid_in_sample",all = TRUE)
    h=h[!is.na(h$p.x) & !is.na(h$p.y),]
    h<-h[order(h$p.x),]
    # d x v y
    
    h[h$p_FDR.x<0.05 & h$p_FDR.y<0.05,"Threshold"]="FDR<0.05 in replication"
    h[h$p_FDR.x < 0.05 & h$p_FDR.y >= 0.05 & h$p_FDR.y <= 1, "Threshold"] <- "FDR<0.05 in discovery"

    h$Threshold=as.factor(h$Threshold)
    cols <- c("FDR<0.05 in replication" = "orange","FDR<0.05 in discovery" = "red","v p_FDR<0.05" = "purple","Discovery p_FDR>=0.05" = "gray80")
    h[h$Threshold=="FDR<0.05 in replication","lab"]=h[h$Threshold=="FDR<0.05 in replication","entrezgenesymbol.x"]
    h[h$Threshold=="FDR<0.05 in discovery","lab"]=h[h$Threshold=="FDR<0.05 in discovery","entrezgenesymbol.x"]
    
    
    
    # plot
    out.plot=paste(out.dir,"plots/Scatter/Scatter_Beta_r vs. d_",exposure[k],"_Model",model[m]," (by p_FDR)_2.tiff",sep="")  
    tiff(filename = out.plot, width = 1500, height = 1000)
    print(
      ggplot(h, aes(beta.x, beta.y, color = Threshold, label = lab)) +
        scale_color_manual(values = cols) +
        geom_text_repel(color = "black", max.overlaps = 25, size = 5, box.padding = 0.2) +
        geom_point(size = 3.5) +  
        theme_classic(base_size = 15) +
        xlab("Beta discovery") +
        ylab("Beta replication") +
        geom_abline(slope = 1, intercept = 0, color = "black") +  # Diagonal line
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +  
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
        theme(
          legend.text = element_text(size = 15),  # Increased legend text size
          legend.title = element_text(size = 17)  # Increased legend title size
        )
    )
    dev.off()
  }
}

#scatter plot - linear - compare discovery vs validation results - density - only display validation results
## 5.1) linear
exposure<-c("density")
Model<-c("+age+female+racegrp+center+elevel02+cigt+bmi+egfrcr+tcal+sprt+ethanl")
model<-c(1)

k=1
m=1

for(k in 1:length(exposure)) {
  for(m in 1:length(model)) {
    
    # beta scatter plots
    d<-read.delim(paste0(out.dir,"linear/Linear_",exposure[k],"_Model",model[m],"_discovery2.txt"), header = TRUE, sep = "\t", stringsAsFactors =FALSE)
    v<-read.delim(paste0(out.dir,"linear/Linear_",exposure[k],"_Model",model[m],"_validation2.txt"), header = TRUE, sep = "\t", stringsAsFactors =FALSE)
    
    h=merge(d,v,by="seqid_in_sample",all = TRUE)
    h=h[!is.na(h$p.x) & !is.na(h$p.y),]
    h<-h[order(h$p.x),]
    # d x v y
    
    # Filter for FDR < 0.05 in replication
    h <- h[h$p_FDR.y < 0.05, ]
    
    h[h$p_FDR.x<0.05 & h$p_FDR.y<0.05,"Threshold"]="FDR<0.05 in replication"
    h[h$p_FDR.x < 0.05 & h$p_FDR.y >= 0.05 & h$p_FDR.y <= 1, "Threshold"] <- "FDR<0.05 in discovery"

    h$Threshold=as.factor(h$Threshold)
    cols <- c("FDR<0.05 in replication" = "orange","FDR<0.05 in discovery" = "red","v p_FDR<0.05" = "purple","Discovery p_FDR>=0.05" = "gray80")
    h[h$Threshold=="FDR<0.05 in replication","lab"]=h[h$Threshold=="FDR<0.05 in replication","entrezgenesymbol.x"]

    
    
    
    # plot
    out.plot=paste(out.dir,"plots/Scatter/Scatter_Beta_r vs. d_",exposure[k],"_Model",model[m]," (by p_FDR)_v.tiff",sep="")  
    tiff(filename = out.plot, width = 1500, height = 1000)
    print(
      ggplot(h, aes(beta.x, beta.y, color = Threshold, label = lab)) +
        scale_color_manual(values = cols) +
        geom_text_repel(color = "black", max.overlaps = 25, size = 5, box.padding = 0.2) +
        geom_point(size = 3.5) +  
        theme_classic(base_size = 15) +
        xlab("Beta discovery") +
        ylab("Beta replication") +
        geom_abline(slope = 1, intercept = 0, color = "black") +  # Diagonal line
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +  
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
        theme(
          legend.text = element_text(size = 15),  # Increased legend text size
          legend.title = element_text(size = 17)  # Increased legend title size
        )
    )
    dev.off()
  }
}




#volcano plot - only for discovery - density
exposure<-c("density")
Model<-c("+age+female+racegrp+center+elevel02+cigt+bmi+egfrcr+tcal+sprt+ethanl")
model<-c(1)

k=1
m=1

for(k in 1:length(exposure)) {
  for(m in 1:length(Model)) {
    h=read.delim(paste(out.dir,"linear/Linear_",exposure[k],"_Model",model[m],"_discovery2.txt",sep=""), header = TRUE, sep = "\t", stringsAsFactors =FALSE)
    sig=h[h$p_FDR<0.05,"seqid_in_sample"]
    h$logp=-log10(h$p)
    h=h[!is.na(h$p),]
    h[h$seqid_in_sample %in% sig,"Threshold"]="FDR<0.05 in discovery"
    h[is.na(h$Threshold),"Threshold"]="Not statistically significant"
    h$Threshold=as.factor(h$Threshold)
    cols <- c("FDR<0.05 in discovery" = "red", "2"= "gray50")
    h[h$Threshold=="FDR<0.05 in discovery","lab"]=h[h$Threshold=="FDR<0.05 in discovery","entrezgenesymbol"]
    
    
    out.plot=paste(out.dir,"plots/linear/","Volcano_linear_",exposure[k],"_Model",model[m],"_2.tiff",sep="")
    tiff(filename = out.plot,width = 1000, height = 1000)
    print(
      ggplot(h, aes(beta, logp, label = lab,color=Threshold)) +
        scale_color_manual(values = cols) +
        guides(colour=FALSE) + 
        geom_point(size=4) +
        geom_text_repel(color="black",size=4.5,max.overlaps=35) +
        theme_classic(base_size = 15) +
        geom_vline(xintercept=0,col="black",lty="dashed") +
        xlim(-0.1, 0.1) +  
        ylab("-log10(p)") +
        xlab("Beta discovery")
    )
    dev.off()
  }
}


########################################################################################################################
########################################################################################################################
#cox models-incident CKD
#stratified
d_nonwhite=soma_v3_SMP[soma_v3_SMP$white==0,] #n=1946
d_white=soma_v3_SMP[soma_v3_SMP$white==1,] #n=8248
d_female=soma_v3_SMP[soma_v3_SMP$female==1,] #n=5583
d_male=soma_v3_SMP[soma_v3_SMP$female==0,] #n=4611
d_diabetes=soma_v3_SMP[soma_v3_SMP$diabetes==1,] #n=1511
d_nodiabetes=soma_v3_SMP[soma_v3_SMP$diabetes==0,] #n=8683
d_hypertens=soma_v3_SMP[soma_v3_SMP$hypertens==1,] #n=4018
d_nohypertens=soma_v3_SMP[soma_v3_SMP$hypertens==0,] #n=6176
#sensitivity
d_5yr_less=soma_v3_SMP[soma_v3_SMP$fu_ckd_def1_v3/365.25<5, ] # n=588, case=150 
d_5yr_above=soma_v3_SMP[soma_v3_SMP$fu_ckd_def1_v3/365.25>=5, ] # n=9606, case=1548
table(d_5yr_less$inc_ckd_def1_v3)
table(d_5yr_above$inc_ckd_def1_v3)


exposure<-c("density")
Model<-c("+age+female+racecen+elevel02+cigt+bmi+egfrcr+tcal+sprt+ethanl")
model<-c(1)

m=1
k=1
h<-read.delim(paste0(out.dir,"linear/Linear_",exposure[k],"_Model",model[m],"_validation.txt"), header = TRUE, sep = "\t", stringsAsFactors =FALSE)
d=soma_v3_SMP
h=h[!is.na(h$p),]
seq_ids_en=h[h$p<0.05/147,"seqid_in_sample"]
length(seq_ids_en)
# N=38
seq_ids_en2=h[h$p_FDR<0.05,"seqid_in_sample"]
length(seq_ids_en2)
#N=85


### 4) cox models
ds<-c("soma_v3_SMP","d_nonwhite","d_white","d_female","d_male","d_diabetes","d_nodiabetes","d_hypertens","d_nohypertens","d_5yr_less","d_5yr_above")
dss<-c("","_nonwhite","_white","_female","_male","_diabetes","_nodiabetes","_hypertens","_nohypertens","_5yr_less","_5yr_above")
outcome<-c("ckd_def1_v3", "hosp_death","ESKD", "ckd_def2_v3")
event=c("inc_ckd_def1_v3", "hosp_dth","eskd", "inc_ckd_def2_v3")
fu=c("fu_ckd_def1_v3","fu_hosp_dth_v3","fu_eskd_v3","fu_ckd_def2_v3")
Model<-c("+age+female+racegrp+center+elevel02+cigt+bmi+egfrcr+tcal+sprt+ethanl+diabetes+hypertens")

model<-c(1)

#using 38 from BFR
length(seq_ids_en)
# N=38

for(k in 1:length(outcome)) {  
  for(m in 1:length(Model)) {
    h<-data.frame(seq_id=character(),beta=double(),se=double(),p=double(),HR=double(),
                  lci=double(),uci=double(),stringsAsFactors=FALSE)
    
    for(i in 1:length(seq_ids_en)) {
      h[i,"seq_id"]<-seq_ids_en[i]
      
      fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~",seq_ids_en[i],Model[m],sep=""))
      f<-coxph(fmla,data=d, na.action=na.exclude)
      h[i,"beta"]<-summary(f)$coefficient[1,1]
      h[i,"se"]<-summary(f)$coefficient[1,3]
      h[i,"p"]<-summary(f)$coefficient[1,5]
      h[i,"HR"]<-summary(f)$coefficient[1,2]
      h[i,"lci"]<-exp(confint(f)[1,1])
      h[i,"uci"]<-exp(confint(f)[1,2])
    }
    
    h$p_FDR<-p.adjust(h$p, method = "fdr")
    write.table(h,file=paste(out.dir,"cox/Cox_",outcome[k],"_Model2",model[m],".txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
  }
}

length(which(h$p<0.05))
# 21


for(j in 1:length(ds)) {
  d=get(ds[j])
  print(dim(d))
  
  for(k in 1:length(outcome)) {
    for(m in 1:length(Model)) {
      h<-data.frame(seq_id=character(),beta=double(),se=double(),z=double(),p=double(),HR=double(),
                    lci=double(),uci=double(),
                    Cstat=double(),
                    p_lrtest_female=double(),p_lrtest_white=double(), p_lrtest_diabetes=double(), p_lrtest_hypertens=double(),
                    stringsAsFactors=FALSE)
      
      for(i in 1:length(seq_ids_en)) {
        h[i,"seq_id"]<-seq_ids_en[i]
        
        fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~",seq_ids_en[i],"+",Model[m],sep=""))
        f<-coxph(fmla,data=d, na.action=na.exclude)
        h[i,"beta"]<-summary(f)$coefficient[1,1]
        h[i,"se"]<-summary(f)$coefficient[1,3]
        h[i,"z"]<-summary(f)$coefficient[1,4]
        h[i,"p"]<-summary(f)$coefficient[1,5]
        h[i,"HR"]<-summary(f)$coefficient[1,2]
        h[i,"lci"]<-exp(confint(f)[1,1])
        h[i,"uci"]<-exp(confint(f)[1,2])
        h[i,"Cstat"]<-concordance(f)$concordance 
        
        fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~",seq_ids_en[i],"+",Model[m],"+female:",seq_ids_en[i],sep=""))
        fx<-coxph(fmla,data=d, na.action=na.exclude)
        lr=lrtest(f, fx)
        h[i,"p_lrtest_female"]<-lr[2,"Pr(>Chisq)"] #Performs a LRT to compare nested models, used here to evaluate the significance of adding interaction terms
        
        fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~",seq_ids_en[i],"+",Model[m],"+white:",seq_ids_en[i],sep=""))
        fx<-coxph(fmla,data=d, na.action=na.exclude)
        lr=lrtest(f, fx)
        h[i,"p_lrtest_white"]<-lr[2,"Pr(>Chisq)"]
      }
      h$p_FDR<-p.adjust(h$p, method = "fdr")
      h<-merge(annot_soma[,c("seqid_in_sample","uniprot_id","entrezgeneid","entrezgenesymbol","targetfullname")],h,by.x="seqid_in_sample",by.y="seq_id")
      h<-h[order(h$p),]
      write.table(h,file=paste(out.dir,"cox/cox_",outcome[k],"_Model2",model[m],dss[j],".txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
    }
  }
}





#using 85 from FDR
length(seq_ids_en2)
# N=85
Model<-c("+age+female+racegrp+center+elevel02+cigt+bmi+egfrcr+tcal+sprt+ethanl+diabetes+hypertens")

for(k in 1:length(outcome)) {  
  for(m in 1:length(Model)) {
    h<-data.frame(seq_id=character(),beta=double(),se=double(),p=double(),HR=double(),
                  lci=double(),uci=double(),stringsAsFactors=FALSE)
    
    for(i in 1:length(seq_ids_en2)) {
      h[i,"seq_id"]<-seq_ids_en2[i]
      
      fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~",seq_ids_en2[i],Model[m],sep=""))
      f<-coxph(fmla,data=d, na.action=na.exclude)
      h[i,"beta"]<-summary(f)$coefficient[1,1]
      h[i,"se"]<-summary(f)$coefficient[1,3]
      h[i,"p"]<-summary(f)$coefficient[1,5]
      h[i,"HR"]<-summary(f)$coefficient[1,2]
      h[i,"lci"]<-exp(confint(f)[1,1])
      h[i,"uci"]<-exp(confint(f)[1,2])
    }
    
    h$p_FDR<-p.adjust(h$p, method = "fdr")
    h<-merge(annot_soma[,c("seqid_in_sample","uniprot_id","entrezgeneid","entrezgenesymbol","targetfullname")],h,by.x="seqid_in_sample",by.y="seq_id")
    write.table(h,file=paste(out.dir,"cox/Cox_",outcome[k],"_Model2_fdr",model[m],".txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
  }
}

length(which(h$p<0.05))
# 45


## interaction with sex, race, diabetes and hypertens
h<-read.delim(paste0(out.dir, "cox/Cox_ckd_def1_v3_Model2_fdr1.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

d=soma_v3_SMP
h=h[!is.na(h$p),]
seq_ids_en3=h[h$p_FDR<0.05,"seqid_in_sample"]
length(seq_ids_en3)
#N=10

ds<-c("soma_v3_SMP")
dss<-c("")
outcome<-c("ckd_def1_v3")
event=c("inc_ckd_def1_v3")
fu=c("fu_ckd_def1_v3")

model<-c(1)

Model<-c("+age+female+white+center+elevel02+cigt+bmi+egfrcr+tcal+sprt+ethanl+diabetes+hypertens") #for white vs. nonwhite, replace racegrp with white in model

for(j in 1:length(ds)) {
  d=get(ds[j])
  print(dim(d))
  
  for(k in 1:length(outcome)) {
    for(m in 1:length(Model)) {
        
        h<-data.frame(seq_id=character(),
                      p_inter=double(),stringsAsFactors=FALSE)
        
        for(i in 1:length(seq_ids_en3)) {
          h[i,"seq_id"]<-seq_ids_en3[i]
          fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~",seq_ids_en3[i],Model[m],sep=""))
          f0<-coxph(fmla,data=d, na.action=na.exclude)
          fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~",seq_ids_en3[i],Model[m],"+",seq_ids_en3[i],"*white",sep=""))
          f1<-coxph(fmla,data=d, na.action=na.exclude)
          lr=lrtest(f1,f0)
          h[i,"p_inter"]<-lr[2,"Pr(>Chisq)"]   
          
      }
      h<-merge(annot_soma[,c("seqid_in_sample","uniprot_id","entrezgeneid","entrezgenesymbol","targetfullname")],h,by.x="seqid_in_sample",by.y="seq_id")
      write.table(h,file=paste(out.dir,"cox/cox (interaction with white)_",outcome[k],"_Model",model[m],dss[j],".txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
   
  }
}
}

      
Model<-c("+age+female+racegrp+center+elevel02+cigt+bmi+egfrcr+tcal+sprt+ethanl+diabetes+hypertens")
      
for(j in 1:length(ds)) {
  d=get(ds[j])
  print(dim(d))
  
  for(k in 1:length(outcome)) {
    for(m in 1:length(Model)) {
      
      # interaction, female
      h<-data.frame(seq_id=character(),
                    p_inter=double(),stringsAsFactors=FALSE)
      
      for(i in 1:length(seq_ids_en3)) {
        h[i,"seq_id"]<-seq_ids_en3[i]
        fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~",seq_ids_en3[i],Model[m],sep=""))
        f0<-coxph(fmla,data=d, na.action=na.exclude)
        fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~",seq_ids_en3[i],Model[m],"+",seq_ids_en3[i],"*female",sep=""))
        f1<-coxph(fmla,data=d, na.action=na.exclude)
        lr=lrtest(f1,f0)
        h[i,"p_inter"]<-lr[2,"Pr(>Chisq)"]
      }
      h<-merge(annot_soma[,c("seqid_in_sample","uniprot_id","entrezgeneid","entrezgenesymbol","targetfullname")],h,by.x="seqid_in_sample",by.y="seq_id")
      write.table(h,file=paste(out.dir,"cox/cox (interaction with female)_",outcome[k],"_Model",model[m],dss[j],".txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)

      
      # interaction, diabetes
      h<-data.frame(seq_id=character(),
                    p_inter=double(),stringsAsFactors=FALSE)
      
      for(i in 1:length(seq_ids_en3)) {
        h[i,"seq_id"]<-seq_ids_en3[i]
        fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~",seq_ids_en3[i],Model[m],sep=""))
        f0<-coxph(fmla,data=d, na.action=na.exclude)
        fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~",seq_ids_en3[i],Model[m],"+",seq_ids_en3[i],"*diabetes",sep=""))
        f1<-coxph(fmla,data=d, na.action=na.exclude)
        lr=lrtest(f1,f0)
        h[i,"p_inter"]<-lr[2,"Pr(>Chisq)"]
      }
      h<-merge(annot_soma[,c("seqid_in_sample","uniprot_id","entrezgeneid","entrezgenesymbol","targetfullname")],h,by.x="seqid_in_sample",by.y="seq_id")
      write.table(h,file=paste(out.dir,"cox/cox (interaction with diabetes)_",outcome[k],"_Model",model[m],dss[j],".txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
      
 
      # interaction, hypertens
      h<-data.frame(seq_id=character(),
                    p_inter=double(),stringsAsFactors=FALSE)
      
      for(i in 1:length(seq_ids_en3)) {
        h[i,"seq_id"]<-seq_ids_en3[i]
        fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~",seq_ids_en3[i],Model[m],sep=""))
        f0<-coxph(fmla,data=d, na.action=na.exclude)
        fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~",seq_ids_en3[i],Model[m],"+",seq_ids_en3[i],"*hypertens",sep=""))
        f1<-coxph(fmla,data=d, na.action=na.exclude)
        lr=lrtest(f1,f0)
        h[i,"p_inter"]<-lr[2,"Pr(>Chisq)"]
      }
      h<-merge(annot_soma[,c("seqid_in_sample","uniprot_id","entrezgeneid","entrezgenesymbol","targetfullname")],h,by.x="seqid_in_sample",by.y="seq_id")
      write.table(h,file=paste(out.dir,"cox/cox (interaction with hypertens)_",outcome[k],"_Model",model[m],dss[j],".txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)

    }
  }
}
    


###############output subgroup hr and ci
h<-read.delim(paste0(out.dir, "cox/Cox_ckd_def1_v3_Model2_fdr1.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

d=soma_v3_SMP
h=h[!is.na(h$p),]
seq_ids_en3=h[h$p_FDR<0.05,"seqid_in_sample"]
length(seq_ids_en3)
#N=10

#####for male - output hr and ci
Model<-c("+age+racegrp+center+elevel02+cigt+bmi+egfrcr+tcal+sprt+ethanl+diabetes+hypertens")

for(k in 1:length(outcome)) {  
  for(m in 1:length(Model)) {
    h<-data.frame(seq_id=character(),beta=double(),se=double(),p=double(),HR=double(),
                  lci=double(),uci=double(),stringsAsFactors=FALSE)
    
    for(i in 1:length(seq_ids_en3)) {
      h[i,"seq_id"]<-seq_ids_en3[i]
      
      fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~",seq_ids_en3[i],Model[m],sep=""))
      f<-coxph(fmla,data=d_male, na.action=na.exclude)
      h[i,"beta"]<-summary(f)$coefficient[1,1]
      h[i,"se"]<-summary(f)$coefficient[1,3]
      h[i,"p"]<-summary(f)$coefficient[1,5]
      h[i,"HR"]<-summary(f)$coefficient[1,2]
      h[i,"lci"]<-exp(confint(f)[1,1])
      h[i,"uci"]<-exp(confint(f)[1,2])
    }
    h<-merge(annot_soma[,c("seqid_in_sample","uniprot_id","entrezgeneid","entrezgenesymbol","targetfullname")],h,by.x="seqid_in_sample",by.y="seq_id")
    h$p_FDR<-p.adjust(h$p, method = "fdr")
    write.table(h,file=paste(out.dir,"cox/Cox_",outcome[k],"_Model2_male",model[m],".txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
  }
}


#####for female - output hr and ci
Model<-c("+age+racegrp+center+elevel02+cigt+bmi+egfrcr+tcal+sprt+ethanl+diabetes+hypertens")

for(k in 1:length(outcome)) {  
  for(m in 1:length(Model)) {
    h<-data.frame(seq_id=character(),beta=double(),se=double(),p=double(),HR=double(),
                  lci=double(),uci=double(),stringsAsFactors=FALSE)
    
    for(i in 1:length(seq_ids_en3)) {
      h[i,"seq_id"]<-seq_ids_en3[i]
      
      fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~",seq_ids_en3[i],Model[m],sep=""))
      f<-coxph(fmla,data=d_female, na.action=na.exclude)
      h[i,"beta"]<-summary(f)$coefficient[1,1]
      h[i,"se"]<-summary(f)$coefficient[1,3]
      h[i,"p"]<-summary(f)$coefficient[1,5]
      h[i,"HR"]<-summary(f)$coefficient[1,2]
      h[i,"lci"]<-exp(confint(f)[1,1])
      h[i,"uci"]<-exp(confint(f)[1,2])
    }
    h<-merge(annot_soma[,c("seqid_in_sample","uniprot_id","entrezgeneid","entrezgenesymbol","targetfullname")],h,by.x="seqid_in_sample",by.y="seq_id")
    h$p_FDR<-p.adjust(h$p, method = "fdr")
    write.table(h,file=paste(out.dir,"cox/Cox_",outcome[k],"_Model2_female",model[m],".txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
  }
}


#####for white - output hr and ci
Model<-c("+age+female+center+elevel02+cigt+bmi+egfrcr+tcal+sprt+ethanl+diabetes+hypertens")

for(k in 1:length(outcome)) {  
  for(m in 1:length(Model)) {
    h<-data.frame(seq_id=character(),beta=double(),se=double(),p=double(),HR=double(),
                  lci=double(),uci=double(),stringsAsFactors=FALSE)
    
    for(i in 1:length(seq_ids_en3)) {
      h[i,"seq_id"]<-seq_ids_en3[i]
      
      fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~",seq_ids_en3[i],Model[m],sep=""))
      f<-coxph(fmla,data=d_white, na.action=na.exclude)
      h[i,"beta"]<-summary(f)$coefficient[1,1]
      h[i,"se"]<-summary(f)$coefficient[1,3]
      h[i,"p"]<-summary(f)$coefficient[1,5]
      h[i,"HR"]<-summary(f)$coefficient[1,2]
      h[i,"lci"]<-exp(confint(f)[1,1])
      h[i,"uci"]<-exp(confint(f)[1,2])
    }
    h<-merge(annot_soma[,c("seqid_in_sample","uniprot_id","entrezgeneid","entrezgenesymbol","targetfullname")],h,by.x="seqid_in_sample",by.y="seq_id")
    h$p_FDR<-p.adjust(h$p, method = "fdr")
    write.table(h,file=paste(out.dir,"cox/Cox_",outcome[k],"_Model2_white",model[m],".txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
  }
}


#####for nonwhite - output hr and ci
Model<-c("+age+female+center+elevel02+cigt+bmi+egfrcr+tcal+sprt+ethanl+diabetes+hypertens")

for(k in 1:length(outcome)) {  
  for(m in 1:length(Model)) {
    h<-data.frame(seq_id=character(),beta=double(),se=double(),p=double(),HR=double(),
                  lci=double(),uci=double(),stringsAsFactors=FALSE)
    
    for(i in 1:length(seq_ids_en3)) {
      h[i,"seq_id"]<-seq_ids_en3[i]
      
      fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~",seq_ids_en3[i],Model[m],sep=""))
      f<-coxph(fmla,data=d_nonwhite, na.action=na.exclude)
      h[i,"beta"]<-summary(f)$coefficient[1,1]
      h[i,"se"]<-summary(f)$coefficient[1,3]
      h[i,"p"]<-summary(f)$coefficient[1,5]
      h[i,"HR"]<-summary(f)$coefficient[1,2]
      h[i,"lci"]<-exp(confint(f)[1,1])
      h[i,"uci"]<-exp(confint(f)[1,2])
    }
    h<-merge(annot_soma[,c("seqid_in_sample","uniprot_id","entrezgeneid","entrezgenesymbol","targetfullname")],h,by.x="seqid_in_sample",by.y="seq_id")
    h$p_FDR<-p.adjust(h$p, method = "fdr")
    write.table(h,file=paste(out.dir,"cox/Cox_",outcome[k],"_Model2_nonwhite",model[m],".txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
  }
}



#####for hypertension - output hr and ci
Model<-c("+age+female+racegrp+center+elevel02+cigt+bmi+egfrcr+tcal+sprt+ethanl+diabetes")

for(k in 1:length(outcome)) {  
  for(m in 1:length(Model)) {
    h<-data.frame(seq_id=character(),beta=double(),se=double(),p=double(),HR=double(),
                  lci=double(),uci=double(),stringsAsFactors=FALSE)
    
    for(i in 1:length(seq_ids_en3)) {
      h[i,"seq_id"]<-seq_ids_en3[i]
      
      fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~",seq_ids_en3[i],Model[m],sep=""))
      f<-coxph(fmla,data=d_hypertens, na.action=na.exclude)
      h[i,"beta"]<-summary(f)$coefficient[1,1]
      h[i,"se"]<-summary(f)$coefficient[1,3]
      h[i,"p"]<-summary(f)$coefficient[1,5]
      h[i,"HR"]<-summary(f)$coefficient[1,2]
      h[i,"lci"]<-exp(confint(f)[1,1])
      h[i,"uci"]<-exp(confint(f)[1,2])
    }
    h<-merge(annot_soma[,c("seqid_in_sample","uniprot_id","entrezgeneid","entrezgenesymbol","targetfullname")],h,by.x="seqid_in_sample",by.y="seq_id")
    h$p_FDR<-p.adjust(h$p, method = "fdr")
    write.table(h,file=paste(out.dir,"cox/Cox_",outcome[k],"_Model2_hypertens",model[m],".txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
  }
}


#####for nohypertens - output hr and ci
Model<-c("+age+female+racegrp+center+elevel02+cigt+bmi+egfrcr+tcal+sprt+ethanl+diabetes")

for(k in 1:length(outcome)) {  
  for(m in 1:length(Model)) {
    h<-data.frame(seq_id=character(),beta=double(),se=double(),p=double(),HR=double(),
                  lci=double(),uci=double(),stringsAsFactors=FALSE)
    
    for(i in 1:length(seq_ids_en3)) {
      h[i,"seq_id"]<-seq_ids_en3[i]
      
      fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~",seq_ids_en3[i],Model[m],sep=""))
      f<-coxph(fmla,data=d_nohypertens, na.action=na.exclude)
      h[i,"beta"]<-summary(f)$coefficient[1,1]
      h[i,"se"]<-summary(f)$coefficient[1,3]
      h[i,"p"]<-summary(f)$coefficient[1,5]
      h[i,"HR"]<-summary(f)$coefficient[1,2]
      h[i,"lci"]<-exp(confint(f)[1,1])
      h[i,"uci"]<-exp(confint(f)[1,2])
    }
    h<-merge(annot_soma[,c("seqid_in_sample","uniprot_id","entrezgeneid","entrezgenesymbol","targetfullname")],h,by.x="seqid_in_sample",by.y="seq_id")
    h$p_FDR<-p.adjust(h$p, method = "fdr")
    write.table(h,file=paste(out.dir,"cox/Cox_",outcome[k],"_Model2_nohypertens",model[m],".txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
  }
}

###sensitivity###
###for follow-up at least 5 year - examine all 85 from validation - seq_ids_en2
Model<-c("+age+female+racegrp+center+elevel02+cigt+bmi+egfrcr+tcal+sprt+ethanl+diabetes+hypertens")

for(k in 1:length(outcome)) {  
  for(m in 1:length(Model)) {
    h<-data.frame(seq_id=character(),beta=double(),se=double(),p=double(),HR=double(),
                  lci=double(),uci=double(),stringsAsFactors=FALSE)
    
    for(i in 1:length(seq_ids_en2)) {
      h[i,"seq_id"]<-seq_ids_en2[i]
      
      fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~",seq_ids_en2[i],Model[m],sep=""))
      f<-coxph(fmla,data=d_5yr_less, na.action=na.exclude)
      h[i,"beta"]<-summary(f)$coefficient[1,1]
      h[i,"se"]<-summary(f)$coefficient[1,3]
      h[i,"p"]<-summary(f)$coefficient[1,5]
      h[i,"HR"]<-summary(f)$coefficient[1,2]
      h[i,"lci"]<-exp(confint(f)[1,1])
      h[i,"uci"]<-exp(confint(f)[1,2])
    }
    h<-merge(annot_soma[,c("seqid_in_sample","uniprot_id","entrezgeneid","entrezgenesymbol","targetfullname")],h,by.x="seqid_in_sample",by.y="seq_id")
    h$p_FDR<-p.adjust(h$p, method = "fdr")
    write.table(h,file=paste(out.dir,"cox/Cox_",outcome[k],"_Model2_5yr_less",model[m],".txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
  }
}




#######################################################################################
### 6) plots for cox
###############
## heatmap
h<-read.delim(paste0(out.dir, "cox/Cox_ckd_def1_v3_Model2_fdr1.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

d=soma_v3_SMP
h=h[!is.na(h$p),]
seq_ids_en3=h[h$p_FDR<0.05,"seqid_in_sample"]
length(seq_ids_en3)
#N=10
cormat=cor(d[,seq_ids_en3],method="spearman")

get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
reorder_cormat <- function(cormat){
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}
# Reorder the correlation matrix
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
lower_tri <- get_lower_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
melted_cormat$i=1:nrow(melted_cormat)
melted_cormat=merge(melted_cormat,annot_soma[,c("seqid_in_sample","entrezgenesymbol")],by.x="Var1",by.y="seqid_in_sample")
colnames(melted_cormat)[colnames(melted_cormat)=="entrezgenesymbol"] <- "entrezgenesymbol1"
melted_cormat=merge(melted_cormat,annot_soma[,c("seqid_in_sample","entrezgenesymbol")],by.x="Var2",by.y="seqid_in_sample")
colnames(melted_cormat)[colnames(melted_cormat)=="entrezgenesymbol"] <- "entrezgenesymbol2"
melted_cormat=melted_cormat[order(melted_cormat$i),]

# Create a ggheatmap
out.plot=paste(out.dir,"plots/heatmaps/","Heatmap_cox_density_model2_fdr.tiff",sep="")
tiff(filename = out.plot,width=4800, height=4800, res=400)

ggheatmap <- ggplot(melted_cormat, aes(x=reorder(entrezgenesymbol2,i), y=reorder(entrezgenesymbol1,i), fill = value))+
  geom_tile()+
  geom_text(aes(x=reorder(entrezgenesymbol2,i), y=reorder(entrezgenesymbol1,i), label = sprintf("%.2f",value)), color = "black", size = 5.5) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="") +
  theme_minimal()+ # minimal theme
  coord_fixed() + 
  theme(axis.text.x=element_text(color = "black", size=16, angle=90, vjust=.8, hjust=0.8)) + 
  theme(axis.text.y=element_text(color = "black", size=16, angle=0, vjust=.8, hjust=0.8)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = c(0.2, 0.8)) +
  ylab("") +
  xlab("") +
  scale_y_discrete(position = "right") 
# Print the heatmap
print(ggheatmap)

dev.off()





####################################
#######cox-pota intake and ckd######
d=soma_v3_SMP

cox_continous <- coxph(Surv(fu_ckd_def1_v3, inc_ckd_def1_v3) ~ density+age+female+racegrp+center+elevel02+cigt+bmi+egfrcr+tcal+sprt+ethanl+diabetes+hypertens, data = d)
summary(cox_continous)






#######################################################################################################################
## 7) mediation analysis
ds=c("soma_v3_SMP")
dss=c("")
model<-c(1)
Model=c("+age+female+racegrp+center+elevel02+cigt+bmi+egfrcr+tcal+sprt+ethanl+diabetes+hypertens")
outcome<-c("ckd_def1_v3")
event=c("inc_ckd_def1_v3")
fu=c("fu_ckd_def1_v3")  


for(j in 1:length(ds)) {
  d=get(paste0(ds[j]))
  annot=annot_soma
  
  for(k in 1:length(outcome)) {
    for(m in 1:length(model)) {
      t=read.delim(paste(out.dir,"cox/Cox_ckd_def1_v3_Model2_fdr1.txt",sep=""), header = TRUE, sep = "\t", stringsAsFactors =FALSE)
      t=t[!is.na(t$p),]
      t=t[t$p_FDR<0.05,]
      met_sel=t$seqid_in_sample
      print(length(met_sel))
      
      h<-data.frame(mid=character(),
                    ACME=double(),p_ACME=double(),ADE=double(),p_ADE=double(),
                    total_effect=double(),p_total_effect=double(),prop_mediated=double(),p_prop_mediated=double(),
                    stringsAsFactors=FALSE)
      
      for(i in 1:length(met_sel)) {
        h[i,"mid"]=met_sel[i]
        
        fmla<-as.formula(paste(event[k],"~density+",met_sel[i],Model[m],sep=""))
        out.fit<-glm(fmla,data=d,family= binomial("probit"))
        #print(summary(out.fit))
        
        fmla<-as.formula(paste(met_sel[i],"~density",Model[m],sep=""))
        med.fit<-lm(fmla,data=d, na.action=na.exclude)
        #print(summary(med.fit))
        
        results = mediate(med.fit, out.fit, treat='density', mediator=met_sel[i], 
                          robustSE = TRUE)  
        #print(summary(results))
        h[i,"ACME"]=results$d.avg
        h[i,"p_ACME"]=results$d.avg.p
        h[i,"ADE"]=results$z.avg
        h[i,"p_ADE"]=results$z.avg.p
        h[i,"prop_mediated"]=results$n.avg
        h[i,"p_prop_mediated"]=results$n.avg.p
        h[i,"total_effect"]=results$tau.coef
        h[i,"p_total_effect"]=results$tau.p
      }
      
      h<-merge(annot[,c("seqid_in_sample","uniprot_id","entrezgeneid","entrezgenesymbol","targetfullname")],h,by.x="seqid_in_sample",by.y="mid")
      write.table(h,file=paste(out.dir,"mediation analysis_",outcome[k],"_model",model[m],"_density",dss[j],".txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
    }
  }
}


#CODE FOR hr VS beta plot
#######################################################################################
#######################################################################################
### 6.1) plotting HR (for prospective analysis of proteins + CKD) vs. beta coefficients (for cross-sectional analysis of potassium + proteins)
h1<-read.delim(paste0(out.dir,"linear/Linear_density_Model1_validation.txt"), header = TRUE, sep = "\t", stringsAsFactors =FALSE)
#h2=read.delim(paste(out.dir,"cox/Cox_ckd_def1_v3_Model2_fdr1.txt",sep=""), header = TRUE, sep = "\t", stringsAsFactors =FALSE)

h2=h2[h2$p_FDR<0.05,]
dim(h2)
# 10

for(k in 1:length(outcome)) {
  h2=read.delim(paste(out.dir,"cox/Cox_",outcome[k],"_Model2_fdr1.txt",sep=""), header = TRUE, sep = "\t", stringsAsFactors =FALSE)
  h2=h2[h2$p_FDR<0.05,]
  print(nrow(h2))
  h=merge(h1[,c("seqid_in_sample","beta","entrezgenesymbol")],h2[,c("seqid_in_sample","HR")],by="seqid_in_sample")
  
  h[(h$beta<0 & h$HR>1) | (h$beta>0 & h$HR<1),"c"]=1
  h[is.na(h$c),"c"]=2
  h$c=as.factor(h$c)
  cols <- c("1" = "red", "2" = "black")
  
  out.plot=paste(out.dir,"plots/HR vs beta_",outcome[k],".tiff",sep="")

tiff(filename = out.plot,width = 800, height = 800)
print(
  ggplot(h, aes(beta, HR, label = entrezgenesymbol,color=c)) +
    scale_color_manual(values = cols) +
    guides(colour=FALSE) +
    geom_text_repel(size = 5,max.overlaps=30) +
    geom_point(size = 3) +
    theme_classic(base_size = 16) +
    xlab("Beta") +
    ylab("HR")  +
    geom_vline(xintercept=0,col="black",lty="dotted") +
    geom_hline(yintercept=1,col="black",lty="dotted") +
    coord_trans(y = "log2")
)
dev.off()
}

########################################


########################################################################################################################
########################################################################################################################
########################################################################################################################
### 6) elastic net, cox models
### Elastic net with cox regression
elasticNetCox <- function(candiateVars, forcedVars = forced_in_Vars, data = d.temp,
                          endpoint, time){
  
  if(length(candiateVars)>0 & length(c(candiateVars,forcedVars))>=2) {
    ## Cross validation for hyperparamter tuning
    cv.fit <- cv.glmnet(as.matrix(data[, c(forcedVars, candiateVars)]),
                        Surv(time, endpoint), family="cox", alpha = 0.5, maxit = 2000, nfolds=5,
                        penalty.factor=c(rep(0, length(forcedVars)), rep(1, length(candiateVars))))
    ## Fit the model
    fit <- glmnet(as.matrix(data[, c(forcedVars, candiateVars)]),
                  Surv(time, endpoint), family="cox", alpha = 0.5, maxit = 2000,
                  penalty.factor=c(rep(0, length(forcedVars)), rep(1, length(candiateVars))))
    #browser()
    Coefficients <- coef(fit, s = cv.fit$lambda.1se)
    nonzeorcoefs <- Coefficients[which(Coefficients != 0)]
    nonzeor_protein = names(nonzeorcoefs) = rownames(Coefficients)[ which(Coefficients != 0)]
    #nonzeor_flags = sapply(nonzeor_protein, function(x) flag_status[which(names(flag_status) == x)] )
    #nonzeor_names = sapply(nonzeor_protein, function(x) names(proteList_all)[which(proteList_all == x)] )
    SurvObj<-paste("Surv(",fu[k],",",event[k],")~",sep="")
    if(length(nonzeorcoefs)>0) {
      cox_formula <- as.formula(paste(SurvObj, paste(c(nonzeor_protein, forcedVars), collapse= "+") ))
      final_cox <- coxph(cox_formula, data = data, method = "efron")
      # Base Model
      hr_fM = summary(final_cox)$coefficients[1:length(nonzeor_protein), 1]
      se_fM = summary(final_cox)$coefficients[1:length(nonzeor_protein), 3]
      p_fM = summary(final_cox)$coefficients[1:length(nonzeor_protein), 5]
      hr_fMs = paste0(round(hr_fM, 2), " (", round(se_fM, 2), ")" )
      p_fMs = sapply(formatC(p_fM, format = "e", digits = 1), toString)
      
      elsNet_results = cbind(unlist(nonzeor_protein), 
                             round(exp(nonzeorcoefs), 2), round(nonzeorcoefs, 2),
                             hr_fMs, p_fMs) %>% data.table()
      names(elsNet_results) <-c("seq_id", "HR(elasticNet)", "logHR(elasticNet)", "logHR(Cox)", "Pvalue(Cox)")
      return(list(elsNet_results = elsNet_results, optimalNet = Coefficients, finalCox = final_cox,fit=fit,cv.fit=cv.fit,flagx=0))
    } else {
      #cox_formula <- as.formula(paste(SurvObj, "1"))
      return(list(elsNet_results = NA,  finalCox = NA,flagx=1))
    }
  } else if(length(forcedVars)!=0) {
    SurvObj<-paste("Surv(",fu[k],",",event[k],")~",sep="")
    cox_formula <- as.formula(paste(SurvObj, paste(c(forcedVars), collapse= "+") ))
    final_cox <- coxph(cox_formula, data = data, method = "efron")
    # Base Model
    hr_fM = summary(final_cox)$coefficients[1:length(forcedVars), 1]
    se_fM = summary(final_cox)$coefficients[1:length(forcedVars), 3]
    p_fM = summary(final_cox)$coefficients[1:length(forcedVars), 5]
    hr_fMs = paste0(round(hr_fM, 2), " (", round(se_fM, 2), ")" )
    p_fMs = sapply(formatC(p_fM, format = "e", digits = 1), toString)
    
    a<-NA
    b<-NA
    elsNet_results = cbind(unlist(forcedVars),a,b,hr_fMs, p_fMs) %>% data.table()
    names(elsNet_results) <-c("seq_id", "HR(elasticNet)", "logHR(elasticNet)", "logHR(Cox)", "Pvalue(Cox)")
    return(list(elsNet_results = elsNet_results, finalCox = final_cox, flagx=0))
  } else {
    return(list(elsNet_results = NA,  finalCox = NA,flagx=1))
  }
}


get.risk.coxph.ex <-
  function(mdl, t0, lp) {
    bash    = basehaz(mdl)
    lambda0 = approx(bash$time, bash$hazard, t0)$y
    risk    = 1 - exp( - lambda0 * exp( lp ) )
    return(risk)
  }





###################################################################################
# cross-validation for cox
ds=c("soma_v3_SMP")
dss=c("")
outcome<-c("ckd_def2_v3")
event=c("inc_ckd_def2_v3")
fu=c("fu_ckd_def2_v3")  
E_covar1=c("age","female","racegrp_B","racegrp_A","racegrp_I",
           "centerF","centerJ","centerW",
           "elevel02_2","elevel02_3","cigt_current","cigt_former",
           "bmi","egfrcr","tcal","sprt","ethanl","diabetes","hypertens")
EN<-c("EN1","EN1C","EN1CF","M1","M1F")



for(j in 1:length(ds)) {
  d<-get(ds[j])
  annot=annot_soma
  
  for(k in 1:length(outcome)) {
    d.temp=d[!is.na(d[,event[k]]),]
    d.temp=d[d[,fu[k]]>0,]
    print(dim(d.temp))
    
    c_validation<-data.frame(model=character(),c_predall=double(),se_predall=double(),
                             sum_c1=double(),sum_c2=double(),
                             sum_c3=double(),sum_c4=double(),sum_c5=double(),sum_c6=double(),
                             sum_c7=double(),sum_c8=double(),sum_c9=double(),sum_c10=double(),
                             stringsAsFactors=FALSE)
    
    
    index<-d.temp[,event[k]]==1
    d.temp.e=d.temp[index,]
    d.temp.ne=d.temp[!index,]
    
    set.seed(2524)
    #set.seed(12345)
    folds.e <- split(sample(nrow(d.temp.e), nrow(d.temp.e),replace=FALSE), as.factor(1:10))
    set.seed(9527)
    #set.seed(23333)
    folds.ne <- split(sample(nrow(d.temp.ne), nrow(d.temp.ne),replace=FALSE), as.factor(1:10))
    
    d.temp2<-d.temp
    
    for(n in 1:length(EN)) {
      
      new.pred.all<-NULL
      cv_var_coef<-NULL
      prisk.all=NULL
      
      if (grepl("EN",EN[n])) {
        num=str_extract(EN[n], "\\-*\\d+\\.*\\d*")
        t<-read.delim(file =paste(out.dir,"cox/Cox_",outcome[k],"_Model2_fdr1.txt",sep=""), header = TRUE, sep = "\t", stringsAsFactors =FALSE)
        t<-t[t$p_FDR<0.05 & !is.na(t$p),]
        t=t[order(t$p),]
        print(nrow(t))
        tops=t$seqid_in_sample
        vars=tops
      } else {
        vars=NULL
      }
      
      if (grepl("C",EN[n]) | grepl("M",EN[n])) {
        num=str_extract(EN[n], "\\-*\\d+\\.*\\d*")
        C=get(paste("E_covar",num,sep=""))
      } else {
        C=NULL
      }
      
      if (grepl("F",EN[n])) {
        forced_in_Vars=C
      } else {
        forced_in_Vars=NULL
      }
      
      if (!grepl("F",EN[n]) & (grepl("C",EN[n]) | grepl("M",EN[n]))) {
        vars=c(vars,C)
      }
      
      
      for(q in 1:10) {
        
        d1<-rbind.data.frame(d.temp.e[folds.e[[q]],],d.temp.ne[folds.ne[[q]],])
        d2<-rbind.data.frame(d.temp.e[-folds.e[[q]],],d.temp.ne[-folds.ne[[q]],])
        
        print(paste(EN[n],"__",q,sep=""))
        
        elsNet_results_part =  elasticNetCox(candiateVars=vars, forcedVars = forced_in_Vars, data = d2,
                                             endpoint = d2[,event[k]], time = d2[,fu[k]])
        
        if(elsNet_results_part$flagx==0) {
          
          print("spot1")
          
          f<-elsNet_results_part$finalCox
          
          print("spot2")
          
          new.pred = predict(f, d1)
          new.pred1=data.frame(new.pred)
          rownames(new.pred1)=d1$SampleId
          
          prisk=get.risk.coxph.ex(f, t0=5,lp=new.pred)
          prisk1=data.frame(prisk)
          rownames(prisk1)=d1$SampleId
          
          print("spot3")
          
          new.pred.all<-rbind.data.frame(new.pred.all,new.pred1)
          prisk.all<-rbind.data.frame(prisk.all,prisk1)
          print(dim(d1))
          print(length(new.pred))
          print(dim(new.pred1))
          print(dim(new.pred.all))
          print(length(prisk))
          print(dim(prisk1))
          print(dim(prisk.all))
          
          print("spot4")
          
          fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~ new.pred",sep=""))
          
          print("spot5")
          
          res.1 = survConcordance(fmla, d1)
          
          print("spot6")
          
          c_validation[n,"model"]=EN[n]
          c_validation[n,paste("sum_c",q,sep="")]=res.1$concordance
          
          
          # save coefficients
          coef<-data.frame(coef(elsNet_results_part$finalCox))
          colnames(coef)<-paste("coef_loop",q,sep="")
          coef$Var<-row.names(coef)
          if (q==1 | is.null(cv_var_coef)) {
            cv_var_coef=coef
          } else if(q>1) {
            cv_var_coef<-merge(cv_var_coef,coef,by="Var",all=TRUE)
          }
          
        } else {
          c_validation[n,"model"]=EN[n]
          if(q==1) {
            cv_var_coef=NULL
          }
        }
        
      }
      
      
      if(!is.null(new.pred.all)) {
        print(dim(d.temp2))
        names(new.pred.all)=paste0("new.pred_",EN[n])
        new.pred.all$SampleId=rownames(new.pred.all)
        d.temp2<-merge(d.temp2,new.pred.all,by="SampleId",all=T)
        print(dim(d.temp2))
        
        names(prisk.all)=paste0("prisk_",EN[n])
        prisk.all$SampleId=rownames(prisk.all)
        d.temp2<-merge(d.temp2,prisk.all,by="SampleId",all=T)
        print(dim(d.temp2))
        
        fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~ new.pred_",EN[n],sep=""))
        res.2 = survConcordance(fmla, d.temp2)
        c_validation[n,"c_predall"]=res.2$concordance
        c_validation[n,"se_predall"]=res.2$std.err
      }
      
      if(!is.null(cv_var_coef)) {
        coefs<-grep("coef_loop", names(cv_var_coef), value=TRUE)
        cv_var_coef$N_select=rowSums(!is.na(cv_var_coef[,coefs]))
        cv_var_coef<-merge(annot[,c("seqid_in_sample","uniprot_id","entrezgeneid","entrezgenesymbol","target","targetfullname")],cv_var_coef,by.x="seqid_in_sample",by.y="Var",all.y=TRUE)
        write.table(cv_var_coef,file=paste(out.dir,"elastic net/cox/10foldCV/Elastic net_cox_coefficients_",EN[n],"_",outcome[k],dss[j],".txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
      }
    }
    
    sumCs<-grep("sum_c", names(c_validation), value=TRUE) 
    c_validation$mean_c=rowMeans(c_validation[,sumCs])
    #c_validation$se_c=rowSds(as.matrix(c_validation[,sumCs]))/sqrt(10)
    write.table(c_validation,file=paste(out.dir,"elastic net/cox/10foldCV/Elastic net_cox_Cstat_",outcome[k],dss[j],".txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
  }
}





### counting how many times selected
for(j in 1:length(ds)) {
  for(k in 1:length(outcome)) {
    
    cv_var_coef_all<-NULL
    
    for(n in 1:length(EN)) {
      cv_var_coef=read.delim(file = paste(out.dir,"elastic net/cox/10foldCV/Elastic net_cox_coefficients_",EN[n],"_",outcome[k],dss[j],".txt",sep=""), header = TRUE, sep = "\t", stringsAsFactors =FALSE)
      cv_var_coef<-cv_var_coef[order(-cv_var_coef$N_select),]
      cv_var_coef$Model<-EN[n]
      cv_var_coef<-cv_var_coef[cv_var_coef$N_select>=7,c("Model","seqid_in_sample","uniprot_id","entrezgeneid","entrezgenesymbol","target","targetfullname","N_select")]
      
      if(n==1) {
        cv_var_coef_all<-cv_var_coef
      } else {
        cv_var_coef_all<-rbind.data.frame(cv_var_coef_all,cv_var_coef)
      }
    }
    
    write.table(cv_var_coef_all,file=paste(out.dir,"elastic net/cox/10foldCV/Elastic net_variables selected 7 times or more_", outcome[k],dss[j],".txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
  }
}



#######################################################################################################################
## 6) c stat comparison using 10 fold CV results
# cox model
ds=c("soma_v3_SMP")
dss=c("")
model<-c(1)
Model=c("age+female+racegrp+center+elevel02+cigt+bmi+egfrcr+tcal+sprt+ethanl+diabetes+hypertens")
outcome<-c("ckd_def2_v3")
event=c("inc_ckd_def2_v3")
fu=c("fu_ckd_def2_v3")  


for(j in 1:length(ds)) {
  d=get(paste0(ds[j]))
  annot=annot_soma
  t=read.delim(file=paste0(out.dir,"elastic net/cox/10foldCV/Elastic net_variables selected 7 times or more_ckd_def2_v3",dss[j],".txt"), header = TRUE, sep = "\t", stringsAsFactors =FALSE)
  t=t[t$Model=="EN1C",]
  t<-t[t$N_select>=7,]
  met_sel=t$seqid_in_sample
  met_sel=grep("SeqId_",met_sel,value=T)
  print(length(met_sel))
  
  for(k in 1:length(outcome)) {
    for(m in 1:length(model)) {
      h<-data.frame(mid=character(),
                    Cstat=double(),Cstat_lci=double(),Cstat_uci=double(),
                    Cstat_noprotein=double(),Cstat_noprotein_lci=double(),Cstat_noprotein_uci=double(),
                    Cstat_diff=double(),Cstat_diff_lci=double(),Cstat_diff_uci=double(),
                    p_Cstat_diff=double(),
                    stringsAsFactors=FALSE)
      
      # covariates only
      fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~",Model[m],sep=""))
      f<-coxph(fmla,data=d, na.action=na.exclude)
      c0<-concordance(f)$concordance
      c0_lci=concordance(f)$concordance-1.96*sqrt(concordance(f)$var)
      c0_uci=concordance(f)$concordance+1.96*sqrt(concordance(f)$var)     
      d[,paste0("new.pred0")] = predict(f, d)
      
      # single protein
      for(i in 1:length(met_sel)) {
        h[i,"mid"]<-met_sel[i]
        fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~",met_sel[i],"+",Model[m],sep=""))
        f<-coxph(fmla,data=d, na.action=na.exclude)
        d[,paste0("new.pred1")] = predict(f, d)
        h[i,"Cstat"]<-concordance(f)$concordance
        h[i,"Cstat_lci"]<-concordance(f)$concordance-1.96*sqrt(concordance(f)$var)
        h[i,"Cstat_uci"]<-concordance(f)$concordance+1.96*sqrt(concordance(f)$var)
        h[i,"Cstat_noprotein"]<-c0
        h[i,"Cstat_noprotein_lci"]<-c0_lci
        h[i,"Cstat_noprotein_uci"]<-c0_uci
        compare<-compareC(d[,fu[k]], d[,event[k]], d[,paste0("new.pred0")],  d[,paste0("new.pred1")])
        h[i,"Cstat_diff"]<-compare$est.diff_c
        h[i,"Cstat_diff_lci"]<-compare$est.diff_c-1.96*sqrt(compare$est.vardiff_c)
        h[i,"Cstat_diff_uci"]<-compare$est.diff_c+1.96*sqrt(compare$est.vardiff_c)
        h[i,"p_Cstat_diff"]<-compare$pval                         
      }
      
      # all proteins together
      i=i+1
      h[i,"mid"]<-"all proteins together"
      fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~",paste(met_sel,collapse="+"),"+",Model[m],sep=""))
      f<-coxph(fmla,data=d, na.action=na.exclude)
      d[,paste0("new.pred1")] = predict(f, d)
      h[i,"Cstat"]<-concordance(f)$concordance
      h[i,"Cstat_lci"]<-concordance(f)$concordance-1.96*sqrt(concordance(f)$var)
      h[i,"Cstat_uci"]<-concordance(f)$concordance+1.96*sqrt(concordance(f)$var)
      h[i,"Cstat_noprotein"]<-c0
      h[i,"Cstat_noprotein_lci"]<-c0_lci
      h[i,"Cstat_noprotein_uci"]<-c0_uci
      compare<-compareC(d[,fu[k]], d[,event[k]], d[,paste0("new.pred0")],  d[,paste0("new.pred1")])
      h[i,"Cstat_diff"]<-compare$est.diff_c
      h[i,"Cstat_diff_lci"]<-compare$est.diff_c-1.96*sqrt(compare$est.vardiff_c)
      h[i,"Cstat_diff_uci"]<-compare$est.diff_c+1.96*sqrt(compare$est.vardiff_c)
      h[i,"p_Cstat_diff"]<-compare$pval  
      
      h<-merge(annot[,c("seqid_in_sample","uniprot_id","entrezgeneid","entrezgenesymbol","targetfullname")],h,by.x="seqid_in_sample",by.y="mid",all.y=T)
      write.table(h,file=paste(out.dir,"cox/Cstats to predict ",outcome[k],"_model",model[m],dss[j],".txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
    }
  }
}



##C-statistics using 10 cox hits
# cox model
ds=c("soma_v3_SMP")
dss=c("")
model<-c(1)
Model=c("age+female+racegrp+center+elevel02+cigt+bmi+egfrcr+tcal+sprt+ethanl+diabetes+hypertens")
outcome<-c("ckd_def1_v3")
event=c("inc_ckd_def1_v3")
fu=c("fu_ckd_def1_v3")  


for(j in 1:length(ds)) {
  d=get(paste0(ds[j]))
  annot=annot_soma
  t<-read.delim(file =paste(out.dir,"cox/Cox_",outcome[k],"_Model2_fdr1.txt",sep=""), header = TRUE, sep = "\t", stringsAsFactors =FALSE)
  t<-t[t$p_FDR<0.05 & !is.na(t$p),]
  met_sel=t$seqid_in_sample
  met_sel=grep("SeqId_",met_sel,value=T)
  print(length(met_sel))
  
  for(k in 1:length(outcome)) {
    for(m in 1:length(model)) {
      h<-data.frame(mid=character(),
                    Cstat=double(),Cstat_lci=double(),Cstat_uci=double(),
                    Cstat_noprotein=double(),Cstat_noprotein_lci=double(),Cstat_noprotein_uci=double(),
                    Cstat_diff=double(),Cstat_diff_lci=double(),Cstat_diff_uci=double(),
                    p_Cstat_diff=double(),
                    stringsAsFactors=FALSE)
      
      # covariates only
      fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~",Model[m],sep=""))
      f<-coxph(fmla,data=d, na.action=na.exclude)
      c0<-concordance(f)$concordance
      c0_lci=concordance(f)$concordance-1.96*sqrt(concordance(f)$var)
      c0_uci=concordance(f)$concordance+1.96*sqrt(concordance(f)$var)     
      d[,paste0("new.pred0")] = predict(f, d)
      
      # single protein
      for(i in 1:length(met_sel)) {
        h[i,"mid"]<-met_sel[i]
        fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~",met_sel[i],"+",Model[m],sep=""))
        f<-coxph(fmla,data=d, na.action=na.exclude)
        d[,paste0("new.pred1")] = predict(f, d)
        h[i,"Cstat"]<-concordance(f)$concordance
        h[i,"Cstat_lci"]<-concordance(f)$concordance-1.96*sqrt(concordance(f)$var)
        h[i,"Cstat_uci"]<-concordance(f)$concordance+1.96*sqrt(concordance(f)$var)
        h[i,"Cstat_noprotein"]<-c0
        h[i,"Cstat_noprotein_lci"]<-c0_lci
        h[i,"Cstat_noprotein_uci"]<-c0_uci
        compare<-compareC(d[,fu[k]], d[,event[k]], d[,paste0("new.pred0")],  d[,paste0("new.pred1")])
        h[i,"Cstat_diff"]<-compare$est.diff_c
        h[i,"Cstat_diff_lci"]<-compare$est.diff_c-1.96*sqrt(compare$est.vardiff_c)
        h[i,"Cstat_diff_uci"]<-compare$est.diff_c+1.96*sqrt(compare$est.vardiff_c)
        h[i,"p_Cstat_diff"]<-compare$pval                         
      }
      
      # all proteins together
      i=i+1
      h[i,"mid"]<-"all proteins together"
      fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~",paste(met_sel,collapse="+"),"+",Model[m],sep=""))
      f<-coxph(fmla,data=d, na.action=na.exclude)
      d[,paste0("new.pred1")] = predict(f, d)
      h[i,"Cstat"]<-concordance(f)$concordance
      h[i,"Cstat_lci"]<-concordance(f)$concordance-1.96*sqrt(concordance(f)$var)
      h[i,"Cstat_uci"]<-concordance(f)$concordance+1.96*sqrt(concordance(f)$var)
      h[i,"Cstat_noprotein"]<-c0
      h[i,"Cstat_noprotein_lci"]<-c0_lci
      h[i,"Cstat_noprotein_uci"]<-c0_uci
      compare<-compareC(d[,fu[k]], d[,event[k]], d[,paste0("new.pred0")],  d[,paste0("new.pred1")])
      h[i,"Cstat_diff"]<-compare$est.diff_c
      h[i,"Cstat_diff_lci"]<-compare$est.diff_c-1.96*sqrt(compare$est.vardiff_c)
      h[i,"Cstat_diff_uci"]<-compare$est.diff_c+1.96*sqrt(compare$est.vardiff_c)
      h[i,"p_Cstat_diff"]<-compare$pval  
      
      h<-merge(annot[,c("seqid_in_sample","uniprot_id","entrezgeneid","entrezgenesymbol","targetfullname")],h,by.x="seqid_in_sample",by.y="mid",all.y=T)
      write.table(h,file=paste(out.dir,"cox/Cstats to predict(10) ",outcome[k],"_model",model[m],dss[j],".txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
    }
  }
}







#######################Calculate protein signature risk score########################
#####################################################################################
#30 proteins from elastic net - extract their beta from linear model to calculate risk score
outcome<-c("ckd_def1_v3")
model<-c(1)
event=c("inc_ckd_def1_v3")
fu=c("fu_ckd_def1_v3")
Model<-c("+age+female+racegrp+center+elevel02+cigt+bmi+egfrcr+tcal+sprt+ethanl+diabetes+hypertens")

d=soma_v3_SMP

g<-read.delim(paste0(out.dir,"elastic net/logistic/10foldCV/Elastic net_logistic_coef_EN1CF_cat_density.txt"), header = TRUE, sep = "\t", stringsAsFactors =FALSE)
g=g[g$N_select>=7,]
g=g[!(g$seqid_in_sample %in% c("(Intercept)")),]
elastic_net_sel <- g$seqid_in_sample


linear_results <- read.delim(paste0(out.dir, "linear/Linear_density_Model1_validation2.txt"),
                             header = TRUE, sep = "\t", stringsAsFactors = FALSE)
filtered_linear_results <- linear_results[linear_results$seqid_in_sample %in% elastic_net_sel, ]
selected_variables <- filtered_linear_results$seqid_in_sample  # Variables
selected_betas <- filtered_linear_results$beta                 # Beta coefficients


d$score <- 0
for (i in seq_along(selected_variables)) {
  d$score <- d$score + d[[selected_variables[i]]] * selected_betas[i]
}

summary(d$score)

# standardize risk score
d[,"score_orig"]=d[,"score"]
sd=sd(d[,"score"],na.rm=T)
d[,"score"]=d[,"score"]/sd

summary(d$score)

#Cox Proportional Hazards Models
for (k in 1:length(outcome)) {
  for (m in 1:length(Model)) {

        h <- data.frame(beta = double(), se = double(), p = double(),
                    HR = double(), lci = double(), uci = double(), stringsAsFactors = FALSE)
    
    i <- 1
    fmla <- as.formula(paste("Surv(", fu[k], ",", event[k], ")~score", Model[m], sep = ""))
    
    f <- coxph(fmla, data = d, na.action = na.exclude)
    
    h[i, "beta"] <- summary(f)$coefficient[1, 1]   
    h[i, "se"] <- summary(f)$coefficient[1, 3]    
    h[i, "p"] <- summary(f)$coefficient[1, 5]     
    h[i, "HR"] <- summary(f)$coefficient[1, 2]    
    h[i, "lci"] <- exp(confint(f)[1, 1])          
    h[i, "uci"] <- exp(confint(f)[1, 2])          
    
    write.table(h, file = paste(out.dir, "cox/Cox (score)_", outcome[k], "_30_", model[m], "_2.txt", sep = ""),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}




#######
#1 proteins from cox elastic net - extract their beta from linear model to calculate risk score#################
outcome<-c("ckd_def1_v3")
model<-c(1)
event=c("inc_ckd_def1_v3")
fu=c("fu_ckd_def1_v3")
Model<-c("+age+female+racegrp+center+elevel02+cigt+bmi+egfrcr+tcal+sprt+ethanl+diabetes+hypertens")

d=soma_v3_SMP

g<-read.delim(paste0(out.dir,"elastic net/cox/10foldCV/Elastic net_cox_coefficients_EN1C_ckd_def1_v3.txt"), header = TRUE, sep = "\t", stringsAsFactors =FALSE)
g=g[g$N_select>=7,]
g=g[!(g$seqid_in_sample %in% c("(Intercept)")),]
elastic_net_sel <- g$seqid_in_sample

linear_results <- read.delim(paste0(out.dir, "linear/Linear_density_Model1_validation2.txt"),
                             header = TRUE, sep = "\t", stringsAsFactors = FALSE)

filtered_linear_results <- linear_results[linear_results$seqid_in_sample %in% elastic_net_sel, ]
selected_variables <- filtered_linear_results$seqid_in_sample  # Variables
selected_betas <- filtered_linear_results$beta                 # Beta coefficients


d$score <- 0
for (i in seq_along(selected_variables)) {
  d$score <- d$score + d[[selected_variables[i]]] * selected_betas[i]
}

summary(d$score)

# standardize risk score
d[,"score_orig"]=d[,"score"]
sd=sd(d[,"score"],na.rm=T)
d[,"score"]=d[,"score"]/sd

summary(d$score)

#Cox Proportional Hazards Models
for (k in 1:length(outcome)) {
  for (m in 1:length(Model)) {
    h <- data.frame(beta = double(), se = double(), p = double(),
                    HR = double(), lci = double(), uci = double(), stringsAsFactors = FALSE)
    
    i <- 1
    fmla <- as.formula(paste("Surv(", fu[k], ",", event[k], ")~score", Model[m], sep = ""))
    
    f <- coxph(fmla, data = d, na.action = na.exclude)
    
    h[i, "beta"] <- summary(f)$coefficient[1, 1]   
    h[i, "se"] <- summary(f)$coefficient[1, 3]    
    h[i, "p"] <- summary(f)$coefficient[1, 5]     
    h[i, "HR"] <- summary(f)$coefficient[1, 2]    
    h[i, "lci"] <- exp(confint(f)[1, 1])          
    h[i, "uci"] <- exp(confint(f)[1, 2])          
    
    write.table(h, file = paste(out.dir, "cox/Cox (score)_", outcome[k], "_1_", model[m], "_2.txt", sep = ""),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}



#10 proteins from cox - extract their beta from linear model to calculate risk score#################
outcome<-c("ckd_def1_v3")
model<-c(1)
event=c("inc_ckd_def1_v3")
fu=c("fu_ckd_def1_v3")
Model<-c("+age+female+racegrp+center+elevel02+cigt+bmi+egfrcr+tcal+sprt+ethanl+diabetes+hypertens")

d=soma_v3_SMP

g<-read.delim(file =paste(out.dir,"cox/Cox_",outcome[k],"_Model2_fdr1.txt",sep=""), header = TRUE, sep = "\t", stringsAsFactors =FALSE)
g<-g[g$p_FDR<0.05 & !is.na(g$p),]
g=g[!(g$seqid_in_sample %in% c("(Intercept)")),]
elastic_net_sel <- g$seqid_in_sample

linear_results <- read.delim(paste0(out.dir, "linear/Linear_density_Model1_validation2.txt"),
                             header = TRUE, sep = "\t", stringsAsFactors = FALSE)

filtered_linear_results <- linear_results[linear_results$seqid_in_sample %in% elastic_net_sel, ]
selected_variables <- filtered_linear_results$seqid_in_sample  # Variables
selected_betas <- filtered_linear_results$beta                 # Beta coefficients


d$score <- 0
for (i in seq_along(selected_variables)) {
  d$score <- d$score + d[[selected_variables[i]]] * selected_betas[i]
}

summary(d$score)

# standardize risk score
d[,"score_orig"]=d[,"score"]
sd=sd(d[,"score"],na.rm=T)
d[,"score"]=d[,"score"]/sd

summary(d$score)


#Cox Proportional Hazards Models
for (k in 1:length(outcome)) {
  for (m in 1:length(Model)) {
    
    h <- data.frame(beta = double(), se = double(), p = double(),
                    HR = double(), lci = double(), uci = double(), stringsAsFactors = FALSE)
    
    i <- 1
    fmla <- as.formula(paste("Surv(", fu[k], ",", event[k], ")~score", Model[m], sep = ""))
    
    f <- coxph(fmla, data = d, na.action = na.exclude)
    
    h[i, "beta"] <- summary(f)$coefficient[1, 1]   
    h[i, "se"] <- summary(f)$coefficient[1, 3]    
    h[i, "p"] <- summary(f)$coefficient[1, 5]     
    h[i, "HR"] <- summary(f)$coefficient[1, 2]    
    h[i, "lci"] <- exp(confint(f)[1, 1])          
    h[i, "uci"] <- exp(confint(f)[1, 2])          
    
    write.table(h, file = paste(out.dir, "cox/Cox (score)_", outcome[k], "_10cox_", model[m], "_2.txt", sep = ""),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}



#6 proteins from cox mediators - extract their beta from linear model to calculate risk score#################
outcome<-c("ckd_def1_v3")
model<-c(1)
event=c("inc_ckd_def1_v3")
fu=c("fu_ckd_def1_v3")
Model<-c("+age+female+racegrp+center+elevel02+cigt+bmi+egfrcr+tcal+sprt+ethanl+diabetes+hypertens")

d=soma_v3_SMP

g<-read.delim(file =paste(out.dir,"mediation analysis_ckd_def1_v3_model1_density.txt",sep=""), header = TRUE, sep = "\t", stringsAsFactors =FALSE)
g<-g[g$p_prop_mediated<=0.05,]
g=g[!(g$seqid_in_sample %in% c("(Intercept)")),]
elastic_net_sel <- g$seqid_in_sample

linear_results <- read.delim(paste0(out.dir, "linear/Linear_density_Model1_validation2.txt"),
                             header = TRUE, sep = "\t", stringsAsFactors = FALSE)

filtered_linear_results <- linear_results[linear_results$seqid_in_sample %in% elastic_net_sel, ]
selected_variables <- filtered_linear_results$seqid_in_sample  # Variables
selected_betas <- filtered_linear_results$beta                 # Beta coefficients


d$score <- 0
for (i in seq_along(selected_variables)) {
  d$score <- d$score + d[[selected_variables[i]]] * selected_betas[i]
}

summary(d$score)

# standardize risk score
d[,"score_orig"]=d[,"score"]
sd=sd(d[,"score"],na.rm=T)
d[,"score"]=d[,"score"]/sd

summary(d$score)


#Cox Proportional Hazards Models
for (k in 1:length(outcome)) {
  for (m in 1:length(Model)) {

        h <- data.frame(beta = double(), se = double(), p = double(),
                    HR = double(), lci = double(), uci = double(), stringsAsFactors = FALSE)
    
    i <- 1
    fmla <- as.formula(paste("Surv(", fu[k], ",", event[k], ")~score", Model[m], sep = ""))
    
    f <- coxph(fmla, data = d, na.action = na.exclude)
    
    h[i, "beta"] <- summary(f)$coefficient[1, 1]   
    h[i, "se"] <- summary(f)$coefficient[1, 3]    
    h[i, "p"] <- summary(f)$coefficient[1, 5]     
    h[i, "HR"] <- summary(f)$coefficient[1, 2]    
    h[i, "lci"] <- exp(confint(f)[1, 1])          
    h[i, "uci"] <- exp(confint(f)[1, 2])          
    
    write.table(h, file = paste(out.dir, "cox/Cox (score)_", outcome[k], "_6mediator_", model[m], "_2.txt", sep = ""),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

