#################################
# Load external functions.
#
# library(seqinr)
# source("detection/detection-V3.R")
# source("detection/getPtt.R")
#################################
# Generic function: from a sequence and a annotation objects, 
# gives a # list with general info and a sublist with genes 
# and attributes.
# fna is from SeqFastadna class, ptt obtained by getPtt() method.
#
getGenes<-function(fna, ptt, stat=FALSE)
{
tit<-titre(fna[[1]])
nwin<-length(ptt$gene)
seq<-1:nwin
countseq<-array(dim=c(nwin,4))
strand<-gsub(" ?\\+","1",ptt$strand)
strand<-gsub(" ?\\-","-1",strand)
strand<-as.numeric(strand)
for(i in 1:nwin){
  sta<-ptt$begin[i]
  end<-ptt$end[i]
  len<- end-sta+1
  if(len <= 0){
    seq[i] <- NA
    countseq[i,] <- NA
  } else {
    frag<-getFrag(getSequence(fna[[1]], as.string=FALSE),sta,end)
    seq[i]<-c2s(getSequence(frag, as.string=FALSE))
    if (strand[i] == -1) {
      seq[i]<-c2s(rev(comp(s2c(seq[i]))))
      }
    for(j in 1:4) countseq[i,j]<-count(s2c(seq[i]),1)[[j]]
  }
}
###
# Basic statistics on genes
A<-countseq[,1]
C<-countseq[,2]
G<-countseq[,3]
T<-countseq[,4]
tot<-A+C+G+T
lis<-list(title=tit,start=ptt$begin,stop=ptt$end,strand=strand,length=ptt$length,pid=ptt$pid,gene=ptt$gene,synonym=ptt$synonym,code=ptt$code,cog=ptt$cog,product=ptt$product,seq=seq,A=A,C=C,G=G,T=T)
if (stat==TRUE) {
  GC<-(G+C)/tot
  PL<-(G+A)/tot
  GCsk<-(G-C)/(G+C)
  TAsk<-(T-A)/(T+A)
  lis<-c(lis,list(GC=GC,PL=PL,GCsk=GCsk,TAsk=TAsk))
  }
lis
}
#################################
# RSCU calculation from reference set of genes.
# df is obtained by the getGenes() method. Vector of relevant set of
# reference is given either by a research with grep method or by
# the vector of indices vect. By default, Bulmer's (1988) correction
# is made
#
RSCUcal<-function(df,crit,vect=FALSE,Bulmer=T)
{
if(vect==FALSE) {
  vect<-grep(crit, df$product)
  }
nrefset<-length(vect)
if (nrefset<1) stop("Number in ref set too low")
s<-""
for(i in 1:nrefset) {
  s<-paste(s,df$seq[vect[i]],sep="")
  }
uscod<-uco(s2c(s),as.data.frame=T)
rn<-as.character(uscod$AA)
ncodperaa<-1:64
mx<-1:64
w<-1:64
for (i in 1:64){
  ncodperaa[i]<-length(grep(substr(rn[i],1,3),rn))
  sub<-grep(substr(rn[i],1,3),rn)
  mx[i]<-max(uscod$eff[sub])
  w[i]<-uscod$eff[i]/mx[i]
  if (w[i]<0.01 && Bulmer) w[i]=0.01
  }
cbind(uscod,ncodperaa,mx,w)
}
#################################
# CAI calculation for a dataframe of genes obtained with getGenes()
# and a RSCU set obtained with RSCU(), following Sharp and Li (1987)
#
CAIcal<-function(df, RSCU) 
{
nwin<-length(df$seq)
CAI<-1:nwin
for (i in 1:nwin){
  uscod<-uco(s2c(df$seq[i]),as.data.frame=T)
  x<-uscod$eff*log(RSCU$w)
  strip<-grep("taa|tga|tag|atg|tgg",as.character(uscod$codon))
  x<-sum(x[-strip])
  CAI[i]<-exp(x/nchar(df$seq[i])*3)
  }
CAI
}
#################################
# Codon usage differences between gene classes (Karlin et al, 1998).
# By default, the ref set (vect2 or crit2) is the whole set in the df.
#
Bcal<-function(df,crit1="each",crit2="all",vect1=FALSE,vect2=FALSE)
{
nwin<-length(df$seq)
###
# Choice of groups
if(vect2==FALSE) {
  if (crit2=="all"){ 
    vect2<-1:nwin
    }
  else{
    vect2<-grep(crit2, c(df$product,df$gene))
    vect2<-vect2%%nwin
    }
  }
###
# Definition of reference group
uscod2<-uco(s2c(df$seq[1]),as.data.frame=T)
uscod2$eff<-rep(0,64)
rn<-as.character(uscod2$AA)
#rn<-row.names(uscod2)
for (i in vect2){
  uscod2$eff<-uscod2$eff+uco(s2c(df$seq[i]),as.data.frame=T)$eff
  }
sm2<-1:64
for (i in 1:64){
  sub<-grep(substr(rn[i],1,3),rn)
  sm2[i]<-sum(uscod2$eff[sub])
  }
freq2<-uscod2$eff/sm2
uscod2<-cbind(uscod2,syn_freq=freq2)
sm1<-1:64
###
# Definition of query group: calculation for all genes
if(crit1=="each"){
  bias<-1:nwin
  for (i in 1:nwin){
    uscod1<-uco(s2c(df$seq[i]),as.data.frame=T)
    for (j in 1:64){
      sub<-grep(substr(rn[j],1,3),rn)
      sm1[j]<-sum(uscod1$eff[sub])
      }
    aaf<-sm1/sum(uscod1$eff)
    freq1<-uscod1$eff/sm1
    uscod1<-cbind(uscod1,syn_freq=freq1)
    biases<-aaf*abs(uscod1$syn_freq-uscod2$syn_freq)
    bias[i]<-sum(biases,na.rm=T)
    }
  bias
  }
###
# Definition of query group: calculation for a group
else{
  if(vect1==FALSE) {
    vect1<-grep(crit1, c(df$product,df$gene))
    vect1<-vect1%%nwin
    }
  uscod1<-uco(s2c(df$seq[1]),as.data.frame=T)
  uscod1$eff<-rep(0,64)
  for (i in vect1){
    uscod1$eff<-uscod1$eff+uco(s2c(df$seq[i]),as.data.frame=T)$eff
    }
  for (i in 1:64){
    sub<-grep(substr(rn[i],1,3),rn)
    sm1[i]<-sum(uscod1$eff[sub])
    }
  aaf<-sm1/sum(uscod1$eff)
  freq1<-uscod1$eff/sm1
  biases<-aaf*abs(uscod1$eff-escod2$eff)
  bias<-sum(biases,na.rm=T)
  bias
  }
}
#################################
# Acal: Calculate the aminoacid bias between each gene of a genome and 
# the average of the genome. In future, should be able to adapt to make
# comparisons between two groups, whatever they are.
#
Acal <- function(df){
  nwin <- length(df$seq)
  stat <- summary.SeqFastaAA(as.SeqFastaAA(translate(s2c(df$seq[1]))))$comp
  names <- names(stat)
  mat <- matrix(nrow=nwin, ncol=21)
  mat <- as.data.frame(mat)
  names(mat)<-names
  for(i in 1:nwin){
     stat <- summary.SeqFastaAA(as.SeqFastaAA(translate(s2c(df$seq[i]))))$comp
     mat[i,] <- t(as.matrix(stat))
  }
  mat<-mat[,-1]
  average <- apply(mat,2,mean)
  bias <- numeric(nwin)
  for(i in 1:nwin){
    bias[i] <- sum(abs(mat[i,]-average))/20
  }
  bias
}
#################################
# Calculates P, or G+C content of any codon position.
#
getP<-function(df,pos=3,strip=F)
{
uscod<-uco(s2c(df$seq[1]),as.data.frame=T)
uscod$eff<-rep(0,64)
cod<-uscod$codon
s<-switch(EXPR=pos,"c..|g..",".c.|.g.","..c|..g")
if (strip){
  stripsub<-grep("taa|tga|tag|atg|tgg|ata",as.character(uscod$codon))
  cod[stripsub]<-NA
  }
sub2<-grep(s,as.character(cod))
nwin<-length(df$seq)
gc<-1:nwin
for (i in 1:nwin){
  uscod<-uco(s2c(df$seq[i]),as.data.frame=TRUE)
  if(strip) {
    gc[i]<-sum(uscod$eff[sub2])/sum(uscod$eff[-stripsub])
    } else gc[i]<-sum(uscod$eff[sub2])/sum(uscod$eff)
  }
gc
}
#################################
# Calculates H-phobicity of a df of genes
#
getHphob<-function(df)
{
hphob<-c(-0.8,-0.8,-0.8,-0.8,3.8,3.8,3.8,3.8,4.2,4.2,-0.4,-0.4,0,2.5,4.5,4.5,-4.5,-4.5,-0.9,2.5,3.8,2.8,3.8,2.8,-4.5,-4.5,-3.5,-3.5,-3.5,-3.9,-3.5,-3.5,-3.2,-3.5,-3.2,-0.4,-0.4,-3.5,-3.9,-3.5,-1.6,-1.6,-1.6,-1.6,-0.7,-0.7,-0.7,-0.7,4.2,4.2,-4.5,-0.8,1.9,1.8,1.8,1.8,1.8,4.5,-4.5,-0.8,0,-1.3,0,-1.3)
nwin<-length(df$seq)
h<-1:nwin
for (i in 1:nwin){
  uscod<-uco(s2c(df$seq[i]),as.data.frame=TRUE)
  h[i]<-sum(uscod$eff*hphob)/sum(uscod$eff)
  }
h
}
#################################
# Calculates GC or TA skew for any codon position
#
getSkewPosCod<-function(df,skew="GC",pos=3)
{
uscod<-uco(s2c(df$seq[1]),as.data.frame=T)
uscod$eff<-rep(0,64)
cod<-uscod$codon
if (skew=="TA"){ 
  s1<-switch(EXPR=pos,"t..",".t.","..t")
  s2<-switch(EXPR=pos,"a..",".a.","..a")
  }
else{
  s1<-switch(EXPR=pos,"g..",".g.","..g")
  s2<-switch(EXPR=pos,"c..",".c.","..c")
  }
sub1<-grep(s1,as.character(cod))
sub2<-grep(s2,as.character(cod))
nwin<-length(df$seq)
sk<-1:nwin
for (i in 1:nwin){
  uscod<-uco(s2c(df$seq[i]),as.data.frame=TRUE)
  sum1<-sum(uscod$eff[sub1])
  sum2<-sum(uscod$eff[sub2])
  sk[i]<-(sum1-sum2)/(sum1+sum2)
  }
sk
}
#################################
# Calculates GC or TA skew for any codon position
#
countNtPosCod<-function(df,pos=3){
  uscod<-uco(s2c(df$seq[1]),as.data.frame=T)
  uscod$eff<-rep(0,64)
  cod<-uscod$codon
  s1<-switch(EXPR=pos,"a..",".a.","..a")
  s2<-switch(EXPR=pos,"c..",".c.","..c")
  s3<-switch(EXPR=pos,"g..",".g.","..g")
  s4<-switch(EXPR=pos,"t..",".t.","..t")
  sub1<-grep(s1,as.character(cod))
  sub2<-grep(s2,as.character(cod))
  sub3<-grep(s3,as.character(cod))
  sub4<-grep(s4,as.character(cod))
  nwin<-length(df$seq)
  sk<-matrix(NA, nrow=nwin, ncol=4)
  for (i in 1:nwin){
    uscod<-uco(s2c(df$seq[i]),as.data.frame=TRUE)
    sk[i,1]<-sum(uscod$eff[sub1])
    sk[i,2]<-sum(uscod$eff[sub2])
    sk[i,3]<-sum(uscod$eff[sub3])
    sk[i,4]<-sum(uscod$eff[sub4])
  }
  sk <- as.data.frame(sk)
  names(sk) <- c("A","C","G","T")
  sk
}
#################################
# Function to count nucleotides in array of sequences
#
countNt<-function(seq,t="t",stat=TRUE,n=1,both=FALSE,row.names=NULL){
  nr<-length(seq)
  #res<-array(dim=c(nr,4^n))
  res<-NULL
  for (i in 1:nr){
    s<-as.character(seq[i])
    s<-tolower(s)
    if (t=="u") s<-gsub("u","t",s)
    if (both) s<-paste(s, c2s(rev(comp(s2c(s)))),sep="x")
    s<-s2c(s)
    c<-count(s,n)
    res<-cbind(res,c)
  }
  res<-t(res)
  row.names(res)<-row.names
  res
}
#################################
# Calculates the dinucleotide bias B, the bias between 
# the observed and the theoretical frequency of diNt.
# Input can be sequences (default) or matrices of nt/dint.
# The Campbell distance is calculated.
#
countB<-function(seq1,seq2,mat1=F,mat2=F,t="t",both=TRUE,row.names=NULL)
{
if (mat1==F){
  seq2<-countNt(seq1,t,n=2,both=both,row.names)
  seq1<-countNt(seq1,t,n=1,both=both,row.names)
  }
n<-length(seq1[,1])
tot1<-apply(seq1,1,sum)
tot2<-apply(seq2,1,sum)
fseq1<-seq1/tot1
fo<-seq2/tot2
ft<-fo
for (i in 1:n){
  ft[i,]<-as.vector(outer(fseq1[i,],fseq1[i,]))
  }
b<-(fo-ft)/ft
row.names(b)<-row.names
list(b=b,dist=dist(b, method="manhattan",diag=T,upper=T)*1000/16)
}
################################
# Split a fasta sequence into an array of windows.
# Lwin must be a multiple of step.
# Takes care of the last window.
#
splitFasta<-function(fasta,lwin=1000,step=lwin){
  if(lwin %% step != 0) stop("lwin must be a multiple of step")
  len <- length(fasta)
  nwin <- ((len-lwin-1) %/% step) + 2
	seq <- splitseq(fasta,word=step)
  if(lwin != step){
    nstep_per_win <- lwin %/% step
    seq2 <- character(nwin)
    for(i in 1:(nwin-1)){
      seq2[i] <- paste(seq[i:(i+nstep_per_win-1)],sep="",collapse="")
    }
    seq <- seq2
  }
	frag <- getFrag(getSequence(fasta, as.string=F),(nwin-1)*step+1,len)
  seq[nwin] <- c2s(getSequence(frag, as.string=F))
  seq
}
################################
# Calculates the skews and cumulative skews in an array of sequences
# 
skews<-function(seq,t="t",norm=TRUE){
	nr<-length(seq)
	res<-array(dim=c(nr,4))
	res<-as.data.frame(res)
	names(res)<-c("gcsk","task","cgcsk","ctask")
	cgcsk<-0
	ctask<-0
	cnt<-countNt(seq,t=t,stat=TRUE,n=1,both=FALSE,row.names=NULL)
	for (i in 1:nr){
		res[i,1] <- cnt[i,3]-cnt[i,2]
    if(norm) {res[i,1] <- res[i,1]/(cnt[i,3]+cnt[i,2])}
		cgcsk <- cgcsk+res[i,1]
		res[i,3] <- cgcsk
		res[i,2] <- (cnt[i,4]-cnt[i,1])
    if(norm) {res[i,2] <- res[i,2] /(cnt[i,4]+cnt[i,1])}
		ctask <- ctask+res[i,2]
		res[i,4] <- ctask
	}
	res
}

