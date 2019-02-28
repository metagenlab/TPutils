####################################
# Recherche du titre dans la sequence
#
titre<-function(x)
{ 
 b<-strsplit(getAnnot(x)[[1]],"\\|")
 titre<-b[[1]][length(b[[1]])]
 titre
}
####################################
# Transformation du ptt si le 1er gene a une adresse de  
# debut superieure a son adresse de fin
ajustptt<-function(ptt,l)
{ 
 ptt$begin<-c(ptt$begin,ptt$begin[1])
 ptt$end<-c(ptt$end,ptt$end[1])
 ptt$strand<-c(ptt$strand,ptt$strand[1])
 ptt$length<-c(ptt$length,ptt$length[1])
 ptt$pid<-c(ptt$pid,ptt$pid[1])
 ptt$gene<-c(ptt$gene,ptt$gene[1])
 ptt$synonym<-c(ptt$synonym,ptt$synonym[1])
 ptt$cod<-c(ptt$cod,ptt$cod[1])
 ptt$cog<-c(ptt$cog,ptt$cog[1])
 ptt$product<-c(ptt$product,ptt$product[1])
 ptt$begin[1]<-1
 ptt$end[length(ptt$end)]<-l
 ajustptt<-ptt
}
####################################
# Calcul de la moyenne des pentes et des ecarts-types sur 1 tableau 
# 
pentesansterm<-function(resid,ctr,inter) 
{
lresid<-length(resid)
lvalpente<-(lresid-inter)
valpente<-(1:lvalpente)
for(i in (inter+1):lresid) {
 valpente[i-inter]<-(resid[i]-resid[i-inter])/(ctr[i]-ctr[i-inter])
}
moy<-mean(valpente)
ecart<-sd(valpente)
#plot(ctr[(1+(inter/2)):(lresid-(inter/2))],valpente)
lim<-list(moy=moy,ecart=ecart)
pentesansterm<-lim
}
####################################
# Calcul de la moyenne des pentes et des ecarts-types sur 2 tableaux 
# 
pente<-function(residot,residto,ctrot,ctrto,inter) 
{
lresidot<-length(residot)
lresidto<-length(residto)
lvalpente<-(lresidot+lresidto-(2*inter))
valpente<-(1:lvalpente)
for(i in (inter+1):lresidot) {
 valpente[i-inter]<-(residot[i]-residot[i-inter])/(ctrot[i]-ctrot[i-inter])
}
for(i in (inter+1):lresidto) {
 valpente[lresidot+i-(2*inter)]<-(residto[i]-residto[i-inter])/(ctrto[i]-ctrto[i-inter])
}
moy<-mean(valpente)
ecart<-sd(valpente)
lim<-list(moy=moy,ecart=ecart)
pente<-lim
}
####################################
# Detection de pentes
#
cherche<-function(resid,ctr,inter,coef,seuil) 
{
lresid<-length(resid)
positif<-(1:lresid)
j<-1;
haut<-seuil$moy+(coef*seuil$ecart)
bas<-seuil$moy-(coef*seuil$ecart)
for(i in (inter+1):lresid) {
 dif<-(resid[i]-resid[(i-inter)])/(ctr[i]-ctr[(i-inter)])
 if ((dif < bas) | (dif > haut))  {
  positif[j]<-(i-inter)
  j<-j+1
 }
}
positif[j]<-0
print(j)
#print(positif[1:j])
k<-0
if (positif[1]>0) {
 mat<-(1:(2*(j-1)))
 dim(mat)<-c((j-1),2)
 k<-1
 mat[1,1]<-positif[1]
 i<-2
 while (positif[i]>0) {
  if (positif[i]>(positif[i-1]+inter)) {
   mat[k,2]<-positif[i-1]+inter
   k<-k+1
   mat[k,1]<-positif[i]
  }
 i<-i+1
 }
 mat[k,2]<-positif[i-1]+inter
 resultats<-(1:(2*k))
 dim(resultats)<-c(k,2)
 print(k)
 for(i in 1:k) {
  resultats[i,1]<-mat[i,1]
  resultats[i,2]<-mat[i,2]
 }
}
else {
 resultats<-(1:2)
 dim(resultats)<-c(1,2)
 resultats[1,]<-1
}
cherche<-resultats
}
####################################
# Calcul des T-A (y=1) ou G-C (y=2)
#
difnuc<-function(x,y)
{ 
s<-summary.SeqFastadna(x)
difnuc<-s$composition[[(5-y)]]-s$composition[[y]]
}
####################################
# Determination des T-A ou G-C, des cumuls, et des centres de fenetres
#
calcul<-function(x,y,z)
{
long<-length(x)
nbwin<-long/y
ite<-trunc(nbwin)
juste<-ite==nbwin
if (juste==FALSE) nbwin<-(ite+1)
tab<-1:nbwin
cumul<-1:nbwin
centres<-1:nbwin
for(i in 1:ite) {
 deb<-y*(i-1)+1
 fin<-y*i
 frag<-getFrag(x,deb,fin)
 tab[i]<-difnuc(frag,z)
 centres[i]<-fin-(y/2)
 if (i==1) cumul[i]<-tab[i] else cumul[i]<-tab[i]+cumul[i-1]
}
if (juste==FALSE) {
 frag<-getFrag(x,(ite*y)+1,long)
 tab[nbwin]<-difnuc(frag,z)
 centres[nbwin]<-long-trunc((long-(ite*y)-1)/2)
 cumul[nbwin]<-tab[nbwin]+cumul[nbwin-1]
}
calcul<-list(tab=tab,cumul=cumul,centres=centres,nbwin=nbwin,long=long)
}
####################################
# Determination des T-A ou G-C, des cumuls, et des centres de g�nes
#
calgen<-function(x,ptt,z) {
nbwin<-length(ptt$gene)
tab<-1:nbwin
cumul<-1:nbwin
centres<-1:nbwin
centrepis<-1:nbwin
long<-0
for(i in 1:nbwin) {
 deb<-ptt$begin[i]
 fin<-ptt$end[i]
 lfrag<-fin-deb+1
###
# Sense analysis
 if (z==3) {
  tab[i]<-senscpt(lfrag,ptt$strand[i])
 }
###
# Sur tous les nucl�otides
 else if (z<3) {
  frag<-getFrag(x,deb,fin)
  tab[i]<-difnuc(frag,z)
 }
###
# Sur nucl�otides 1 � 3 
 else if (z<10) {
  tysk<-trunc((z-1)/3)
  posnuc<-(z-(3*(tysk-1))-3)
  if (ptt$strand[i]==" -"|ptt$strand[i]=="-") posnuc<-(4-posnuc)
  frag<-(1:(lfrag/3))
  for (j in 1:(lfrag/3)) frag[j]<-x[((deb-1)+(j-1)*3 +posnuc)]
  tab[i]<-difnuc(frag,tysk)
 }
 else {
  print("Erreur parametre sktype")
 }
 centres[i]<-deb+trunc((fin-deb)/2)
 centrepis[i]<-long+1+trunc((fin-deb)/2)
 long<-lfrag+long
 if (i==1) cumul[i]<-tab[i] else cumul[i]<-tab[i]+cumul[i-1]
}
calgen<-list(tab=tab,cumul=cumul,centres=centres,centrepis=centrepis,nbwin=nbwin,long=long)
}
####################################
# Determination de l'origine et du terminus
# termtype 1 : Origine au maximum du GC skew, Terminus au minimum
# termtype 2 : Origine au minimum du GC skew, Terminus au maximum
# termtype 3 : Origine au maximum du TA skew, Terminus au minimum
# termtype 4 : Origine au minimum du TA skew, Terminus au maximum
# x : sequence
# y : taille de fenetre de decoupage du genome
#
outerm<-function(x,y,termtype) 
{
if (termtype<3) rescal<-calcul(x,y,2) else rescal<-calcul(x,y,1)
if ((termtype==1)|(termtype==3)) {
 Terminus<-which.min(rescal$cumul)
 Origine<-which.max(rescal$cumul)
}
else {
 Terminus<-which.max(rescal$cumul)
 Origine<-which.min(rescal$cumul)
}
outerm<-list(Origine=Origine,Orcen=rescal$centres[Origine],Terminus=Terminus,Termcen=rescal$centres[Terminus])
}
####################################
# Determination de l'origine et du terminus
# termtype 5 : Origine au minimum du Sense Analysis, Terminus au maximum
# x   : sequence
# ptt : detail par gene issu du getPtt
#
outermgen<-function(x,ptt,termtype) 
{
rescal<-calgen(x,ptt,3)
Terminus<-which.max(rescal$cumul)
Origine<-which.min(rescal$cumul)
outerm<-list(Origine=Origine,Orcen=rescal$centres[Origine],Terminus=Terminus,Termcen=rescal$centres[Terminus])
}
####################################
# R�gression lin�aire et r�cup�ration des r�sidus - Detection des pentes
#
detect<-function(centresbis,cumulbis,nbwin,Origine,Terminus,interval,coef)
{
res<-1:nbwin
###
# Cas Terminus avant Origine
if (Terminus < Origine) {
 regot<-lsfit(centresbis[c(Origine:nbwin,1:Terminus)],cumulbis[c(Origine:nbwin,1:Terminus)])
 regto<-lsfit(centresbis[(Terminus+1):(Origine-1)],cumulbis[(Terminus+1):(Origine-1)])
 limite<-pente(regot$residuals,regto$residuals,centresbis[c(Origine:nbwin,1:Terminus)],centresbis[(Terminus+1):(Origine-1)],interval)
 trouvot<-cherche(regot$residuals,centresbis[c(Origine:nbwin,1:Terminus)],interval,coef,limite)
 trouvto<-cherche(regto$residuals,centresbis[(Terminus+1):(Origine-1)],interval,coef,limite)
 res[1:Terminus]<-regot$residuals[(nbwin-Origine+2):(nbwin-Origine+Terminus+1)]
 res[(Terminus+1):(Origine-1)]<-regto$residuals
 res[Origine:nbwin]<-regot$residuals[1:(nbwin-Origine+1)]
 lot<-length(trouvot[,1])
 lto<-length(trouvto[,1])
 trouv<-(1:(2*(lot+lto)))
 dim(trouv)<-c((lot+lto),2)
 for (i in 1:lot) {
  for (j in 1:2) {
   if (trouvot[i,j] > (nbwin-Origine+1)) trouv[i,j]<-trouvot[i,j]-nbwin+Origine-1
   else                                  trouv[i,j]<-trouvot[i,j]+Origine-1
  }
 }
 for (i in 1:lto) {
  for (j in 1:2) {
   trouv[i+lot,j]<-trouvto[i,j]+Terminus  
  }
 }
}
###
# Cas Origine avant Terminus
else {
 regot<-lsfit(centresbis[Origine:Terminus],cumulbis[Origine:Terminus])
 if ((Terminus!=nbwin) & (Origine!=1)) {
  regto<-lsfit(centresbis[c((Terminus+1):nbwin,1:(Origine-1))],cumulbis[c((Terminus+1):nbwin,1:(Origine-1))])
  limite<-pente(regot$residuals,regto$residuals,centresbis[Origine:Terminus],centresbis[c((Terminus+1):nbwin,1:(Origine-1))],interval)
  trouvot<-cherche(regot$residuals,centresbis[Origine:Terminus],interval,coef,limite)
  trouvto<-cherche(regto$residuals,centresbis[c((Terminus+1):nbwin,1:(Origine-1))],interval,coef,limite)
 }
 else if (Terminus!=nbwin) {
  regto<-lsfit(centresbis[(Terminus+1):nbwin],cumulbis[(Terminus+1):nbwin])
  limite<-pente(regot$residuals,regto$residuals,centresbis[Origine:Terminus],centresbis[(Terminus+1):nbwin],interval)
  trouvot<-cherche(regot$residuals,centresbis[Origine:Terminus],interval,coef,limite)
  trouvto<-cherche(regto$residuals,centresbis[(Terminus+1):nbwin],interval,coef,limite)
 }
 else if (Origine!=1) {
  regto<-lsfit(centresbis[1:(Origine-1)],cumulbis[1:(Origine-1)])
  limite<-pente(regot$residuals,regto$residuals,centresbis[Origine:Terminus],centresbis[1:(Origine-1)],interval)
  trouvot<-cherche(regot$residuals,centresbis[Origine:Terminus],interval,coef,limite)
  trouvto<-cherche(regto$residuals,centresbis[1:(Origine-1)],interval,coef,limite)
 }
 res[Origine:Terminus]<-regot$residuals
 if (Origine!=1) {
  res[1:(Origine-1)]<-regto$residuals[(nbwin-Terminus+1):(nbwin-Terminus+Origine-1)]
 }
 if (Terminus!=nbwin) {
  res[(Terminus+1):nbwin]<-regto$residuals[1:(nbwin-Terminus)]
 }
 lot<-length(trouvot[,1])
 lto<-length(trouvto[,1])
 trouv<-(1:(2*(lot+lto)))
 dim(trouv)<-c((lot+lto),2)
 for (i in 1:lot) {
  for (j in 1:2) {
   trouv[i,j]<-trouvot[i,j]+Origine-1
  }
 }
 for (i in 1:lto) {
  for (j in 1:2) {
   if (trouvto[i,j] > (nbwin-Terminus)) trouv[i+lot,j]<-trouvto[i,j]+Terminus-nbwin
   else                                 trouv[i+lot,j]<-trouvto[i,j]+Terminus  
  }
 }
}
detect<-list(res=res,trouv=trouv)
}
####################################
# Detection de regions sur tout le g�nome avec determination de l'origine et du terminus
#
regg<-function(x,pttnom,y,z,termtype,orforce,teforce,interval,coef)
{
tit<-titre(x)
###signification explicite du type de skew
if      (z==1) expsktype<-"T-A skew"
else if (z==2) expsktype<-"G-C skew"
else           expsktype<-"Erreur parametre sktype"
###signification explicite du type de terminus
if      (termtype==1) exptermtype<-"Minimum du G-C skew"
else if (termtype==2) exptermtype<-"Maximum du G-C skew"
else if (termtype==3) exptermtype<-"Minimum du T-A skew"
else if (termtype==4) exptermtype<-"Maximum du T-A skew"
else if (termtype==5) exptermtype<-"Maximum du Sense Analysis"
else if (termtype==6) exptermtype<-"Saisi manuellement"
else                  exptermtype<-"Erreur parametre termtype"
###
# Calculs des cumuls
rescal<-calcul(x,y,z)
nbwin<-rescal$nbwin
long<-rescal$long
###
# Recherche Terminus/Origine
if (termtype < 5) {
 chercheterm<-outerm(x,y,termtype)
 termin<-chercheterm$Termcen
 orig<-chercheterm$Orcen
 Terminus<-chercheterm$Terminus
 Origine<-chercheterm$Origine
}
else {
 if (termtype == 5) {
  ptt<-getPtt(pttnom)
  if (ptt$begin[1]>ptt$end[1]) ptt<-ajustptt(ptt,length(x))
  chercheterm<-outermgen(x,ptt,termtype)
  termin<-chercheterm$Termcen
  orig<-chercheterm$Orcen
 }
 else if (termtype == 6) {
  termin<-teforce
  orig<-orforce
 }
 Origine=trunc((orig/y)+1)
 Terminus=trunc((termin/y)+1)
}
###
# D�calage sur l'origine
cumulbis<-1:nbwin
centresbis<-1:nbwin
if (Origine==1) {
 centrebis<-rescal$centres
 cumulbis<-rescal$cumul
}
else {
 centresbis[Origine:nbwin]<-rescal$centres[Origine:nbwin]-(rescal$centres[Origine-1]+(y/2))
 centresbis[1:(Origine-1)]<-rescal$centres[1:(Origine-1)]+(long-rescal$centres[Origine-1]-(y/2))
 cumulbis[Origine:nbwin]<-rescal$cumul[Origine:nbwin]-rescal$cumul[Origine-1]
 cumulbis[1:(Origine-1)]<-rescal$cumul[1:(Origine-1)]+cumulbis[nbwin]
}
plot(rescal$centres,cumulbis,"l",main=tit,sub=expsktype,ylab="Cumul",xlab="Genome")
if (Origine!=1) lines(rescal$centres[(Origine-1):Origine],cumulbis[(Origine-1):Origine],type="l",col="white")
#write.table(cbind(rescal$centres,centresbis,rescal$tab,rescal$cumul,cumulbis),"data.txt",col.name=FALSE)
###
# R�gression lin�aire et r�cup�ration des r�sidus - Detection des pentes
resdet<-detect(centresbis,cumulbis,nbwin,Origine,Terminus,interval,coef)
trouv<-resdet$trouv
lt<-length(trouv[,1])
regions<-(1:(2*lt))
dim(regions)<-c(lt,2)
for (i in 1:lt) {
 for (j in 1:2) {
  regions[i,j]<-rescal$centres[trouv[i,j]]
 }
}
###
# Graphique des residus avec pentes detect�es en rouge
windows()
plot(rescal$centres,resdet$res,"l",main=tit,sub=paste(expsktype,"(Residual)"),ylab="Residus",xlab="Genome")
for (i in 1:lt) {
 if (trouv[i,2] >= trouv[i,1]) lines(rescal$centres[trouv[i,1]:trouv[i,2]],resdet$res[trouv[i,1]:trouv[i,2]],type="l",col="red")
 else                        { lines(rescal$centres[trouv[i,1]:nbwin],resdet$res[trouv[i,1]:nbwin],type="l",col="red")
                               lines(rescal$centres[1:trouv[i,2]],resdet$res[1:trouv[i,2]],type="l",col="red") } 
}
if (Origine !=1) lines(rescal$centres[(Origine-1):Origine],resdet$res[(Origine-1):Origine],type="l",col="white")
if (Terminus !=nbwin) lines(rescal$centres[Terminus:(Terminus+1)],resdet$res[Terminus:(Terminus+1)],type="l",col="white")
#write.table(regions,"toskout.txt",col.name=FALSE)
res<-cbind(rescal$centres,rescal$tab,resdet$res)
regg<-list(organism=tit,sktype=expsktype,termtype=exptermtype,window=y,ori=orig,term=termin,interval=interval,coef=coef,regions=regions,residuals=res)
}
####################################
# Detection de regions sur tout le g�nome sans determination de l'origine et du terminus
#
difg<-function(x,y,z,interval,coef) 
{
tit<-titre(x)
###signification explicite du type de skew
if      (z==1) expsktype<-"T-A skew"
else if (z==2) expsktype<-"G-C skew"
else           expsktype<-"Erreur parametre sktype"
###
exptermtype<-"Pas de terminus"
###
# Graphique des cumuls
rescal<-calcul(x,y,z)
nbwin<-rescal$nbwin
plot(rescal$centres,rescal$cumul,"l",main=tit,sub=expsktype,ylab="Cumul",xlab="Genome")
###
# Graphique des residus avec detection de pente
reg<-lsfit(rescal$centres,rescal$cumul)
res<-reg$residuals
limite<-pentesansterm(res,rescal$centres,interval)
trouv<-cherche(res,rescal$centres,interval,coef,limite)
lt<-length(trouv[,1])
regions<-(1:(2*lt))
dim(regions)<-c(lt,2)
for (i in 1:lt) {
 for (j in 1:2) {
  regions[i,j]<-rescal$centres[trouv[i,j]]
 }
}
windows()
plot(rescal$centres,res,"l",main=tit,sub=paste(expsktype,"(Residual)"),ylab="Residus",xlab="Genome")
for (i in 1:lt) {
 lines(rescal$centres[trouv[i,1]:trouv[i,2]],res[trouv[i,1]:trouv[i,2]],type="l",col="red")
 }
#write.table(regions,"toskout.txt",col.name=FALSE)
res<-cbind(rescal$centres,rescal$tab,res)
difg<-list(organism=tit,sktype=expsktype,termtype=exptermtype,window=y,interval=interval,coef=coef,regions=regions,residuals=res)
}
####################################
# Detection de regions sur tout le g�nome
#
tosk<-function(seq,pttnom="",sktype=1,termtype=2,fenetre=100,origine=1,terminus=1000,interval=100,coef=2) 
{
#seqlue<-read.fasta(File=seq)
#seqlue<-getSequence(seqlue[[1]])
if (termtype==0) retour<-difg(seq,fenetre,sktype,interval,coef)
else             retour<-regg(seq,pttnom,fenetre,sktype,termtype,origine,terminus,interval,coef)
tosk<-retour
}
####################################
# Traduction num�rique du brin : +1 ou -1
senscpt<-function(l,brin)
{
if (brin== " +"||brin== "+") senscpt <- l else senscpt <- -l
}
####################################
# Detection de regions par g�nes sans determination de l'origine et du terminus
#
genoosk<-function(x,pttnom,z,interval,coef) 
{
tit<-titre(x)
###signification explicite du type de skew
if      (z==1) expsktype<-"T-A skew"
else if (z==2) expsktype<-"G-C skew"
else if (z==3) expsktype<-"Sense Analysis"
else if (z==4) expsktype<-"T-A skew Codon 1"
else if (z==5) expsktype<-"T-A skew Codon 2"
else if (z==6) expsktype<-"T-A skew Codon 3"
else if (z==7) expsktype<-"G-C skew Codon 1"
else if (z==8) expsktype<-"G-C skew Codon 2"
else if (z==9) expsktype<-"G-C skew Codon 3"
else           expsktype<-"Erreur parametre sktype"
###
exptermtype<-"Pas de terminus"
###
# Graphique des cumuls
ptt<-getPtt(pttnom)
rescal<-calgen(x,ptt,z)
nbwin<-rescal$nbwin
plot(rescal$centres,rescal$cumul,"l",main=tit,sub=expsktype,ylab="Cumul",xlab="Genome")
###
# Graphique des residus avec detection de pente
reg<-lsfit(rescal$centrepis,rescal$cumul)
res<-reg$residuals
limite<-pentesansterm(res,rescal$centrepis,interval)
trouv<-cherche(res,rescal$centrepis,interval,coef,limite)
lt<-length(trouv[,1])
regions<-(1:(2*lt))
dim(regions)<-c(lt,2)
genesregion<-1:nbwin
genesregion[1:nbwin]<-0
for (i in 1:lt) {
 deb<-trouv[i,1]
 fin<-trouv[i,2]
 regions[i,1]<-ptt$begin[deb]
 regions[i,2]<-ptt$end[fin]  
 for (j in deb:fin) {
  genesregion[j]<-i
 }
}
genes<-data.frame(region=I(genesregion),begin=I(ptt$begin),end=I(ptt$end),gene=I(ptt$gene),product=I(ptt$product))
windows()
plot(rescal$centres,res,"l",main=tit,sub=paste(expsktype,"(Residual)"),ylab="Residus",xlab="Genome")
for (i in 1:lt) {
 lines(rescal$centres[trouv[i,1]:trouv[i,2]],res[trouv[i,1]:trouv[i,2]],type="l",col="red")
 }
genoosk<-list(organism=tit,sktype=expsktype,termtype=exptermtype,interval=interval,coef=coef,regions=regions,genes=genes)
}
####################################
# Detection de regions par g�nes avec determination de l'origine et du terminus
#
genotosk<-function(x,pttnom,z,termtyp,termwin,orforce,teforce,interval,coef)
{
tit<-titre(x)
###signification explicite du type de skew
if      (z==1) expsktype<-"T-A skew Genes"
else if (z==2) expsktype<-"G-C skew Genes"
else if (z==3) expsktype<-"Sense Analysis"
else if (z==4) expsktype<-"T-A skew Codon 1"
else if (z==5) expsktype<-"T-A skew Codon 2"
else if (z==6) expsktype<-"T-A skew Codon 3"
else if (z==7) expsktype<-"G-C skew Codon 1"
else if (z==8) expsktype<-"G-C skew Codon 2"
else if (z==9) expsktype<-"G-C skew Codon 3"
else           expsktype<-"Erreur parametre sktype"
###signification explicite du type de terminus
if      (termtyp==1) exptermtype<-"Minimum du G-C skew"
else if (termtyp==2) exptermtype<-"Maximum du G-C skew"
else if (termtyp==3) exptermtype<-"Minimum du T-A skew"
else if (termtyp==4) exptermtype<-"Maximum du T-A skew"
else if (termtyp==5) exptermtype<-"Maximum du Sense Analysis"
else if (termtyp==6) exptermtype<-"Saisi manuellement"
else                 exptermtype<-"Erreur parametre termtype"
###
# Calculs des cumuls
ptt<-getPtt(pttnom)
if (ptt$begin[1]>ptt$end[1]) ptt<-ajustptt(ptt,length(x))
rescal<-calgen(x,ptt,z)
nbwin<-rescal$nbwin
long<-rescal$long
###
# Recherche Terminus/Origine dans le genome
if (termtyp < 5) {
 chercheterm<-outerm(x,termwin,termtyp)
 termin<-chercheterm$Termcen
 orig<-chercheterm$Orcen
}
else if (termtyp == 5) {
 chercheterm<-outermgen(x,ptt,termtyp)
 termin<-chercheterm$Termcen
 orig<-chercheterm$Orcen
}
else if (termtyp == 6) {
 termin<-teforce
 orig<-orforce
}
###
# Recherche Terminus/Origine - Position dans le tableau
if (termtyp == 5) {
 Terminus<-chercheterm$Terminus
 Origine<-chercheterm$Origine
}
else {
 Origine<-1
 Terminus<-0
 i<-1
 while ((rescal$centres[i]<termin)&(i<nbwin)) {
  Terminus <- Terminus+1
  i <- i+1
 }
 if ((Terminus==0) | (rescal$centres[nbwin]<termin)) Terminus <- nbwin
 i<-1
 while ((rescal$centres[i]<orig)&(i<nbwin)) {
  Origine <- Origine+1
  i <- i+1
  if (rescal$centres[nbwin]<orig) Origine <- 1
 }
}
print(nbwin)
print(long)
print(orig)
print(Origine)
print(termin)
print(Terminus)
###
# D�calage sur l'origine
cumulbis<-1:nbwin
centresbis<-1:nbwin
if (Origine==1) {
 centresbis<-rescal$centrepis
 cumulbis<-rescal$cumul
}
else {
 ajustori<-trunc((ptt$end[Origine-1]-ptt$begin[Origine-1])/2)
 centresbis[Origine:nbwin]<-rescal$centrepis[Origine:nbwin]-(rescal$centrepis[Origine-1]+ajustori)
 centresbis[1:(Origine-1)]<-rescal$centrepis[1:(Origine-1)]+(long-rescal$centrepis[Origine-1]-ajustori)
 cumulbis[Origine:nbwin]<-rescal$cumul[Origine:nbwin]-rescal$cumul[Origine-1]
 cumulbis[1:(Origine-1)]<-rescal$cumul[1:(Origine-1)]+cumulbis[nbwin]
}
plot(rescal$centres,cumulbis,"l",main=tit,sub=expsktype,ylab="Cumul",xlab="Genome")
if (Origine!=1) lines(rescal$centres[(Origine-1):Origine],cumulbis[(Origine-1):Origine],type="l",col="white")
write.table(cbind(rescal$centres,centresbis,rescal$tab,rescal$cumul,cumulbis),"data.txt",col.name=FALSE)
###
# R�gression lin�aire et r�cup�ration des r�sidus - Detection des pentes
resdet<-detect(centresbis,cumulbis,nbwin,Origine,Terminus,interval,coef)
trouv<-resdet$trouv
lt<-length(trouv[,1])
regions<-(1:(2*lt))
dim(regions)<-c(lt,2)
genesregion<-1:nbwin
genesregion[1:nbwin]<-0
for (i in 1:lt) {
 deb<-trouv[i,1]
 fin<-trouv[i,2]
 if (deb != fin) {
  regions[i,1]<-ptt$begin[deb]
  regions[i,2]<-ptt$end[fin]  
  for (j in deb:fin) {
   genesregion[j]<-i
  }
 }
 else regions[i,]<-0 
}
genes<-data.frame(region=I(genesregion),begin=I(ptt$begin),end=I(ptt$end),gene=I(ptt$gene),product=I(ptt$product))
###
# Graphique des residus avec pentes detect�es en rouge
windows()
plot(rescal$centres,resdet$res,"l",main=tit,sub=paste(expsktype,"(Residual)"),ylab="Residus",xlab="Genome")
for (i in 1:lt) {
 if (trouv[i,2] >= trouv[i,1]) lines(rescal$centres[trouv[i,1]:trouv[i,2]],resdet$res[trouv[i,1]:trouv[i,2]],type="l",col="red")
 else                        { lines(rescal$centres[trouv[i,1]:nbwin],resdet$res[trouv[i,1]:nbwin],type="l",col="red")
                               lines(rescal$centres[1:trouv[i,2]],resdet$res[1:trouv[i,2]],type="l",col="red") } 
}
if (Origine !=1) lines(rescal$centres[(Origine-1):Origine],resdet$res[(Origine-1):Origine],type="l",col="white")
if (Terminus !=nbwin) lines(rescal$centres[Terminus:(Terminus+1)],resdet$res[Terminus:(Terminus+1)],type="l",col="white")
write.table(regions,"geskout.txt",col.name=FALSE)
genotosk<-list(organism=tit,sktype=expsktype,termtype=exptermtype,window=termwin,ori=orig,term=termin,interval=interval,coef=coef,regions=regions,genes=genes)
}
####################################
# Detection de regions par g�nes 
#
gesk<-function(seq,pttnom,sktype=1,termtype=2,fenetre=100,origine=1,terminus=1000,interval=10,coef=2)
{
#seqlue<-read.fasta(File=seq)
#seqlue<-getSequence(seqlue[[1]])
if (termtype==0) retour<-genoosk(seq,pttnom,sktype,interval,coef)
else             retour<-genotosk(seq,pttnom,sktype,termtype,fenetre,origine,terminus,interval,coef)
gesk<-retour
}

