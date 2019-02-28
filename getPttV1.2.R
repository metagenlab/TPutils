getPtt<-function(x,strip.white=T) 
{
	testfic<-read.table(x,skip=3,nrows=1,sep="\t")
	if (is.numeric(testfic[1,3])) {
		ajust <- 0
		lectit<-try(read.table(x,nrows=1,sep="\t"))
	}
	else {
		ajust <- 1
		lectit<-try(read.table(x,skip=1,nrows=1,sep="\t"))
	}
	titre<-as.character(lectit[1,1])
	lect<-try(read.table(x,skip=3,sep="\t",quote=""))
	l<-(length(lect[,1])-ajust)
	champs<-(1:l)
	pos<-(1:l)
	debgen<-(1:l)
	fingen<-(1:l)
	brin<-(1:l)
	long<-(1:l)
	pid<-(1:l)
	nomgen<-(1:l)
	syn<-(1:l)
	cod<-(1:l)
	cog<-(1:l)
	prod<-(1:l)
	for (i in 1:l) {
		champs[i]<-as.character(lect[(i+ajust),1])
		if(strip.white){
			pos[i]<-strsplit(champs[i],split="\\..")
			debgen[i]<-as.numeric(pos[[i]][1])
			fingen[i]<-as.numeric(pos[[i]][2])
		}
		else{
			debgen[i] <- as.numeric(substr(champs[i],1,9))
			fingen[i] <- as.numeric(substr(champs[i],12,21))
		}
		brin[i]<-as.character(lect[(i+ajust),2])
		long[i]<-as.numeric(as.character(lect[(i+ajust),3]))
		pid[i]<-as.character(lect[(i+ajust),4])
		nomgen[i]<-as.character(lect[(i+ajust),5])
		syn[i]<-as.character(lect[(i+ajust),6])
		cod[i]<-as.character(lect[(i+ajust),7])
		cog[i]<-as.character(lect[(i+ajust),8])
		prod[i]<-as.character(lect[(i+ajust),9])
	}
	getPtt<-list(title=titre,begin=debgen,end=fingen,strand=brin,length=long,pid=pid,gene=nomgen,synonym=syn,code=cod,cog=cog,product=prod)
}



