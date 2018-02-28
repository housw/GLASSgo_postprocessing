
library(reshape)
require(seqinr)
require(phangorn)
mafft<-function(filename="ncrna.fa", outname="ncrna_aligned.fa"){
	command<-paste("mafft --maxiterate 1000 --localpair --quiet --inputorder ", filename, " > ", outname, sep="" )
	system(command)
	#fas<-read.fasta(outname)
	#fas
} 


closest<-function(x,y){
	out<-c()
	for(i in 1:length(x)){
		temp<- which.min(abs(y-x[i]))
		out<-c(out,temp)
	}
out
}


extend<-function(x){
	while(length(x)>length(unique(x))){
		x[duplicated(x)]<-x[duplicated(x)]+1
	}
x	
}

wd<-getwd()
d<-list.dirs( recursive=FALSE)
dis<-list()
for(i in 1:length(d)){
	setwd(d[i])
	
	load("refined_GLASSgo_table.Rdata")
	coor<-coor2
	fasta<-c()
	positions<-seq(1,nrow(coor))
	for(j in 1:length(positions)){
		fasta<-c(fasta, paste(">",coor[positions[j],"fin"],sep=""))
		fasta<-c(fasta, as.character(coor[positions[j],"sequence"]))
	}
	write.table(fasta, file="temp_fasta", row.names=F, col.names=F, quote=F)
	mafft(filename="temp_fasta")
	
	dat<-read.phyDat("ncrna_aligned.fa", format="fasta", type="DNA")
	dm <- dist.ml(dat, model="F81")
	treeNJ <- NJ(dm)
	fit = pml(treeNJ, data=dat)
	fitJC <- optim.pml(fit, model = "GTR")
	fit2<-fitJC
	fitJC<-fit2
	lab<-fitJC$tree$tip.label
	lab2<-lab
	#m <- as.matrix(dm)
	p<-cophenetic(fitJC$tree)
	temp<-p[,1]
	#names(temp)<-colnames(m)
	dis[[i]]<-temp
	
	num2<-25
	
	su<-summary(p[,1])
	
	su2<-p[which(p[,1]<su[4])]
	names(su2)<-rownames(p)[which(p[,1]<su[4])]
	su3<-p[which(p[,1]>=su[4])]
	names(su3)<-rownames(p)[which(p[,1]>=su[4])]
	
	su2<-sort(su2)
	clos<-su2[1:4]
	su2<-su2[5:length(su2)]
	num<-length(su2)*1.5
	if(num>60){
		num<-60
	}
	su2<-(log(su2))
	su2<-su2+abs(min(su2))
	su2<-sqrt(su2)
	#plot(sort(su2))
	#md<-max(dis[[i]])
	md<-max(su2)-min(su2)
	cu<-md/num
	rcu<-rep(cu,num)
	rcu<-min(su2)+cumsum(rcu)
	#rcu<-(rcu)
	#plot(sort(su2))
	#abline(h=(rcu))

	#rcu<-(rcu)^(1/10)
	
	#out<-c(clos,dis[[i]][closest(rcu,su2)])
	out<-c(clos,su2[closest(rcu,su2)])
	out<-unique(names(out))
	#out<-unique(out)
	out<-match(out,names(dis[[i]]))
	
	
	#out<-match(names(out), lab2)
	#out<-unique(out)
	m<-min(num2,length(out))
	out<-out[1:m]
	
	outn<-names(dis[[i]])[out]
	#out<-extend(out)
	#out<-out[1:num2]
	print(length(out))
	lab<-match(lab,coor2[,"fin"])
	lab<-coor2[lab,"nam2"]
	fitJC$tree$tip.label<-lab
	out2<-match(names(su3), lab2)
	ooi<-match(rownames(p)[1],lab2)
	
	colo<-rep("1",length(lab))
	colo[out]<-"deepskyblue2"
	colo[out2]<-"orangered2"
	colo[ooi]<-"olivedrab3"
	pdf("tree.pdf")
	plot(midpoint(fitJC$tree), tip.color=colo,cex=0.30, main=gsub("\\./","",d[i]))
	add.scale.bar()
	dev.off()
	
	
	# command<-paste("clustalo -i ", "temp_fasta", " --distmat-out=distmatout.txt --full --output-order=input-order --use-kimura --force --max-hmm-iterations=-1", sep="")
	# system(command)
	# na<-grep(">", fasta)
	# na<-gsub(">","",fasta[na])
	# temp<-read.delim("distmatout.txt",sep="",header=F, , skip=1)
	# unlink("distmatout.txt")
	# unlink("temp_fasta")
	# temp<-temp[,2:ncol(temp)]
	# colnames(temp)<-na
	# rownames(temp)<-na
	setwd(wd)
}
	

names(dis)<-d

