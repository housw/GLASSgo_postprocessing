# author Jens Georg
# selects candidates for copraRNA based on a phylogenetic tree of the input sRNAs
# 

#call: R --slave -f  GLASSgo2CopraRNA_balanced_selection.r --args wildcard=NC_000913,NC_000911,NC_003197 exclude=NZ_CP009781.1,NZ_LN681227.1 max_number=20 outfile_prefix=sRNA ooi=NC_000913

args <- commandArgs(trailingOnly = TRUE)

ooi<-"NC_000913"
wildcard<-c("NC_000913","NC_000911","NC_003197","NC_016810","NC_000964","NC_002516","NC_003210","NC_007795","NC_003047")
max_number<-30
outfile_prefix<-"sRNA"
exclude<-c("NZ_CP009781.1","NZ_LN681227.1")


for(i in 1:length(args)){
	temp<-strsplit(args[i],"=")
	temp<-temp[[1]]
	temp1<-temp[1]
	temp2<-temp[2]
	assign(as.character(temp1),temp2)
 }

wil<-grep(",",wildcard)
if(length(wil)>0){
	wildcard<-strsplit(wildcard,",")[[1]]
} 
 
max_number<-as.numeric(max_number)

require(ape)
load("refined_GLASSgo_table.Rdata")


temp<-coor2
if(length(exclude)>0){
	temp_ex<-c()
	for(i in 1:length(exclude)){
		temp_ex1<-grep(exclude[i], coor2[,"fin"])
		if(length(temp_ex1)>0){
			temp_ex<-c(temp_ex,temp_ex1)
		}
	}
	if(length(temp_ex)>0){
	
	temp<-temp[-temp_ex,]
	}
}
coor2<-temp



# clustal omega call with kimura distance matrix for tree generation
clustalo3<-function(coor, positions){
	fasta<-c()
	for(i in 1:length(positions)){
		fasta<-c(fasta, paste(">",coor[positions[i],"fin"],sep=""))
		fasta<-c(fasta, as.character(coor[positions[i],"sequence"]))
	}
	write.table(fasta, file="temp_fasta", row.names=F, col.names=F, quote=F)
	wd<-getwd()
	command<-paste("clustalo -i ", "temp_fasta", " --distmat-out=distmatout.txt --full --output-order=input-order --use-kimura --force --max-hmm-iterations=-1", sep="")
	system(command)
	na<-grep(">", fasta)
	na<-gsub(">","",fasta[na])
	temp<-read.delim("distmatout.txt",sep="",header=F, , skip=1)
	unlink("distmatout.txt")
	unlink("temp_fasta")
	temp<-temp[,2:ncol(temp)]
	colnames(temp)<-na
	rownames(temp)<-na
	temp
}

# clustal omega call with percent identity matrix 
clustalo4<-function(coor, positions){
	wd<-getwd()
	fasta<-c()
	for(i in 1:length(positions)){
		fasta<-c(fasta, paste(">",coor[positions[i],"fin"],sep=""))
		fasta<-c(fasta, as.character(coor[positions[i],"sequence"]))
	}
	write.table(fasta, file="temp_fasta", row.names=F, col.names=F, quote=F)
	command<-paste("clustalo -i ", "temp_fasta", " --distmat-out=distmatout.txt --full --percent-id --output-order=input-order --force --max-hmm-iterations=-1", sep="")
	system(command)
	na<-grep(">", fasta)
	na<-gsub(">","",fasta[na])
	temp<-read.delim("distmatout.txt",sep="",header=F, , skip=1)
	unlink("distmatout.txt")
	temp<-temp[,2:ncol(temp)]
	colnames(temp)<-na
	rownames(temp)<-na
	temp
}

if(nrow(coor2)<max_number){
	fasta<-c()
	for(i in 1:nrow(coor2)){
		fasta<-c(fasta, paste(">",coor2[i,"fin"],sep=""))
		fasta<-c(fasta, as.character(coor2[i,"sequence"]))
	}
	nam<-paste(outfile_prefix,"CopraRNA_input_balanced.fasta", sep="_" )
	write.table(fasta, file=nam, row.names=F, col.names=F, quote=F)
}

if(nrow(coor2)>max_number){
	wildcard<-unique(c(ooi,wildcard))
	wild<-c()
	for(i in 1:length(wildcard)){
		wild<-c(wild,grep(wildcard[i],coor2[,"fin"]))
	}
	max_number<-max_number-length(wild)
	coortemp<-coor2
	if(length(wild)>0){
		coortemp<-coor2[-wild,]
	}
	dis<-clustalo3(coortemp, seq(1,nrow(coortemp)))
	dis<-as.dist(dis)
	clus<-(hclust(dis,method="average"))
	plot(clus)
	knum<-min(max_number,length(clus$labels)-1)
	if(knum<2){
		knum<-2
	}
	clus2<-rect.hclust(clus,k=knum)
	
	out<-c()
	for(i in 1:length(clus2)){
		temp<-clus2[[i]]
		temp_ooi<-grep(ooi,names(temp))
		if(length(temp_ooi)==0){
			temp2<-sample(length(temp),1)
			out<-c(out, names(temp)[temp2])
		}
	}
	out<-c(coor2[wild,"fin"],out)
	out<-match(out,coor2[,"fin"])
	fasta<-c()
	for(i in 1:length(out)){
		fasta<-c(fasta, paste(">",coor2[out[i],"fin"],sep=""))
		fasta<-c(fasta, as.character(coor2[out[i],"sequence"]))
	}
	nam<-paste(outfile_prefix,"CopraRNA_input_balanced.fasta", sep="_" )
	fasta<-gsub("\\..*","",fasta)
	write.table(fasta, file=nam, row.names=F, col.names=F, quote=F)


	dis<-clustalo3(coor2, seq(1,nrow(coor2)))
	dis<-as.dist(dis)
	clus<-(hclust(dis,method="average"))
	clus<-as.phylo(clus)
	lab<-clus$tip.label

	nam_selected<-match(coor2[out,"fin"],lab)
	nam_wildcard<-match(coor2[wild,"fin"],lab)
	nam_ooi<-match(coor2[wild[1],"fin"],lab)

	lab<-match(lab,coor2[,"fin"])
	lab<-coor2[lab,"nam2"]
	clus$tip.label<-lab
	nam<-paste(outfile_prefix,"tree_coprarna_candidates_balanced.pdf", sep="_" )
	pdf(nam)
	colo<-rep("1",length(lab))
	colo[nam_selected]<-"dodgerblue1"
	colo[nam_wildcard]<-"olivedrab2"
	colo[nam_ooi]<-"purple1"
	par(mar=c(3, 1, 1, 1), xpd=TRUE)
	plot(clus,tip.color=colo, cex=0.5 )

	legend("bottom",  inset=c(-0.05),bty="n", legend=c("organism of interst (ooi)","pre-selected organisms","selected organisms"), text.col=c("purple1","olivedrab2","dodgerblue1"),cex=0.6)
	par(xpd=FALSE)
	dev.off()

}

unlink("Rplots.pdf")