
# author Jens Georg
# The script extracts all sequences which are available for CopraRNA and writes them in a fasta file
# Dependencies: "taxid_to_refseq" and "CopraRNA_available_organisms.txt"

#call: R --slave -f  GLASSgo_2_copraRNA_fasta_generation.r --args filename=4083138.result refpath=taxid_to_refseq cop_path=copra_refseq_positivelist.txt outfile=coprarna_candidates.txt


args <- commandArgs(trailingOnly = TRUE)

refpath<-"taxid_to_refseq"
cop_path<-"CopraRNA_available_organisms.txt"
filename<-"4083138.result"
outfile<-"coprarna_candidates.txt"

 for(i in 1:length(args)){
	temp<-strsplit(args[i],"=")
	temp<-temp[[1]]
	temp1<-temp[1]
	temp2<-temp[2]
	assign(as.character(temp1),temp2)
 }




to_refseq2<-function(result){
		load(refpath)
		result<-gsub("\\..*","",result)
		temp<-c()
		
		for(i in 1:length(result)){
		temp1<-grep(result[i],ref[,"full_genome_entry"])
		temp<-c(temp,temp1[1])
		}
		temp<-ref[temp,"Chromosomes.RefSeq"]
		refseq<-temp
		refseq
	}
	
	fill<-function(x){
		temp<-which(x[,ncol(x)]=="no_name")
		out<-x[,ncol(x)]
		count<-1
		while(length(temp)>0){
			out[temp]<-x[temp,(ncol(x)-count)]
			temp<-which(out=="no_name")
		}
		out
	}

export_accesions<-function(x){ # x = copy pasted text from GLASSgo fasta output
  temp<-grep(">", x)
  x<-x[temp]
  x<-strsplit(x,"\\|")
  out<-c()
  for(i in 1:length(x)){
    #if(length(x[[i]])>2){
      #out<-c(out,x[[i]][2]) # for(RNAlien)
	  out<-c(out,x[[i]][4]) # forGLASgo
    #}
  }
  if(is.matrix(out)==F){
	out<-t(as.matrix(out))
  }
 
 
  out<-matrix(out,length(out),1)
  colnames(out)<-"Accesion_number"
  out
}

export_ncRNA_coordinates<-function(x){ # x = copy pasted text from GLASSgo fasta output
  temp<-grep(">", x)
  header<-as.character(x[temp[1:length(temp)]])
  seqs<-as.character(x[temp[1:length(temp)]+1])
  x<-x[temp]
  x<-x[1:length(x)]
  y<-x
  z<-x
  z<-gsub(">gb\\|","",z)
  z<-gsub("\\|.*","",z)
  y<-gsub(".*[0-9]{1,}-[0-9]{1,} ","", y)
  y<-gsub(",.*","",y)
  y<-gsub("genome assembly","",y)
  y<-gsub("complete genome","",y)
  y<-gsub("chromosome","",y)
  x<-gsub(".*\\|[:ABCDEFGHIJKLMNOPQRSTUVWXYZ:]*[:0123456789:]*\\:","",x)
  x<-gsub(".*\\|\\:","",x)
  x<-gsub(" .*","", x)
  out<-matrix(,length(x),6)
  colnames(out)<-c("Strand","start","end","name","Full_header","sequence")
  out[,1]<-"+"
  for(i in 1:length(x)){
    temp<-grep("c",x[i])
    if(length(temp)==1){
      out[i,1]<-"-"
    }
    temp<-gsub("c","",x[i])
    temp<-strsplit(temp,"-")
    out[i,2]<-temp[[1]][1]
    out[i,3]<-temp[[1]][2]
  }
  out[,4]<-y
  out[,5]<-as.character(header)
  out[,6]<-as.character(seqs)
  
  out<-list(out,y,z)
  out
}




require(seqinr)
x<-read.fasta(filename)
write.fasta(x, file.out=filename, nbchar=100000, names=names(x))

fasta<-read.delim(filename,header=F)[,1]
	
direct<-grep("taxID:", fasta[1])
if(length(direct)==0){
	fasta<-fasta[-c(1,2)]
}
	

nn<-grep(">",fasta)
nn1<-grep("n",fasta[nn+1], ignore.case=T)
if(length(nn1)>0){
fasta<-fasta[-c(nn[nn1],nn[nn1]+1)]
}
ac<-export_accesions(as.character(fasta))  

co<-export_ncRNA_coordinates(fasta)
coor<-cbind(ac,co[[1]])


#if more than 1 homolog is detected for one organism or the same Refseq ID, keep only the homolog with the highest identity to the input

iden<-gsub(".*VAL:","",coor[,"Full_header"])
iden<-as.numeric(gsub("%.*","",iden))
coor<-cbind(coor,iden)
coor<-coor[order(iden, decreasing=T), ]

dup<-which(duplicated(coor[,1]))
if(length(dup)>0){
coor<-coor[-dup,]
}
coor2<-coor
ac<-coor[,1]
coor[,1]<-paste(coor[,1],coor[,3], sep="_")




taxi<-gsub(".*taxID:","",coor[,"Full_header"])

fin<-to_refseq2(ac)
coor<-cbind(coor,taxi,fin)



dup<-which(duplicated(fin))

if(length(dup)>0){
coor<-coor[-dup,]
coor2<-coor2[-dup,]
}


	if(is.matrix(coor)==T){
		na<-which((coor[,"fin"])=="-")
		if(length(na)>0){
			coor<-coor[-na,]
		}
	}
	
	
	if(is.matrix(coor)==T){
		na<-which(is.na(coor[,"fin"]))
		if(length(na)>0){
			coor<-coor[-na,]
		}
	}


	
	
# keep only homologs represented in the coprarna reference file	
	
copref<-read.delim(cop_path, sep="\t", header=T,comment.char = "#")	
	se<-function(x){
		out<-grep(x, copref[,1])[1]
		if(length(out)==0){
			out<-NA
		}
		out
	}
notinlist<-which(is.na(unlist(lapply(gsub("\\..*","",coor[,"fin"]),se))))

if(length(notinlist)>0){
	coor<-coor[-notinlist,]	
}


fasta3<-c()
		for(i in 1:nrow(coor)){
			fasta3<-c(fasta3, paste(">",coor[i,"fin"],"|",gsub(">","",coor[i,"Full_header"]),sep=""))
			fasta3<-c(fasta3, as.character(coor[i,"sequence"]))
			
		}
		
write.table(fasta3, file=outfile, row.names=F, col.names=F, quote=F)

copref<-read.delim(cop_path, sep="\t", header=T,comment.char = "#")	
nam<-c()
for(i in 1:nrow(coor)){
	tnam<-grep(gsub("\\..*","",coor[i,"fin"]),copref[,1])
	nam<-c(nam,as.character(copref[tnam,2]))
	
	
}

nam2<-c()
for(i in 1:length(nam)){
	temp1<-substr(nam[i],1,3)
	temp2<-strsplit(nam[i],"_")[[1]]
	temp1<-paste(temp1,"_",temp2[2], sep="")
	if(length(temp2)>2){
		temp1<-paste(temp1, temp2[length(temp2)], sep="_")
	}
	nam2<-c(nam2,temp1)
}
nam2<-paste(nam2,coor[,"fin"], sep="_")

coor<-cbind(coor,nam2)


save(coor, file="full_GLASSgo_table.Rdata")
