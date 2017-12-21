
# script to convert GLASSgo output to CopraRNA input
# author Jens Georg
# version 02-10-2017 4pm
# edited by Martin Raden 02-10-2017 5pm

#R --slave -f  ../GLASSgo2CopraRNA2.r --args filename=4083138.result sim=3 sim2=3 ooi=gdfb refpath=taxid_to_refseq cop_path=copra_refseq_positivelist.txt maxnumber=20 outfile_prefix=sRNA exclude=NZ_CP009513,NZ_CP010327 wildcard=NC_4678292092 full_evo=TRUE balanced=TRUE balanced_ooi=TRUE

# possible arguments:
# filename: name of the GLASSgo input File
# ooi: Refseq ID of organism of interest
# wildcard: list of always included genomes e.g. wildcard=NC_000913.3,NC_016810.1
# refpath: path to taxid to refseq file
# cop_path: path to copraRNA positive list
# maxnumber: maximal number of sequences included in the output fasta
# sim: maximal number of closely related organisms to the ooi for the  CopraRNA_input_ooi.fasta file
# sim2: maximal number of closely related organisms to the initial organisms for the CopraRNA_input_balanced_w_neighbours.fasta file
# outfile_prefix: prefix of the output files generated

args <- commandArgs(trailingOnly = TRUE)

wildcard<-c("NC_000913","NC_000911","NC_003197","NC_016810","NC_000964","NC_002516","NC_003210","NC_007795","NC_003047")
maxnumber<-20
sim<-3
sim2<-3
ooi<-c()

refpath<-"taxid_to_refseq"
cop_path<-"CopraRNA_available_organisms.txt"
outfile_prefix<-"sRNA"
exclude<-c()
full_evo<-TRUE
balanced<-TRUE
balanced_ooi<-TRUE
no_sims_for_preselected<-TRUE
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
 
ex<-grep(",",exclude)
if(length(ex)>0){
	exclude<-strsplit(exclude,",")[[1]]
}  
 
maxnumber=as.numeric(maxnumber)
sim<-as.numeric(sim)
sim2<-as.numeric(sim2)
 
#print(wildcard)
 require(ape)
 
glassgo2copra<-function(filename="results.txt",ooi="NC_000913.3",sim=4,sim2=3,wildcard=wildcard,refpath="taxid_to_refseq",cop_path="copra_refseq_positivelist.txt",maxnumber=20,outfile_prefix="sRNA", exclude=c(), full_evo=T, balanced=T, balanced_ooi=T, no_sims_for_preselected=T){


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
 # write.table(out, file="Accesions.txt", col.names = F, row.names = F, quote=F) # list of all Accesion numbers
 
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
  x<-gsub(".*\\:","",x)
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
  #write.table(y, file="organism_names.txt", col.names = F, row.names = F, quote=F) # list of all names
  #write.table(z, file="GI_list.txt", col.names = F, row.names = F, quote=F) # list of GIs
  out<-list(out,y,z)
  out
}






	#d<-dir()
	
	
	
	#d1<-grep("_verified$", d)
	
	#d3<-grep("_results$",d)
	#id<-gsub("_results.*","",d[d3])
	#fasta1<-read.delim(id,header=F)[,1]
	fasta<-read.delim(filename,header=F)[,1]
	
	direct<-grep("taxID:", fasta[1])
	if(length(direct)==0){
		fasta<-fasta[-c(1,2)]
	}
	#fasta<-c(as.character(fasta1[2]),as.character(fasta))
	##fasta<-as.character(read.delim(d[d3],header=F)[,1])


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

if(length(ooi)==0){
	ooi<-coor[1,"fin"]
}

ooi2<-grep(ooi,coor[,"fin"])

if(length(ooi2)==0){
	ooi<-coor[1,"fin"]
}

wildcard<-unique(c(ooi,wildcard))
	
	
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


# extract organisms from the wildcard (which should be in any case in the output)
pre<-c()
for(i in 1:length(wildcard)){
	pre_temp<-grep(wildcard[i],coor[,"fin"] )
	if(length(pre_temp)>0){
	pre<-c(pre,pre_temp)
	}
}


		cop_pre<-c()
		prestring<-c()
		if(length(pre)>0){
			cop_pre<-coor[pre,"fin"]
			prestring<-coor[pre,]
			#coor<-coor[-pre,]
		}

	if(length(prestring)>0){	
	if(is.matrix(prestring)==F){
		prestring<-t(as.matrix(prestring))
	}
	}


if(maxnumber<length(cop_pre)){
	maxnumber<-length(cop_pre)
}

fasta3<-c()
		for(i in 1:nrow(coor)){
			fasta3<-c(fasta3, paste(">",coor[i,"fin"],sep=""))
			fasta3<-c(fasta3, as.character(coor[i,"sequence"]))
			
		}

#fasta3<-c(as.character(fasta[2]),fasta3)
#fasta3<-c(as.character(fasta[1]),fasta3)
##fasta3<-c(as.character(fasta[2]),fasta3)
#fasta3<-c(as.character(fasta[1]),fasta3)		
write.table(fasta3, file="temp_fasta", row.names=F, col.names=F, quote=F)
#write.table(fasta3, file="temp_fasta2", row.names=F, col.names=F, quote=F)	


# calculate percent identity matrix
	
clustalo2<-function(coor){
	wd<-getwd()
	
	command<-paste("clustalo -i ", "temp_fasta", " --distmat-out=distmatout.txt --full --percent-id --output-order=input-order --force --max-hmm-iterations=-1", sep="")
	system(command)
	na<-c(coor[,1])
	
	#na<-c("query",coor[,1])
	temp<-read.delim("distmatout.txt",sep="",header=F, , skip=1)
	unlink("distmatout.txt")
	#temp<-temp[,2:ncol(temp)]
	colnames(temp)<-na
	rownames(temp)<-na
	temp
}
dis<-clustalo2(coor)


dis3<-dis
	
temp<-na.omit(match(colnames(dis)[2:ncol(dis)],coor[,1]))

# fasta4<-c()
		# for(i in 1:length(temp)){
			# fasta4<-c(fasta4, paste(">",coor[temp[i],"fin"],sep=""))
			# fasta4<-c(fasta4, as.character(coor[temp[i],"sequence"]))
			
		# }
# #write.table(fasta4, file="temp_fasta2", row.names=F, col.names=F, quote=F)				
	


### exclude 100% identicals
if(length(prestring)>0){
for(i in 1:length(pre)){
	g<-match(prestring[,1],colnames(dis))
	temp<-prestring[i,1]
	temp<-grep(temp, colnames(dis))
	temp1<-which(dis[temp,]==100)
	ex<-na.omit(match(g,temp1))
	if(length(ex)>0){
		temp1<-temp1[-ex]
	}
	
	
	if(length(temp1)>0){
		dis<-dis[-temp1,-temp1]
		
	}
}
}


i<-1
while(nrow(dis)>i){

	g<-c()
	if(length(prestring)>0){
		g<-match(prestring[,1],colnames(dis))
	}
	temp<-which(dis[i,]==100)
	ex<-na.omit(match(c(i,g), temp))
	if(length(ex)>0){
		temp<-temp[-ex]
	}
	
	if(length(temp)>0){
		if(i==1){
		
			if(colnames(dis)[1]=="query"){
				temp<-sort(temp)
				if(length(temp)>1){
					temp<-temp[2:length(temp)]
				}
				if(length(temp)==1){
					temp<-c()
				}
			}
		}
		
	}
	if(length(temp)>0){
		dis<-dis[-temp,-temp]
		
	}
	i<-i+1
}
temp<-na.omit(match(colnames(dis),coor[,1]))
if(length(exclude)>0){
	temp_ex<-c()
	for(i in 1:length(exclude)){
		temp_ex1<-grep(exclude[i], coor[,"fin"])
		if(length(temp_ex1)>0){
			temp_ex<-c(temp_ex,temp_ex1)
		}
	}
	if(length(temp_ex)>0){
	temp_ex<-intersect(temp_ex, temp)
	temp_ex<-match(temp_ex,temp)
	temp<-temp[-temp_ex]
}
}


temp_cop_ref<-na.omit(match(cop_pre, coor[temp,"fin"]))
cop_pre_w_ooi<-c()
if(length(cop_pre)>1){
	cop_pre_w_ooi<-cop_pre[2:length(cop_pre)]
}
temp_cop_ref_w_ooi<-na.omit(match(cop_pre_w_ooi, coor[temp,"fin"]))




temp2<-temp
if(length(temp_cop_ref)>0){
	temp2<-temp[-temp_cop_ref]
}
temp_5<-temp
if(length(temp_cop_ref_w_ooi)>0){
	temp_5<-temp[-temp_cop_ref_w_ooi]
}

fasta_wo_wild_w_ooi<-c()
		for(i in 1:length(temp_5)){
			fasta_wo_wild_w_ooi<-c(fasta_wo_wild_w_ooi, paste(">",coor[temp_5[i],"fin"],sep=""))
			fasta_wo_wild_w_ooi<-c(fasta_wo_wild_w_ooi, as.character(coor[temp_5[i],"sequence"]))
			
		}



fasta3<-c()
		for(i in 1:length(temp)){
			fasta3<-c(fasta3, paste(">",coor[temp[i],"fin"],sep=""))
			fasta3<-c(fasta3, as.character(coor[temp[i],"sequence"]))
			
		}
	
fasta_wo_wild<-c()
		for(i in 1:length(temp2)){
			fasta_wo_wild<-c(fasta_wo_wild, paste(">",coor[temp2[i],"fin"],sep=""))
			fasta_wo_wild<-c(fasta_wo_wild, as.character(coor[temp2[i],"sequence"]))
			
		}
write.table(fasta_wo_wild_w_ooi, file="temp_fasta_wo_wild_w_ooi.fasta", row.names=F, col.names=F, quote=F)
write.table(fasta_wo_wild, file="temp_fasta_wo_wild.fasta", row.names=F, col.names=F, quote=F)
write.table(fasta3, file="temp_fasta", row.names=F, col.names=F, quote=F)

fasta_unique<-fasta3 # unique CopraRNA "positive" sequences 
write.table(fasta_unique, file="un_fasta", row.names=F, col.names=F, quote=F)

copref<-read.delim(cop_path, sep="\t", header=T,comment.char = "#")	
nam<-c()
for(i in 1:nrow(coor[temp,])){
	tnam<-grep(gsub("\\..*","",coor[temp[i],"fin"]),copref[,1])
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
nam2<-paste(nam2,coor[temp,"fin"], sep="_")


# clustal omega call with kimura distance matrix for tree generation
clustalo3<-function(coor,fasta3, name="temp_fasta"){
	wd<-getwd()
	
	command<-paste("clustalo -i ", name, " --distmat-out=distmatout.txt --full --output-order=input-order --use-kimura --force --max-hmm-iterations=-1", sep="")
	system(command)
	na<-grep(">", fasta3)
	na<-gsub(">","",fasta3[na])
	temp<-read.delim("distmatout.txt",sep="",header=F, , skip=1)
	unlink("distmatout.txt")
	temp<-temp[,2:ncol(temp)]
	colnames(temp)<-na
	rownames(temp)<-na
	temp
}
# clustal omega call with percent identity matrix 
clustalo4<-function(coor){
	wd<-getwd()
	
	command<-paste("clustalo -i ", "temp_fasta", " --distmat-out=distmatout.txt --full --percent-id --output-order=input-order --force --max-hmm-iterations=-1", sep="")
	system(command)
	na<-grep(">", fasta3)
	na<-gsub(">","",fasta3[na])
	temp<-read.delim("distmatout.txt",sep="",header=F, , skip=1)
	unlink("distmatout.txt")
	temp<-temp[,2:ncol(temp)]
	colnames(temp)<-na
	rownames(temp)<-na
	temp
}
dis3<-clustalo4(coor)	



dis_uni<-clustalo3(coor, fasta_unique, name="un_fasta")
unlink("un_fasta")
name_uni<-colnames(dis_uni)
colnames(dis_uni)<-nam2
rownames(dis_uni)<-nam2
dis_uni<-as.dist(dis_uni)

clus_uni<-(hclust(dis_uni,method="average"))
clus_uni<-as.phylo(clus_uni)
lab<-clus_uni$tip.label



# if the number of "usable" homologs is smaller than the defined maximal number of homologs, all sequences are useed for the CopraRNA input
if(ncol(dis3)<=maxnumber){	
	out<-colnames(dis3)
	
	
	out<-match(paste(">",out,sep=""),fasta3)
	fasta2<-c()
		for(i in 1:length(out)){
			fasta2<-c(fasta2, as.character(fasta3[out[i]]))
			fasta2<-c(fasta2, as.character(fasta3[out[i]+1]))
			
		}
	nam<-paste(outfile_prefix,"CopraRNA_input.fasta", sep="_" )
	write.table(fasta2, file=nam, row.names=F, col.names=F, quote=F)
	
}

# if the number of "usable" homologs is bigger than the defined maximal number of homologs, sequences are selected based on a UPGMA tree
if(ncol(dis3)>maxnumber){
	
# (MARTIN) ensure 'dis2' is global to enable access below
dis2<<-clustalo3(coor, fasta_wo_wild, name="temp_fasta_wo_wild.fasta")
unlink("temp_fasta_wo_wild.fasta")

ooil<-sort(dis3[grep(ooi, colnames(dis3)),], decreasing=T)

# (MARTIN) ensure 'dis2' is global to enable access below
dis2<<-as.dist(dis2)
clus<-(hclust(dis2,method="average"))

plot(clus)

knum<-min(maxnumber-length(cop_pre),length(clus$labels)-1)
if(knum<2){
	knum<-2
}
clus2<-rect.hclust(clus,k=knum)







# wclus<-c()
 # for(i in 1:length(clus2)){
	
	 # temp<-clus2[[i]]
	 # #wil<-na.omit(match(cop_pre,names(temp)))
	 # #temp2<-sample(length(temp),1)
	 # #if(length(wil)>0){
		# # wclus<-c(wclus,i)
	 # #}
 # }



# while(maxnumber-(length(cop_pre))!=knum){
# #print(knum)
# knum<-knum-1
# pdf("temp.pdf")
# plot(clus)

# #knum<-min(maxnumber,length(clus$labels)-1)
# clus2<-rect.hclust(clus,k=knum)

# dev.off()
# wclus<-c()
# for(i in 1:length(clus2)){
	
	# temp<-clus2[[i]]
	# wil<-na.omit(match(cop_pre,names(temp)))
	# #temp2<-sample(length(temp),1)
	# if(length(wil)>0){
		# wclus<-c(wclus,i)
	# }
	
# }


# }





##########


out<-c()
for(i in 1:length(clus2)){
	
	temp<-clus2[[i]]
	wil<-na.omit(match(cop_pre,names(temp)))
	temp2<-sample(length(temp),1)
	#if(length(wil)>0){
	#	temp2<-wil
	#}
	out<-c(out, names(temp)[temp2])
}

# wil<-match(cop_pre,out)
# if(length(wil)>0){
	# out<-c(cop_pre,out[-wil])
# }
# if(length(out)>maxnumber){
# out<-out[1:maxnumber]
# }
# out1<-out
# }

out<-c(out,cop_pre)
out1<-out
}

if(balanced==T){



out<-match(paste(">",out,sep=""),fasta3)
fasta2<-c()
		for(i in 1:length(out)){
			fasta2<-c(fasta2, as.character(fasta3[out[i]]))
			fasta2<-c(fasta2, as.character(fasta3[out[i]+1]))
			
		}
		write.table(fasta2, file="temp_fasta", row.names=F, col.names=F, quote=F)
dis_balanced<-clustalo3(coor, fasta2)
fasta2<-gsub("\\..*","",fasta2)
unlink("temp_fasta")

nam_overlap<-c()
for(i in 1:length(colnames(dis_balanced))){
	t1<-grep(colnames(dis_balanced)[i],lab)
	if(length(t1)>0){
		nam_overlap<-c(nam_overlap,t1)
	}
}
nam_overlap2<-c()
for(i in 1:length(cop_pre)){
	t1<-grep(cop_pre[i],lab)
	if(length(t1)>0){
		nam_overlap2<-c(nam_overlap2,t1)
	}
}
nam_ooi<-c()
for(i in 1:length(ooi)){
	t1<-grep(ooi[i],lab)
	if(length(t1)>0){
		nam_ooi<-c(nam_ooi,t1)
	}
}
dis_balanced<-as.dist(dis_balanced)
clus<-(hclust(dis_balanced,method="average"))

nam<-paste(outfile_prefix,"tree_coprarna_candidates_balanced.pdf", sep="_" )




pdf(nam)
colo<-rep("1",length(name_uni))
colo[nam_overlap]<-"dodgerblue1"
colo[nam_overlap2]<-"olivedrab2"
colo[nam_ooi]<-"purple1"
par(mar=c(3, 1, 1, 1), xpd=TRUE)
plot(clus_uni,tip.color=colo, cex=0.5 )

legend("bottom",  inset=c(-0.05),bty="n", legend=c("organism of interst (ooi)","pre-selected organisms","selected organisms"), text.col=c("purple1","olivedrab2","dodgerblue1"),cex=0.6)
par(xpd=FALSE)
dev.off()


nam<-paste(outfile_prefix,"CopraRNA_input_balanced.fasta", sep="_" )
write.table(fasta2, file=nam, row.names=F, col.names=F, quote=F)
}

############### balanced with close organisms for all start organisms

# (MARTIN) ensure dis2 exists

if(full_evo==T){
if(ncol(dis3)>maxnumber){
sel<-out1

if(no_sims_for_preselected==T){
	tempooi<-match(ooi,cop_pre)
	cop_pre2<-cop_pre[-tempooi]
	
	t2<-na.omit(match(cop_pre2,sel))
	if(length(t2)>0){
		sel<-sel[-t2]
	}
	
	
}


write.table(fasta3, file="temp_fasta", row.names=F, col.names=F, quote=F)
dis2<<-clustalo3(coor, fasta3, name="temp_fasta")
unlink("temp_fasta")
dis2<-as.dist(dis2)

clus<-(hclust(dis2,method="average"))
for(i in 1:length(sel)){
	#print(i)
	sim3<-sim2
	ooi_t<-sel[i]
	ooi_tl<-sort(dis3[grep(ooi_t, colnames(dis3)),], decreasing=T)
	ident<-which(ooi_tl==100)
	if(length(ident)>0){
		ooi_tl<-ooi_tl[-ident]
	}
	
close_orgs<-c()
iii<-0
while((length(close_orgs)+1)<sim3){
plot(clus)

knum2<-min(maxnumber,length(clus$labels)-1)-iii
if(knum2<2){
break
}
#print(knum2)
clus3<-rect.hclust(clus,k=knum2)

ooi_t1<-which(names(unlist(clus3))==ooi_t)
len<-as.numeric(summary(clus3)[,1])
su<-0
ii<-1
while(su<ooi_t1){
	su<-su+len[ii]
	ii<-ii+1
}
ii<-ii-1



close_orgs<-sort(ooi_tl[intersect(names(clus3[[ii]]),names(ooi_tl))],decreasing=T)
iii<-iii+1
}

# temp_full_evo<-match(names(close_orgs), coor[,"fin"])
# fasta_temp<-c()
		# for(jj in 1:length(temp_full_evo)){
			# fasta_temp<-c(fasta_temp, paste(">",coor[temp_full_evo[jj],"fin"],sep=""))
			# fasta_temp<-c(fasta_temp, as.character(coor[temp_full_evo[jj],"sequence"]))
			
		# }
# write.table(fasta_temp, file="temp_fasta", row.names=F, col.names=F, quote=F)
# dist_temp<-clustalo3(coor,fasta_temp, name="temp_fasta")	
# unlink("temp_fasta")
# dist_temp<-as.dist(dist_temp)	
if(length(close_orgs)>sim3){
	 n<-length(close_orgs)%/%sim3
	 n<-seq(1,n*sim3,by=n)
	 n<-names(close_orgs)[n]
	 sel<-unique(c(sel,n))
 }
# clus_temp<-hclust(dist_temp,method="average")
# plot(clus_temp)
# clus_temp<-rect.hclust(clus,k=sim2)
# out_temp<-c()
# for(jj in 1:length(clus_temp)){
	
	# temp4<-clus_temp[[jj]]
	# wil<-na.omit(match(cop_pre,names(temp)))
	# temp5<-sample(length(temp4),1)
	# if(length(wil)>0){
		# temp5<-wil
	# }
	# out_temp<-c(out_temp, names(temp4)[temp5])
# }
# sel<-unique(c(sel,out_temp))
if(length(close_orgs)<=sim3){
sel<-unique(c(sel,names(close_orgs)))

}
}


out<-unique(c(sel,cop_pre))


out<-match(paste(">",out,sep=""),fasta3)
fasta2<-c()
		for(i in 1:length(out)){
			fasta2<-c(fasta2, as.character(fasta3[out[i]]))
			fasta2<-c(fasta2, as.character(fasta3[out[i]+1]))
			
		}
write.table(fasta2, file="temp_fasta", row.names=F, col.names=F, quote=F)
dis_balanced<-clustalo3(coor, fasta2,nam="temp_fasta")
unlink("temp_fasta")


nam_overlap<-c()
for(i in 1:length(colnames(dis_balanced))){
	t1<-grep(colnames(dis_balanced)[i],lab)
	if(length(t1)>0){
		nam_overlap<-c(nam_overlap,t1)
	}
}
nam_overlap2<-c()
for(i in 1:length(out1)){
	t1<-grep(out1[i],lab)
	if(length(t1)>0){
		nam_overlap2<-c(nam_overlap2,t1)
	}
}
nam_overlap3<-c()
for(i in 1:length(cop_pre)){
	t1<-grep(cop_pre[i],lab)
	if(length(t1)>0){
		nam_overlap3<-c(nam_overlap3,t1)
	}
}
nam_ooi<-c()
for(i in 1:length(ooi)){
	t1<-grep(ooi[i],lab)
	if(length(t1)>0){
		nam_ooi<-c(nam_ooi,t1)
	}
}
fasta2<-gsub("\\..*","",fasta2)




#dis_balanced<-as.dist(dis_balanced)
#clus<-(hclust(dis_balanced,method="average"))

nam<-paste(outfile_prefix,"tree_coprarna_candidates_balanced_w_neighbours.pdf", sep="_" )
pdf(nam)
colo<-rep("1",length(name_uni))
colo[nam_overlap]<-"orangered"
colo[nam_overlap2]<-"dodgerblue1"
colo[nam_overlap3]<-"olivedrab2"
colo[nam_ooi]<-"purple1"
par(mar=c(3, 1, 1, 1), xpd=TRUE)
plot.phylo(clus_uni,tip.color=colo, cex=0.5 )

legend("bottom",  inset=c(-0.05),bty="n", legend=c("organism of interst (ooi)","pre-selected organisms","selected organisms","similar_organisms"), text.col=c("purple1","olivedrab2","dodgerblue1","orangered"),cex=0.6)
par(xpd=FALSE)
dev.off()








nam<-paste(outfile_prefix,"CopraRNA_input_balanced_w_neighbours.fasta", sep="_" )
write.table(fasta2, file=nam, row.names=F, col.names=F, quote=F)


} # (MARTIN) ensure dis2 exists

}



#############





if(balanced_ooi==T){
if(ncol(dis3)>maxnumber){



dis2<<-clustalo3(coor, fasta_wo_wild_w_ooi, name="temp_fasta_wo_wild_w_ooi.fasta")
dis2<<-as.dist(dis2)
unlink("temp_fasta_wo_wild_w_ooi.fasta")


sel<-c()
ident<-which(ooil==100)
	if(length(ident)>0){
		ooil<-ooil[-ident]
	}

clus<-(hclust(dis2,method="average"))

plot(clus)
close_orgs<-c()
iii<-0
while(length(close_orgs)<sim){




knum2<-min(maxnumber,length(clus$labels)-1)-iii
if(knum2<2){
break
}
#print(knum2)
clus3<-rect.hclust(clus,k=knum2)

ooi1<-grep(ooi, names(unlist(clus3)))
#print(ooi1)
len<-as.numeric(summary(clus3)[,1])
su<-0
ii<-1
while(su<ooi1){
	su<-su+len[ii]
	ii<-ii+1
}
ii<-ii-1


#close_orgs<-names(sort(ooil[intersect(names(clus3[[ii]]),names(ooil))],decreasing=T))
close_orgs<-sort(ooil[intersect(names(clus3[[ii]]),names(ooil))],decreasing=T)

	
iii<-iii+1
}


if(length(close_orgs)>sim){
	n<-length(close_orgs)%/%sim
	n<-seq(1,n*sim,by=n)
	n<-names(close_orgs)[n]
	sel<-unique(c(sel,n))
}


if(length(close_orgs)<=sim){
sel<-unique(c(sel,names(close_orgs)))

}






wildcard2<-unique(c(wildcard,sel))


pre<-c()
for(i in 1:length(wildcard2)){
	pre<-c(pre, grep(wildcard2[i], coor[,"fin"]))
	
}




		cop_pre<-c()
		prestring<-c()
		if(length(pre)>0){
			cop_pre<-coor[pre,"fin"]
			prestring<-coor[pre,]
			#coor<-coor[-pre,]
		}

	if(length(prestring)>0){	
	if(is.matrix(prestring)==F){
		prestring<-t(as.matrix(prestring))
	}
	}


	
maxnumber2<-maxnumber-length(pre)
knum<-min(maxnumber2,length(clus$labels)-1)

if(knum<2){
	knum<-2
}



temp_cop<-na.omit(match(wildcard2, coor[temp,"fin"]))






if(length(temp_cop)>0){
	temp_5<-temp[-temp_cop]
}

fasta_wo_wild<-c()
		for(i in 1:length(temp_5)){
			fasta_wo_wild<-c(fasta_wo_wild_w_ooi, paste(">",coor[temp_5[i],"fin"],sep=""))
			fasta_wo_wild<-c(fasta_wo_wild_w_ooi, as.character(coor[temp_5[i],"sequence"]))
			
		}

write.table(fasta_wo_wild, file="temp_fasta_wo_wild.fasta", row.names=F, col.names=F, quote=F)




dis2<<-clustalo3(coor, fasta_wo_wild, name="temp_fasta_wo_wild.fasta")
dis2<<-as.dist(dis2)
unlink("temp_fasta_wo_wild.fasta")



clus<-hclust(dis2, method="average")

plot(clus)

clus2<-rect.hclust(clus,k=knum)




out<-c()
for(i in 1:length(clus2)){
	
	temp<-clus2[[i]]
	#wil<-na.omit(match(cop_pre,names(temp)))
	temp2<-sample(length(temp),1)
	# if(length(wil)>0){
		# temp2<-wil
	# }
	out<-c(out, names(temp)[temp2])
}

 # wil<-match(cop_pre,out)
 # if(length(wil)>0){
	 # out<-c(cop_pre,out[-wil])
 # }

 out<-c(out,prestring[,"fin"])
 

out<-match(paste(">",out,sep=""),fasta3)
fasta2<-c()
		for(i in 1:length(out)){
			fasta2<-c(fasta2, as.character(fasta3[out[i]]))
			fasta2<-c(fasta2, as.character(fasta3[out[i]+1]))
			
		}


write.table(fasta2, file="temp_fasta", row.names=F, col.names=F, quote=F)
dis_ooi<-clustalo3(coor, fasta2)
nam_overlap<-c()
for(i in 1:length(colnames(dis_ooi))){
	t1<-grep(colnames(dis_ooi)[i],lab)
	if(length(t1)>0){
		nam_overlap<-c(nam_overlap,t1)
	}
}
nam_overlap2<-c()
for(i in 1:length(wildcard)){
	t1<-grep(wildcard[i],lab)
	if(length(t1)>0){
		nam_overlap2<-c(nam_overlap2,t1)
	}
}
nam_overlap3<-c()
for(i in 1:length(sel)){
	t1<-grep(sel[i],lab)
	if(length(t1)>0){
		nam_overlap3<-c(nam_overlap3,t1)
	}
}
nam_ooi<-c()
for(i in 1:length(ooi)){
	t1<-grep(ooi[i],lab)
	if(length(t1)>0){
		nam_ooi<-c(nam_ooi,t1)
	}
}
fasta2<-gsub("\\..*","",fasta2)

dis_ooi<-clustalo3(coor, fasta2)
unlink("temp_fasta")

#dis_ooi<-as.dist(dis_ooi)
#clus<-(hclust(dis_ooi,method="average"))
nam<-paste(outfile_prefix,"tree_coprarna_candidates_ooi.pdf", sep="_" )

pdf(nam)
colo<-rep("1",length(name_uni))

colo[nam_overlap]<-"dodgerblue1"
colo[nam_overlap3]<-"orangered"       

colo[nam_overlap2]<-"olivedrab2"
colo[nam_ooi]<-"purple1"
par(mar=c(3, 1, 1, 1), xpd=TRUE)
plot(clus_uni,tip.color=colo, cex=0.5 )

legend("bottom",  inset=c(-0.05),bty="n", legend=c("organism of interst (ooi)","pre-selected organisms","selected organisms","similar_organisms"), text.col=c("purple1","olivedrab2","dodgerblue1","orangered"),cex=0.6)
par(xpd=FALSE)
dev.off()













nam<-paste(outfile_prefix,"tree_coprarna_candidates_ooi.fasta", sep="_" )

write.table(fasta2, file=nam, row.names=F, col.names=F, quote=F)
}
}
unlink("temp.pdf")
unlink("Rplots.pdf")
}











############################################
# call function with given default values (from top of the file) or overwritten values from commandline call
############################################

glassgo2copra(filename=filename,sim=sim,sim2=sim2,ooi=ooi,wildcard=wildcard, refpath=refpath, maxnumber=maxnumber, cop_path=cop_path,outfile_prefix=outfile_prefix, exclude=exclude, full_evo=full_evo, balanced=balanced, balanced_ooi=balanced_ooi,no_sims_for_preselected=no_sims_for_preselected)



# exit script (shouldnt be needed)
q();





























