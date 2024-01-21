options(stringsAsFactors=FALSE)
library(stringr)
library(plyr)
library(pheatmap)
library(ComplexHeatmap)
library(zoo)

# set work directory
wor_dire=
setwd(wor_dire)
for (prot in c('SRRM2'))
{

# Average the amino acid charge every: 
aver=20

# Display amino acid position every;
disp=100

# Read the fasta file (saved in the same dire)

fast=read.table(paste(prot,'.fasta',sep=""),header=TRUE,sep="\t",check.names=FALSE);colnames(fast)='f'
seri=as.data.frame(unlist(strsplit(fast$f,"")))
colnames(seri)=c('a')
seri=as.data.frame(seri[seri$a!="*",])
colnames(seri)=c('a')

amin=unique(seri)

# Assign the charge to each residue in the fasta file
# charges from Emboss http://www.bioinformatics.nl/cgi-bin/emboss

acharge=vector()
acharge["A"]=0.000
acharge["C"]=0.000
acharge["D"]=-1.000
acharge["E"]=-1.000
acharge["F"]=0.000
acharge["G"]=0.000
acharge["H"]=0.500
acharge["I"]=0.000
acharge["K"]=1.000
acharge["L"]=0.000
acharge["M"]=0.000
acharge["N"]=0.000
acharge["P"]=0.000
acharge["Q"]=0.000
acharge["R"]=1.000
acharge["S"]=0.000
acharge["T"]=0.000
acharge["V"]=0.000
acharge["W"]=0.000
acharge["Y"]=0.000

# Assign the charge to each amino acid in the fasta file
chfi=cbind.data.frame(as.numeric(rownames(seri)),seri,acharge[seri$a])
colnames(chfi)=c('Position','Residue','Charge')

# Compute the running average
ma=rollmean(chfi$Charge, aver,fill = list(NA, NULL, NA))

# Order of the residues in the heatmap

orde=c("C","W","Y","F","M","L","I","V","A","G","P","Q","N","T","S","E","D","K","H","R")

# Build the heatmap of the residues

w=1
slid=data.frame()
for(i in seq(1,nrow(seri),w))
{
 comp=as.data.frame(as.data.frame(seri[i:((i+w)-1),]))
 d=count(comp)
 colnames(d)[1]=c("x")
 c=merge(amin$a,d,by.y='x',all.x=TRUE)
 if(nrow(c[is.na(c$freq),])>0){c[is.na(c$freq),]$freq=0}
 e=as.data.frame(t(c))
 f=e[2,];colnames(f)=e[1,];rownames(f)=i
 slid=rbind.data.frame(slid,f)
}

r=as.matrix(sapply(slid[,], as.numeric)); rownames(r)=rownames(slid)
r=r[,match(orde, colnames(r))] #,nomatch=0)]
labe=as.data.frame(rownames(r));colnames(labe)="c"
labe[as.numeric(labe$c)%%(disp)!=0,]=""
rownames(r)=labe$c
colnames(r)=orde
if(length(r[is.na(r)])>0){r[is.na(r)]=0}


# charge of the amino acid at the side 
sa=as.data.frame(acharge[colnames(r)])
colnames(sa)=c("c")
# Printing session

library(circlize)
col_fun = colorRamp2(c(0,1), c("white", "black"))


titl=prot
if (prot=='AF9_Freiburg_ORF'){titl='AF9'}
if (prot=='ELOA2'){titl='Elongin-A2'}
if (prot=='huMeCP2'){titl='MeCP2'}

ha=HeatmapAnnotation(charge=anno_barplot(ma,gp=gpar(col=ifelse(ma>0, "red", "blue")),size=unit(2,"points")),height=unit(4, "cm")); 
ht=Heatmap(column_title=titl, t(r),cluster_rows=FALSE,cluster_columns=FALSE,show_heatmap_legend=FALSE,
	bottom_annotation = ha, 
	#right_annotation=hs, 
	col=col_fun, border =1);

pdf(paste('single.',prot,'.complot.pdf',sep=""),width=20)
draw(ht)
dev.off()

}