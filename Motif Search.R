library(anopheles.db0)
library(org.Ag.eg.db)
library(AnnotationDbi)
library(Biostrings)
library(MotifDb)
library(GenomicFeatures)
library(motifStack)
library(seqLogo)
library(MotifDb)
library(MotIV)
library(seqinr)
library("Biostrings")



seqs<-readDNAStringSet('/Users/vickyingham/Downloads/mart_export.txt')
seqs.null<-readDNAStringSet('/Users/vickyingham/Documents/PhD Year 3/Array/Maf-S 2 - GFP/5_arrays/non-significant up stream 100bp.txt')


numbers.null<-sample(1:length(seqs.null), length(seqs))

seqs.2<-seqs.null[numbers.null]

#Transcription factor binding sites of interest
C<-(c(113,20,0,0,46,292,0,215,33,0,474,0,128,63,12)/474)
T<-(c(60,0,474,0,0,102,225,91,142,28,0,0,114,190,248)/474)
G<-(c(159,114,0,474,21,80,169,90,163,422,0,106,99,0,0)/474)
A<-c(142/474,340/474,0,0,407/474,0,80/474,78/474,136/474,24/474,0,368/474,133/474,221/474,214/474)
mafs<-rbind(A,C,G,T)
sequence.for.flank.search<-mafs

A<-(c(726,0,0,1090,16,110,80,978,8,0,931,383,361,292,211)/1090)
C<-(c(26,0,0,0,871,1,937,0,4,1063,17,203,100,102,283)/1090)
G<-(c(334,0,1082,0,166,12,22,51,1041,22,73,283,130,91,106)/1090)
T<-(c(4,1090,8,0,37,967,51,61,37,5,69,221,499,605,490)/1090)
mafs.homo<-rbind(A,C,G,T)
sequence.for.flank.search<-mafs.homo

A<-c(0,0,1,0,.25,.25,.25,0,0)
C<-c(0,0,0,0.5,0.25,0.25,0.25,0,1)
G<-c(0,1,0,0,0.25,0.25,0.25,1,0)
T<-c(1,0,0,0.5,0.25,0.25,0.25,0,0)
ARE<-rbind(A,C,G,T)
sequence.for.flank.search<-ARE

A<-c(0,0.5,1,0.25,0.25,0.5,0,0,1,0,0.25,0.25,0.25,0,0,0.5,0.5,0.5,0.5,0.5)
C<-c(0,0.5,0,0.25,0.25,0,0,0,0,0.5,0.25,0.25,0.25,0,1,0,0,0,0,0)
G<-c(0,0,0,0.25,0.25,0.5,0,1,0,0,0.25,0.25,0.25,1,0,0.5,0,0,0,0)
T<-c(1,0,0,0.25,0.25,0,1,0,0,0.5,0.25,0.25,0.25,0,0,0,0.5,0.5,0.5,0.5)
ARE.full<-rbind(A,C,G,T)
sequence.for.flank.search<-ARE.full

C<-c(1,0,1,0,.5,0,0,0,0)
T<-c(0,0,0,0,.5,0,0,.5,0)
G<-c(0,0,0,1,0,1,.5,0,1)
A<-c(0,1,0,0,0,0,0.5,0.5,0)
met.aedes<-rbind(A,C,G,T)
sequence.for.flank.search<-met.aedes

C<-c(1,0,1,0,0,0)
T<-c(0,0,0,0,1,0)
G<-c(0,0,0,1,0,1)
A<-c(0,1,0,0,0,0)
met.CACGTG<-rbind(A,C,G,T)
sequence.for.flank.search<-met.CACGTG

C<-c(1,0,0,0,0,0)
T<-c(0,0,0,0,1,1)
G<-c(0,1,1,0,0,0)
A<-c(0,0,0,1,0,0)
met.CGGATT<-rbind(A,C,G,T)
sequence.for.flank.search<-met.CGGATT

C<-c(1,0,1,0,1,0)
T<-c(0,0,0,0,0,0)
G<-c(0,0,0,1,0,1)
A<-c(0,1,0,0,0,0)
met.CACGCG<-rbind(A,C,G,T)
sequence.for.flank.search<-met.CACGCG

par(mfrow=c(2,2))
mdb<-MotifDb
query.result<-query(query(MotifDb, 'Dmelanogaster'),"met") ##Change to the transcription factor of interest
print(query.result)

seqLogo(met.aedes) ##this depends on the number of results from the print query.results, change accordingly
seqLogo(query.result[[3]])


sequence.for.flank.search2<-query.result[[1]] ##Again change according the number of different putative binding sites
sequence.for.flank.search<-query.result[[1]] 

#How many times does motif appear in each sequence?
vector.out<-c()
for(i in 1:length(seqs))
{
  matched<-matchPWM(sequence.for.flank.search,unname(substring(seqs[i],1,1000)))
  if(length(matched)>0)
  {
    matched2<-toString(matched)
    strings<-c()
    for(j in 1:length(matched2))
    {
      strings<-c(strings,matched2[j])
    }
    vector.out<-rbind(vector.out,c(substr(names(seqs[i]),12,nchar(names(seqs[i]))),length(matched),strings))
  }
}
sum(as.numeric(vector.out[,2]))
write.table(vector.out,'/Users/vickyingham/Documents/PhD Year 3/Array/Maf-S 2 - GFP/5_arrays/Actual significant probes/1000bp cnc overexpressed ARE-full binding sites.txt',row.names=F,sep='\t')

#Iteratively search flank for selected motif. To change to null change all occurences of seq to seq.2 in the j in: for loop below. If only searching for none null, use only inner for loop
vector.out.2<-c()
 for(k in 1:10)
 {
   numbers.null<-sample(1:8000, length(seqs))
   seqs.2<-seqs.null[numbers.null]
for(j in 1:(1000-ncol(sequence.for.flank.search)))
{
  starts<-seq(1,1000-ncol(sequence.for.flank.search),by=1)
  ends<-seq(ncol(sequence.for.flank.search),1000,by=1)
  vector.out<-c()
  sections<-subseq(seqs,start=starts[j],end=ends[j])
  for(i in 1:length(sections))
  {
    sequence.length<-as.character(sections[i[1]])
    flank.vector<-paste(sequence.length,sep="",collapse="")
    matched<-matchPWM(sequence.for.flank.search,flank.vector)
    if(length(matched)>=1)
    {
      vector.out<-rbind(vector.out,length(matched))
      print(j)
    }
    else
    {
      vector.out<-rbind(vector.out,0)
    }
  }
  vector.out.2<-c(vector.out.2,sum(vector.out))
  }
}
vector.out.2


#Plot the results for iterative search

number.of.motifs<-read.delim('/Users/vickyingham/Desktop/CACGCG signficant.txt',header=T)
number.of.motifs<-as.matrix(number.of.motifs)
number.of.motifs<-t(number.of.motifs)
number.of.motifs<-as.vector(number.of.motifs)

null.motifs<-read.delim('/Users/vickyingham/Desktop/CACGTG Non Sig.txt',header=T)
null.motifs<-as.matrix(null.motifs)
null.motifs<-t(null.motifs)
null.motifs<-as.vector(null.motifs)

null.motifs.2<-cbind(null.motifs[1:995],null.motifs[996:1991])

vector.high<-c()
vector.low<-c()
for(i in 1:nrow(null.motifs.2))
{
  numbers<-as.numeric(null.motifs.2[i,])
  vector.high<-rbind(vector.high,(mean(numbers)))
}

plot((number.of.motifs),type='l',col='skyblue',ylim=c(min(number.of.motifs),max(number.of.motifs)),pch=15,xlab='Nucleotide Position',ylab='Number of Motifs')
lines((vector.high),type='l',pch=22,lty=2,col='black')
lines((vector.low),type='l',pch=22,lty=2,col='black')



plot(met1[1,],type="o",col='gray',ylim=c(min(met1),max(met1)),pch=22,lty=2,axes=F,ann=F)
axis(1,at=1:10,lab=c('1-100','101-200','201-300','301-400','401-500','501-600','601-700','701-800','801-900','901-1000'),lwd=2)
axis(2,las=1,at=2*0:10,lwd=2)
lines(met1[3,],type="o",col='gray',pch=22,lty=2)
lines(met1[2,],type="o",col='gray',pch=15)

#sequence.for.flank<-read.delim('/Users/vickyingham/Documents/PhD Year 2/TF information/Flanking Region Search/Upstream/A Priori/A priori upstream.txt',header=F) #File from http://genome.ucsc.edu/cgi-bin/hgTables
#sequence.for.flank<-as.matrix(sequence.for.flank)


#Search for transcription factor using databases
par(mfrow=c(2,2))
query.result<-query(query(MotifDb, 'Dmelanogaster'),"met") ##Change to the transcription factor of interest
print(query.result)

seqLogo(query.result[[1]]) ##this depends on the number of results from the print query.results, change accordingly
seqLogo(query.result[[3]])


sequence.for.flank.search<-query.result[[1]] ##Again change according the number of different putative binding sites
sequence.for.flank.search.2<-query.result[[3]] 

vector.out.2<-c()

positions<-grep('>',sequence.for.flank)

for(i in 1:length(positions))
{
  flank.vector<-paste(sequence.for.flank[(positions[i]+1):(positions[i]+5),],sep="",collapse="")
  matched<-matchPWM(sequence.for.flank.search,flank.vector) ## Change again for the number of binding sites
  matched2<-matchPWM(sequence.for.flank.search.2,flank.vector)
  if(length(matched)>=1)
  {
    vector.out<-cbind(sequence.for.flank[positions[i],1],as.character(matched))
    vector.out.2<-rbind(vector.out.2,vector.out)
  }
 if(length(matched2)>=1)
{
  vector.out<-cbind(sequence.for.flank[positions[i],1],as.character(matched2))
  vector.out.2<-rbind(vector.out.2,vector.out)
}
}

description.vector<-c()
descriptions<-read.delim('/Users/vickyingham/Documents/PhD Year 1/GaGa /Ortholog_Searches/anopheles_description_file.txt',header=T)
descriptions<-as.matrix(descriptions)

for(i in 1:nrow(vector.out.2))
{
  substring<-substr(vector.out.2[i,1],2,11)
  position<-grep(substring,descriptions[,2])
  out.vector<-cbind(substring,as.character(descriptions[position[1],5]))
  description.vector<-rbind(description.vector,out.vector)
}

vector.out.3<-cbind(description.vector,vector.out.2[,2])
  
write.table(vector.out.3,'/Users/vickyingham/Documents/PhD Year 2/TF information/Flanking Region Search/Upstream/A Priori/a priori cyc.txt',row.names=F,sep='\t')


