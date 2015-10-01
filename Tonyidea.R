#settings
datadir <- '/home/heather/Data/RAD/Quercus'
setwd('/home/heather/Projects/RAD/Quercus/pairwise')
system(paste('ls',datadir,'> temp.txt'))
temp <- read.table("temp.txt", quote="\"")
#temp$command1<-paste("gunzip /home/heather/Data/Mosquito/2010/notag/",temp$V1,"/*.fasta.gz",sep="")
#write.table(temp$command1,'gunzip.sh',row.names=F,col.names=F,quote=F)
k=21 #k used for reads selection
Quercus_CB <- read.delim("Quercus_Cavender-Bares_committeemeetingversion.txt", header=F)
# #Quercus_CB$kmercountcommand<-paste("kmer_count -l ",k," -G 6 -o /home/heather/Data/Mosquito/2010/notag/",
#                                    Quercus_CB$V1,".pkdat.gz -f FA -i /home/heather/Data/Mosquito/2010/notag/",
#                                    Quercus_CB$V1,'/',Quercus_CB$V2,'_notag.fasta > ', Quercus_CB$V1,'_k',
#                                      k,".wc &",sep='')
#write.table(Quercus_CB$kmercountcommand,"kmer_count.sh",row.names=F,col.names=F,quote=F)
#Construct phylogeny before selection
phylokmer1.call<-paste('python ~/bin/aaf_phylokmer.py -k', k, '-t 10 -G 10 -d',datadir)
phylokmer1.call
system(phylokmer1.call)

combm<-combn(temp$V1,2)
pk=21 #k used for distance calculation
kmer_merge=NULL
ReadsSelector=NULL
mkdir=NULL
phylokmer=NULL
distance=NULL
for (col in 1:ncol(combm)) {
  kmer_file<-paste(combm[1,col],"__",combm[2,col],"_k",k,".kmer",sep='')
  kmer_merge<-c(kmer_merge,(paste("kmer_merge -k s -c -d '0' -A A -B A ",
                                  combm[1,col],".pkdat.gz ",combm[2,col],
                                  ".pkdat.gz | cut -f 1 > ",kmer_file,sep='')))
  mkdir<-c(mkdir,paste("mkdir ",combm[1,col],combm[2,col],sep=''))
  mkdir<-c(mkdir,paste("mkdir ",combm[1,col],combm[2,col],"/",combm[1,col],sep=''))
  mkdir<-c(mkdir,paste("mkdir ",combm[1,col],combm[2,col],"/",combm[2,col],sep=''))
  ReadsSelector<-c(ReadsSelector,(paste("ReadsSelector -k ",kmer_file, " -o ",combm[1,col],combm[2,col],"/",combm[1,col],
                                        "/",combm[1,col],combm[2,col],"_",combm[1,col],
                                        "_k",k,".fa", " -s ",datadir,'/',combm[1,col],
                                        "/SRR*.fa",sep='')))
  ReadsSelector<-c(ReadsSelector,(paste("ReadsSelector -k ",kmer_file, " -o ",combm[1,col],combm[2,col],"/",combm[2,col],
                                        "/",combm[1,col],combm[2,col],"_",combm[2,col],
                                        "_k",k,".fa", " -s ",datadir,'/',combm[2,col],
                                        "/SRR*.fa",sep='')))
  phylokmer<-c(phylokmer,paste("python ~/bin/aaf_phylokmer.py -k ",pk," -t 2 -n 2 -f FA -d ",combm[1,col],combm[2,col]," -G 100",sep=''))
  distance<-c(distance,paste("python ~/bin/aaf_distance.py -i ",combm[1,col],combm[2,col],"/phylokmer.dat.gz -t 15 -G 5 -o ",
                             combm[1,col],"__",combm[2,col],"_k",k,"_pk",pk," -f ",combm[1,col],combm[2,col],"/kmer_diversity.wc",sep=''))
}
write.table(kmer_merge,"kmer_merge.sh",row.names=F,col.names=F,quote=F)
write.table(mkdir,"mkdir.sh",row.names=F,col.names=F,quote=F)
write.table(ReadsSelector,"ReadsSelector.sh",row.names=F,col.names=F,quote=F)
write.table(phylokmer,"phylokmer.sh",row.names=F,col.names=F,quote=F)
write.table(distance,"distance.sh",row.names=F,col.names=F,quote=F)

system('sh kmer_merge.sh')
system('sh mkdir.sh')
system('sh ReadsSelector.sh')
system('sh phylokmer.sh')
system('sh distance.sh')


