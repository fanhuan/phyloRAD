#One simulation
#dependencies
#source("http://bioconductor.org/biocLite.R")
#biocLite("Biostrings")


if(!require("readr")){
  install.packages("readr")
}
library(readr)
library(phylosim)
library(ape)
library(Biostrings)
#rank_fh is a function that takes a dataframe and add a rank column according to one column (given the column number)
#such that the numbers above 0 are ranked from high to low as 1,2,3,4 while numbers below 0 are ranked from
#low to high as -1, -2, -3, -4. 0 is ranked as 0.
rank_fh <- function(df,number){
  df_ord<-df[order(df[,number],decreasing = TRUE),]
  i <- 1
  for (a in df_ord[,number]) {
    if (a > 0) {
      df_ord$rank[i]<-i
    }
    else if (a == 0) {
      df_ord$rank[i]<-0
    }
    else {
      df_ord$rank[i]<-i-nrow(df_ord)-1
    }
    i <- i+1
  }
  return(df_ord)
}

setwd("/home/heather/R/GWAS_TB")
#name this simulation by assigning repn
repn=091101
#Enable the "fast and careless mode"
PSIM_FAST <- TRUE

#codon freqs for flanking regions are calculated based on 11 rpoB genes in Mycobacterium genus, stored in rpoB_alignment.fa 
codon.freqs.flanking <- c(0.0017830839599969, 0.0276765640747345, 7.75253895650826e-05, 0.0100007752538957, 
                          0.00217071090782231, 0.017830839599969, 0.00186060934956198, 0.0233351422590899, 0.00139545701217149, 
                          0.0186836188851849, 0.00139545701217149, 0.00565935343825103, 0.00767501356694317, 0.0031010155826033, 
                          0.0128692146678037, 0.000930304674780991, 0.0560508566555547, 0.00248081246608264, 0.013101790836499, 
                          0.00240328707651756, 0.0368245600434142, 0.0013179316226064, 0.0174432126521436, 0.00240328707651756, 
                          0.0269788355686487, 0.0106209783704163, 0.0382975424451508, 0.00248081246608264, 0.0179083649895341, 
                          0.00449647259477479, 0.0512442825025196, 0.000387626947825413, 0.0259710055043027, 0.0031010155826033, 
                          0.0369796108225444, 0.00155050779130165, 0.0133343670051942, 0.00186060934956198, 0.0323280874486394, 
                          0.00317854097216838, 0.0396929994573223, 0.0035661679199938, 0.00922552135824482, 0, 0.000697728506085743, 
                          0.00627955655477169, 0.0421738119234049, 0.000930304674780991, 0.0434142181564462, 0.00465152337390495, 
                          0.0372121869912396, 0.00348864253042872, 0.0255058531669122, 0.00814016590433367, 0.062020311652066, 
                          0.0146522986278006, 0.0747344755407396, 0.0240328707651756, 0.0464377083494845, 0.0035661679199938, 
                          0.0103884022017211)

p.flanking <- GY94(codon.freqs = codon.freqs.flanking)
p.flanking$kappa = 2 #kappa is the transition/transversion rate
summary(p.flanking)

##set up hotspot
#6 codons for serine
serine.codons <- c("TCT", "TCC", "TCA", "TCG", "AGT", "AGC")

# from Kumar&Jena: The single amino acid mutation at codon 450 (i.e., Ser to Leu [TCG to TTG]) of the rpoB gene is reported to be the most widespread mutation 
# setting up the code frequencies this way means that transitions from TCG (or any other codon) to TTG is high. Since the ancestral codon is TCG, this sets the rate of loss of this codon.
codon.freqs.S450 <- array(1,dim=c(1,61))# no stop codon
codon.freqs.S450[4] <- 20 #TTG, twenty times higher, selected
codon.freqs.S450 <- codon.freqs.S450/sum(codon.freqs.S450)
p.S450 <- GY94(codon.freqs = codon.freqs.S450, rate.multiplier = 10)
p.S450$kappa = 2

#construct the root sequence (M.canettii)
string="TTGGCAGATTCCCGCCAGAGCAAAACAGCCGCTAGTCCTAGTCCGAGTCGCCCGCAAAGTTCCTCGAATAACTCCGTACCCGGAGCGCCAAACCGGGTCTCCTTCGCTAAGCTGCGCGAACCACTTGAGGTTCCGGGACTCCTTGACGTCCAGACCGATTCGTTCGAGTGGCTGATCGGTTCGCCGCGCTGGCGCGAATCCGCCGCCGAGCGGGGTGATGTCAACCCAGTGGGTGGCCTGGAAGAGGTGCTCTACGAGCTGTCTCCGATCGAGGACTTCTCCGGGTCGATGTCGTTGTCGTTCTCTGACCCTCGTTTCGACGATGTCAAGGCACCCGTCGACGAGTGCAAAGACAAGGACATGACGTACGCGGCTCCACTGTTCGTCACCGCCGAGTTCATCAACAACAACACCGGTGAGATCAAGAGTCAGACGGTGTTCATGGGTGACTTCCCGATGATGACCGAGAAGGGCACGTTCATCATCAACGGGACCGAGCGTGTGGTGGTCAGCCAGCTGGTGCGGTCGCCCGGGGTGTACTTCGACGAGACCATTGACAAGTCCACCGACAAGACGCTGCACAGCGTCAAGGTGATCCCGAGCCGCGGCGCGTGGCTCGAGTTTGACGTCGACAAGCGCGACACCGTCGGCGTGCGCATCGACCGCAAACGCCGGCAACCGGTCACCGTGCTGCTCAAGGCGCTGGGCTGGACCAGCGAGCAGATTGTCGAGCGGTTCGGGTTCTCCGAGATCATGCGATCGACGCTGGAGAAGGACAACACCGTCGGCACCGACGAGGCGCTGTTGGACATCTACCGCAAGCTGCGTCCGGGCGAGCCCCCGACCAAAGAGTCAGCGCAGACGCTGTTGGAAAACTTGTTCTTCAAGGAGAAGCGCTACGACCTGGCCCGCGTCGGTCGCTATAAGGTCAACAAGAAGCTCGGGCTGCATGTCGGCGAGCCCATCACGTCGTCGACGCTGACCGAAGAAGACGTCGTGGCCACCATCGAATATCTGGTCCGCTTGCACGAGGGTCAGACCACGATGACCGTTCCGGGCGGCGTCGAGGTGCCGGTGGAAACCGACGACATCGACCACTTCGGCAACCGCCGCCTGCGTACGGTCGGCGAGCTGATCCAAAACCAGATCCGGGTCGGCATGTCGCGGATGGAGCGGGTGGTCCGGGAGCGGATGACCACCCAGGACGTGGAGGCGATCACACCGCAGACGTTGATCAACATCCGGCCGGTGGTCGCCGCGATCAAGGAGTTCTTCGGCACCAGCCAGCTGAGCCAATTCATGGACCAGAACAACCCGCTGTCGGGGTTGACCCACAAGCGCCGACTGTCGGCGCTGGGGCCCGGCGGTCTGTCACGTGAGCGTGCCGGGCTGGAGGTCCGCGACGTGCACCCGTCGCACTACGGCCGGATGTGCCCGATCGAAACCCCTGAGGGGCCCAACATCGGTCTGATCGGCTCGCTGTCGGTGTACGCGCGGGTCAACCCGTTCGGGTTCATCGAAACGCCGTACCGCAAGGTGGTCGACGGCGTGGTTAGCGACGAGATCGTGTACCTGACCGCCGACGAGGAGGACCGCCACGTGGTGGCACAGGCCAATTCGCCGATCGATGCGGACGGTCGCTTCGTCGAGCCGCGCGTGCTGGTCCGCCGCAAGGCGGGCGAGGTGGAGTACGTGCCCTCGTCTGAGGTGGACTACATGGACGTCTCGCCCCGCCAGATGGTGTCGGTGGCCACCGCGATGATTCCCTTCCTGGAGCACGACGACGCCAACCGTGCCCTCATGGGGGCAAACATGCAGCGCCAGGCGGTGCCGCTGGTCCGTAGCGAGGCCCCGCTGGTGGGCACCGGGATGGAGCTGCGCGCGGCGATCGACGCCGGCGACGTCGTCGTCGCCGAAGAAAGCGGCGTCATCGAGGAGGTGTCGGCCGACTACATCACTGTGATGCACGACAACGGCACCCGGCGTACCTACCGGATGCGCAAGTTTGCCCGGTCCAACCACGGCACTTGCGCCAACCAGTGCCCCATCGTGGACGCGGGCGACCGAGTCGAGGCCGGTCAGGTGATCGCCGACGGTCCCTGTACTGACGACGGCGAGATGGCGCTGGGCAAGAACCTGCTGGTGGCCATCATGCCGTGGGAGGGCCACAACTACGAGGACGCGATCATCCTGTCCAACCGCCTGGTCGAAGAGGACGTGCTCACCTCGATCCACATCGAGGAGCATGAGATCGATGCTCGCGACACCAAGCTGGGTGCGGAGGAGATCACCCGCGACATCCCGAACATCTCCGACGAGGTGCTCGCCGACCTGGATGAGCGGGGCATCGTGCGCATCGGTGCCGAGGTTCGCGACGGGGACATCCTGGTCGGCAAGGTCACCCCGAAGGGTGAGACCGAGCTGACGCCGGAGGAGCGGCTGCTGCGTGCCATCTTCGGTGAGAAGGCCCGCGAGGTGCGCGACACTTCGCTGAAGGTGCCGCACGGCGAATCCGGCAAGGTGATCGGCATTCGGGTGTTTTCCCGCGAGGACGAGGACGAGTTGCCGGCCGGTGTCAACGAGCTGGTGCGTGTGTATGTGGCTCAGAAACGCAAGATCTCCGACGGTGACAAGCTGGCCGGCCGGCACGGCAACAAGGGCGTGATCGGCAAGATCCTGCCGGTTGAGGACATGCCGTTCCTTGCCGACGGCACCCCGGTGGACATTATTTTGAACACCCACGGCGTGCCGCGACGGATGAACATCGGCCAGATTTTGGAGACCCACCTGGGTTGGTGTGCCCACAGCGGCTGGAAGGTCGACGCCGCCAAGGGGGTTCCGGACTGGGCCGCCAGGCTGCCCGACGAACTGCTCGAGGCGCAGCCGAACGCCATTGTGTCGACGCCGGTGTTCGACGGCGCCCAGGAGGCCGAGCTGCAGGGCCTGTTGTCGTGCACGCTGCCCAACCGCGACGGTGACGTGCTGGTCGACGCCGACGGCAAGGCCATGCTCTTCGACGGGCGCAGCGGCGAGCCGTTCCCGTACCCGGTCACGGTTGGCTACATGTACATCATGAAGCTGCACCACCTGGTGGACGACAAGATCCACGCCCGCTCCACCGGGCCGTACTCGATGATCACCCAGCAGCCGCTGGGCGGTAAGGCGCAGTTCGGTGGCCAGCGGTTCGGGGAGATGGAGTGCTGGGCCATGCAGGCCTACGGTGCTGCCTACACCCTGCAGGAGCTGTTGACCATCAAGTCCGATGACACCGTCGGCCGCGTCAAGGTGTACGAGGCGATCGTCAAGGGTGAGAACATCCCGGAGCCGGGCATCCCCGAGTCGTTCAAGGTGCTGCTCAAAGAACTGCAGTCGCTGTGCCTCAACGTCGAGGTGCTATCGAGTGACGGTGCGGCGATCGAACTGCGCGAAGGTGAGGACGAGGACCTGGAGCGGGCCGCGGCCAACCTGGGAATCAATCTGTCCCGCAACGAATCCGCAAGTGTCGAGGATCTTGCG"
rpoB <- CodonSequence(string=string)
substr(rpoB,1348,1350) #Check the S450 

attachProcess(rpoB, p.flanking, 1:449) #left hand side
attachProcess(rpoB, p.S450, 450) #susceptive spot
attachProcess(rpoB, p.flanking, 451:1172) #right hand side

#Construct a deletion process proposing deletions with rate 0.00083333 according to a discrete length distribution
d <- DiscreteDeletor(rate = 0.00083333, sizes = c(1, 2, 3, 4, 5, 6, 7), probs = c(0.2, 0.2, 0.2, 
                                                                                  0.1, 0.1, 0.1, 0.1))

#Construct an insertion process proposing insertions with rate 0.25 according to a discrete length distribution
i <- DiscreteInsertor(rate = 0.00083333, sizes = c(1, 2, 3, 4, 5, 6, 7), probs = c(0.2, 0.2, 0.2, 
                                                                                   0.1, 0.1, 0.1, 0.1))

#Set the templete sequence for the insertion process (what to insert)
i$templateSeq <- NucleotideSequence(length = 7, processes = list(list(p.flanking)))

#Attaching the indel processes to the flanking regions.
attachProcess(rpoB,d,1:449)
attachProcess(rpoB,d,451:1172)
attachProcess(rpoB,i,1:449)
attachProcess(rpoB,i,451:1172)

# Sample omegas from a discrete model. Omegas are dN/dS. 1 is neutral. >1 if favored, <1 if deleterious
# Omegas are only useful to measure a general selection process on codon, not specific selection on specific codon. I have just set the omegas to 1, although this could be changed.
omegaVarM0(rpoB, p.flanking, omega = 1, index = c(1:449, 
                                                  451:1172))
omegaVarM0(rpoB, p.S450, omega = 1, index = 450)

setRateMultipliers(rpoB, p.flanking, value = .2, index = c(1:449, 
                                                           451:1172))
setRateMultipliers(rpoB, p.S450, value = .2, index = 450)

# Input tree
cat("(((((MAF_11821_03:0.412,(MAF_GM_0981:0.582,MTB_95_0545:0.550):0.184):0.212,((MTB_K21:0.242,(MTB_K67:0.280,MTB_K93:0.234):0.096):0.108,(MTB_T17:0.624,MTB_T92:0.584):0.150):0.114):0.096,(((((MTB_00_1695:0.510,MTB_T67:0.518):0.182,MTB_T85:0.278):0.154,(MTB_98_1833:0.392,MTB_M4100A:0.072):0.142):0.134,(MTB_4783_04:0.574,MTB_GM_1503:0.596):0.166):0.118,MTB_91_0079:0.576):0.094):0.074,(MTB_K37:0.518,MTB_K49:0.288):0.070):0.138,MTB_H37Rv:0.256);",file="TBSimulation.nwk")
phy.base <- read.tree(file = "TBSimulation.nwk")
phy <- phy.base
p <- Ntip(phy) #number of tips/species 

# This can be used to rescale the branch lengths of the tree
# phy$edge.length <- phy$edge.length

# check with a single example
system.time(sim <- Simulate(PhyloSim(root.seq = rpoB, phylo = phy), quiet = T))
saveRDS(sim, file=paste("sim_",repn,".rds",sep=''))
alignment.names <- names(sim$alignment[1, ])

# This identifies the location (index) of 450
index <- (1:1172)[alignment.names == 450 & !is.na(alignment.names)]
#plot(sim, aln.xlim=index + c(-13.5, 13.5), num.pages=1)

# This makes a data.frame for 450 containing codons and whether they code for serine 
w <- data.frame(tip = phy$tip.label, codon = NA, serine = T)
selected.codons <- sim$alignment[,index]
for (kk in phy$tip.label) {
  w$codon[w$tip == kk] <- selected.codons[names(selected.codons) == kk]
  w$serine[w$tip == kk] <- is.element(selected.codons[names(selected.codons) == kk], serine.codons)
}

write.csv(w, file = paste("rpoB_GY84_trait_", repn, ".csv", sep = ""))

# Function for simulating and saving alignments. It also saves as a separate file 
#the phenotype (serine or not at 450) and the codons at 450. Only simulations with
#at least 3 serines and 3 non-serines are saved.

if((sum(w$serine) >= 3) & (sum(w$serine) <= length(w$serine)-3)){
  saveAlignment(sim, file = paste("rpoB_GY84_", repn, ".fas", sep = ""), skip.internal = T) 
} else {
  print("This simulation ended up with two many or too few seriens, abort")
  stop()
}

#Run aaf_phylosim.py and kmer_pattern.py
#TBSimulation Huan$ aaf_phylosim.py -k 9 -o phylokmer.dat -i rpoB_GY84_102.fas -W
k<-9
phylokmer.call <- paste("python ~/bin/aaf_phylosim.py -k ",as.character(k)," -o phylokmer.dat -i rpoB_GY84_",
                        as.character(repn), ".fas -W", sep="")
system(phylokmer.call, ignore.stdout = F)
#TBSimulation Huan$ python kmer_pattern.py -i rpoB_GY84_102/phylokmer.dat -o rpoB_GY84_102_
kmerpatter.call <- paste("python ~/bin/kmer_pattern.py -i rpoB_GY84_", as.character(repn), 
                         "/phylokmer.dat -p 1 -o rpoB_GY84_", as.character(repn),"_", sep="")
system(kmerpatter.call, ignore.stdout = F)
#Load in the shared kmer table
phylokmer <- read.delim(paste("rpoB_GY84_", repn, "/phylokmer.dat", sep = ""), header=FALSE)
colnames(phylokmer)<-c("kmer",sort(phy$tip.label))
#Grep the patterns for kmers invloving S450 (kmer_list)

#generate kmers including S450 (we do need to run aaf_phylosim.py first)
kmer<-sim$alignment[,(index-k/3+1):(index+k/3-1)]
kmer_list<-NULL #This does not consider possible deletions ('NA' in some codon)
for (x in 1:p) {
  kmers<-paste(kmer[x,],collapse='')
  for (j in 1:(nchar(kmers)-k+1)){
    kmer_list<-c(kmer_list, substr(kmers, j, j+k-1))
  } 
}
kmer_df<-data.frame(kmer_list)
kmer_df<-rename(kmer_df,c("kmer_list"="kmer"))
kmer_df$position<-rep(c((nchar(kmers)-k+1):1),p)
#one kmer won't have different positions, too short
kmer_df<-unique(kmer_df)
#add their reverse compliment conterparts in (kmer_count only uses whoever that is alphabetically
#first, the original kmer or it's rc.)
rc<-NULL
for (i in 1:nrow(kmer_df)) {
  kmer_df<-rbind(kmer_df,data.frame(kmer=as.character(reverseComplement(DNAString(kmer_df$kmer[i]))), 
                                    position=-kmer_df$position[i]))
}
kmer_df<-unique(kmer_df)
kmer_list_rc<-kmer_df$kmer
S450_kmers<-phylokmer[phylokmer$kmer %in% kmer_list_rc,]
#Note that less than half of the kmers in kmer_list_rc ends up in S450_kmers because
#1. only the original OR the rc kmer is in phylokmer
#2. only kmers shared at least by two species are in phylokmer

#score ONCE. We need information from simulations, not only the fas output.
C <- as.matrix(vcv(phy))
C <- C/det(C)^(1/p) #some transformation of C to make the determinant = 1, which should make calculations easier
#n <- nr.replicates
n<-1
iC <- solve(C) #inverse of C
ones <- array(1, c(p, 1)) # array of ones for the intercept

#output <- array(0, c(n, 1)) # collects the results for the n patterns
#threshold <- 3 # minimum number of zeros or ones for inclusion of pattern #already taken care of in kmer_pattern.py
#y should be kmer pattern
#X should be the trait
#read in the kmer pattern as Y, and calculate score for each y.
Y <- read.fwf(file=paste("rpoB_GY84_", repn, "_kmerPattern.stats", sep = ""), widths=array(1,c(1,p)),
              header=F)
colnames(Y)<-sort(phy$tip.label)

#sort the dataframe by tip name, get the serine column with Trues and falses and convert it into 0/1(by *1 or +0)
trait<-w[order(w$tip),]$serine
X<-t(t(trait*1))
threshold<-3
#set up dataframe for scores
score<-array(0, c(nrow(Y), 1))
pattern<-array(NA,c(nrow(Y),1))
output<-data.frame(pattern,score)
for (i in 1:nrow(Y)) {
  y<-Y[i,]
  output$pattern[i]<-paste(y,collapse='')
  sumy <- sum(y)#sum of the pattern
  if (sumy < threshold | sumy > (p - threshold)) {
    output$score[i] <- 0 #all patterns with only 1 one or zero are assigned score 0.
  } else {
    xx <- cbind(ones, X)
    
    XiCX <- t(xx) %*% iC %*% xx
    XiCY <- t(xx) %*% iC %*% t(y)
    
    b <- solve(XiCX, XiCY)
    
    h <- t(y) - (xx %*% b)
    
    MSE <- t(h) %*% iC %*% h/(p - 2)
    
    iXiCX <- solve(XiCX)
    bSE <- (MSE * iXiCX[2, 2])^0.5
    output$score[i] <- b[2]/bSE
  }
}
output_ranked<-rank_fh(output,2)
#read in the kmer pattern again, focus on the frequency this time.
kmerPattern.stats<-read.table(file=paste("rpoB_GY84_", repn, "_kmerPattern.stats", sep = ""),
                              colClasses=c('character', 'integer'),
                              col.names=c('pattern','freq'))
total_score<-merge(output_ranked,kmerPattern.stats,by='pattern',all.x=T)
total_score<-total_score[order(total_score$score,decreasing=T),]
write.csv(total_score,paste("rpoB_GY84_", as.character(repn), "_scores.csv", sep=""),row.names = FALSE, )

#Get scores for S450 containing kmers.
kft<-1 #kmer frequency threshold
pattern450<-array(0,nrow(S450_kmers))
for (i in 1:nrow(S450_kmers)){
  pattern450_j<-array(1,ncol(S450_kmers)-1)
  for (j in 2:ncol(S450_kmers)){
    if (S450_kmers[i,][j]<kft){
      pattern450_j[j-1]<-0
    }
  }
  pattern450[i]<-paste(pattern450_j,collapse='' )
}
S450_kmers$pattern<-pattern450
S450_kmers_score<-merge(S450_kmers,total_score,by="pattern",all.x=T)
#Note that some of the kmers does not have scores because only patterns with more than two 1 or 0 are scored.
#Merge S450_kmer_score with kmer_df
output_450<-merge(S450_kmers_score, kmer_df, by="kmer",all.x=T)
output_450_light<-subset(output_450,select=c("kmer","position","pattern",'freq',"score",'rank'))
write.csv(output_450_light,paste("rpoB_GY84_", as.character(repn), "_450summary.csv", sep=""),row.names = FALSE, )
#output_450<-output[output$pattern %in% pattern450,]