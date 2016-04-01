#Simulation 2
setwd("/home/heather/R/RAD/simulation")
Sys.setenv(PATH = "/sbin:/usr/sbin:/bin:/usr/bin:/home/heather/bin")
rep<-100
infile<-'/home/heather/R/RAD/shallow30M'
enzyme <- 'SbfI'
rate<-c(0,0.05,0.1)
kc<-c(13,15,17)
intree<-'/home/heather/R/RAD/shallow.tre'
system(paste('cp',intree,'intree'))
mistake<-2    #topological mistakes with AAF on alignment
results<-data.frame(dropRate=numeric(0),rep=integer(0),k=integer(0),
                    percentb=numeric(0), #percentage of SBR before selection
                    percenta=numeric(0), #percentage of SBR after selection
                    mistakeb=numeric(0), #mistakes before selection
                    mistakea=numeric(0), #mistakes after selection
                    mistakesb=numeric(0)) #mistakes made by SBR only

# prepare setting files for treedist, treedistsetting contains everything but intree2
system(paste('cp /home/heather/R/RAD/treedistsetting treedistsetting_test'))
echo.call <- 'echo test.tre1 | cat >> treedistsetting_test'
system(echo.call)
system(paste('cp /home/heather/R/RAD/treedistsetting treedistsetting_all'))
echo.call <- 'echo test_all.tre1 | cat >> treedistsetting_all'
system(echo.call)
system(paste('cp /home/heather/R/RAD/treedistsetting treedistsetting_selected'))
echo.call <- 'echo test_selected.tre1 | cat >> treedistsetting_selected'
system(echo.call)

# this is one for one rose simulation, god! I think we do need condor!
row = 0
for (r in rate) {
  for (i in (1:rep)) {
    system('rm test -f -r')
    system('mkdir test')
    simRad.call<-paste('python ~/bin/simRadAlignment.py -i',infile,'-d test -a','-e',enzyme)
    simRad.call
    simRad.out=as.numeric(system(simRad.call, intern=TRUE))
    percental<-simRad.out[1] #percentage of SBR in the alignment
    percentb<-simRad.out[2]
    for (k in kc) {
      system('rm -r test_all -f'))
      system('rm test_selected -f -r')
      system('mkdir test_all test_selected')
      #aaf for SBA only
      phylokmer1.call<-paste('python ~/bin/aaf_phylokmer.py -k', k, '-t 10 -d test_all -G 10')
      phylokmer1.call
      system(phylokmer1.call)
      
      distance1.call<-paste('python ~/bin/aaf_distance.py -i phylokmer.dat.gz -t 4 -f kmer_diversity.wc -o test_all')
      distance1.call
      system(distance1.call)
      
      # aaf for before seletion
      phylokmer2.call<-paste('python ~/bin/aaf_phylokmer.py -k', k, '-t 10 -d test -G 10')
      phylokmer2.call
      system(phylokmer2.call)
      
      distance2.call<-paste('python ~/bin/aaf_distance.py -i phylokmer.dat.gz -t 4 -f kmer_diversity.wc -o test')
      distance2.call
      system(distance2.call)
      
      #reads selection
      merge.call<-paste("kmer_merge -k s -c -d '0' -a A 'SP10_RAD.pkdat.gz' 'SP11_RAD.pkdat.gz' 'SP12_RAD.pkdat.gz' 'SP13_RAD.pkdat.gz' 'SP14_RAD.pkdat.gz' 'SP15_RAD.pkdat.gz' 'SP16_RAD.pkdat.gz' 'SP17_RAD.pkdat.gz' 'SP18_RAD.pkdat.gz' 'SP19_RAD.pkdat.gz' 'SP1_RAD.pkdat.gz' 'SP2_RAD.pkdat.gz' 'SP3_RAD.pkdat.gz' 'SP4_RAD.pkdat.gz' 'SP5_RAD.pkdat.gz' 'SP6_RAD.pkdat.gz' 'SP7_RAD.pkdat.gz' 'SP8_RAD.pkdat.gz' 'SP9_RAD.pkdat.gz' | cut -f 1 > test.kmer")
      merge.call
      system(merge.call)
      
      system(paste('sh /home/heather/R/RAD/readselection_test.sh'))
      
      selection.call<-paste('python ~/bin/selectedSites.py -d test_selected')
      percenta<-percenta,as.numeric(system(selection.call,intern=TRUE))
      
      #aaf on after selection
      phylokmer3.call<-paste('python ~/bin/aaf_phylokmer.py -k', k, '-t 10 -d test_selected -G 10')
      phylokmer3.call
      system(phylokmer3.call)
      
      distance3.call<-paste('python ~/bin/aaf_distance.py -i phylokmer.dat.gz -t 4 -f kmer_diversity.wc -o test_selected')
      distance3.call
      system(distance3.call)
      
      #mistakes made in each step
      sed1.call<-'sed s/_RAD//g test.tre > test.tre1'
      sed1.call
      system(sed1.call)
      
      system('rm outfile -f')
      system2("treedist", stdin = paste('treedistsetting_test')) 
      mistakeb<-as.numeric(scan('outfile', skip = 6, what=list('character','character','character','numeric'))[[4]])
      
      sed2.call<-'sed s/_RAD//g test_all.tre > test_all.tre1'
      sed2.call
      system(sed2.call)
      
      system('rm outfile -f')
      system2("treedist", stdin = paste('treedistsetting_all')) 
      mistakesba<-as.numeric(scan('outfile', skip = 6, what=list('character','character','character','numeric'))[[4]])
      
      sed3.call<-'sed s/_RAD//g test_selected.tre > test_selected.tre1'
      sed3.call
      system(sed3.call)
      
      system('rm outfile -f')
      system2("treedist", stdin = paste('treedistsetting_selected')) 
      mistakea<-as.numeric(scan('outfile', skip = 6, what=list('character','character','character','numeric'))[[4]])
      
      # Gathering information for each iteration
      row<-row + 1
      results[row,] = c(r,i,k,percentb,percenta,mistakeb,mistakea,mistakesba)
    }
  }
}

