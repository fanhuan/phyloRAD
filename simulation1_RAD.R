#Set up parameters
setwd("/home/heather/R/RAD")
Sys.setenv(PATH = "/sbin:/usr/sbin:/bin:/usr/bin:/home/heather/bin")
rep<-100
infile<-'shallow2M'
enzyme <- 'EcoRI'
rate<-0.1
intree<-'shallow.tre'
system(paste('cp',intree,'intree'))
mistakesb=NULL
mistakesa=NULL
for (i in (1:rep)) {
  simRad.call<-paste('python ~/bin/simRadAlignment.py -i',infile,'-d',paste(enzyme,'_',i,sep=''),'-e',enzyme )
  system(simRad.call, ignore.stdout = T)
  sed1.call<-paste('sed s/',enzyme,'/',enzyme,'_',i,'/g EcoRI_1.sh > step1_',i,'.sh', sep='')
  system(sed1.call) #should not ignore stdout
  sed2.call<-paste('sed s/',enzyme,'/',enzyme,'_',i,'/g readselection_',enzyme,'.sh > readselection_',
                  enzyme,'_',i,'.sh', sep='')
  system(sed2.call)
  sh1.call<-paste('sh step1_',i,'.sh',sep='')
  system(sh1.call)
  sed3.call<-paste('sed s/',enzyme,'/',enzyme,'_',i,'/g EcoRI_2.sh > step2_',i,'.sh', sep='')
  system(sed3.call) #should not ignore stdout
  sh2.call<-paste('sh step2_',i,'.sh',sep='')
  system(sh2.call)
  system(paste('cp treedistsetting treedistsetting_',i,sep=''))
  echo.call <- paste('echo ',enzyme,'_',i,'.tre1 | cat >> treedistsetting_',i, sep='')
  system(echo.call)
  system('rm outfile -f')
  system2("treedist", stdin = paste('treedistsetting_',i,sep='')) 
  mistakesb[i]<-as.numeric(scan('outfile', skip = 6, what=list('character','character','character','numeric'))[[4]])
  system(paste('cp treedistsetting treedistsetting_',i,sep=''))
  echo.call <- paste('echo ',enzyme,'_',i,'_selected.tre1 | cat >> treedistsetting_',i, sep='')
  system(echo.call)
  system('rm outfile -f')
  system2("treedist", stdin = paste('treedistsetting_',i,sep='')) 
  mistakesa[i]<-as.numeric(scan('outfile', skip = 6, what=list('character','character','character','numeric'))[[4]])
}

