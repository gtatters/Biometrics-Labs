wd<-("~/Dropbox/Biometrics Labs")
setwd(paste0(wd,"/Lab 5/SPSS files"))
l.files<-list.files(pattern=".sav")
l.files

f<-l.files[6]
f
for(f in l.files){
  library(foreign)
  d<-read.spss(f)
  # library(psych)
  # pairs.panels(d)
  # str(d)
  #d$chromosome<-gsub(" ", "", d$chromosome)
  #d$onX<-gsub(" ", "", d$onX)
  #d$frySource<-gsub(" ", "", d$frySource)
  #d$survival<-gsub(" ", "", d$survival)
  filename<-gsub(".sav", ".csv", f)
  write.csv(d, filename, row.names=F)
}

