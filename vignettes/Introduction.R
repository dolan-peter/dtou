## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(dtou)
library(kableExtra)

## ------------------------------------------------------------------------
str<-"AAAAACCCGACTGGGCTCA"
(results<-dtou(str,rc=TRUE)[[1]])

## ----echo=FALSE----------------------------------------------------------
tmp=strsplit(str,'')[[1]]
names(results)=tmp
kable(t(results))

## ------------------------------------------------------------------------
(results<-dtou(str,rc=FALSE)[[1]])

## ----echo=FALSE----------------------------------------------------------
tmp=strsplit(str,'')[[1]]
names(results)=tmp
kable(t(results))

## ----echo=FALSE,output=FALSE,message=FALSE-------------------------------
#Here is where we actually read the data into the file
library(Biostrings)
fn1=system.file("extdata","ecoli_str_k_12_substr_mg1655.fa",package="dtou")
fn2=system.file("extdata","ecoli_O157_H7_str_Sakai.fa",package="dtou")
ecoli=readDNAStringSet(fn1)
ecoli2=readDNAStringSet(fn2)
#ecoli=readDNAStringSet("../_data/ecoli_str_k_12_substr_mg1655.fa")
#ecoli2=readDNAStringSet("../data/ecoli_O157_H7_str_Sakai.fa")

## ----eval=FALSE----------------------------------------------------------
#  library(Biostrings)
#  ecoli=readDNAStringSet("data/ecoli_str_k_12_substr_mg1655.fa")
#  ecoli2=readDNAStringSet("data/ecoli_O157_H7_str_Sakai.fa")

## ------------------------------------------------------------------------
str<-as.character(ecoli)
str2<-as.character(ecoli2)
system.time({ecoli.dtou<-dtou(str)[[1]]})
system.time({ecoli2.dtou<-dtou(str)[[1]]})

## ----fig.width=7---------------------------------------------------------
plot(ecoli.dtou,
		 pch=".",
		 main="Distance to Uniqueness for base-pairs in Ecoli sequence",
		 xlab="index",ylab="dtou"
)

## ----fig.width=7---------------------------------------------------------
window.length=7000
i=which.max(ecoli.dtou)
lower=i-window.length
upper=i+window.length
zoom=lower:upper
head(zoom)
plot(zoom,
		 ecoli.dtou[zoom],
		 pch=".",
		 main="Distance to Uniqueness for base-pairs in Ecoli sequence (zoomed)",
		 xlab="index",ylab="dtou"
)

## ------------------------------------------------------------------------
quantile(ecoli.dtou,1:9*0.1)

## ------------------------------------------------------------------------
sum(ecoli.dtou>20)
sum(ecoli.dtou>20)/length(ecoli.dtou)

## ----fig.width=7---------------------------------------------------------
hist(log10(ecoli.dtou[ecoli.dtou>20]),main="E.coli Upper Tail Distance to Uniqueness",xlab="log10 dtou",sub="dtou > 20")

## ------------------------------------------------------------------------
system.time({ecoli_both.dtou<-dtou(c(str,str2),optimizeForSpeed = TRUE)})

## ----fig.width=7---------------------------------------------------------
difference=ecoli_both.dtou[[1]]-ecoli.dtou
plot(log10(difference+1),pch=".", ylab="log10 dtou",xlab="sequence index",main="E. coli K12 MG1655 vs E. coli O157:H7 Sakai")

## ------------------------------------------------------------------------
sum(difference!=0)/length(difference)

## ----fig.width=7---------------------------------------------------------
hist(log10(difference),main="Comparison of E. coli K12 MG1655 to E. coli O157:H7 Sakai",xlab="log10 of diffrence in dtou")

## ----eval=FALSE----------------------------------------------------------
#  chr1<-readDNAStringSet("data/hs_ref_GRCh38.p12_chr1.fa",format="fasta")

## ----eval=FALSE----------------------------------------------------------
#  chr1<-chr1[!grepl("alternate",names(chr1))]
#  chr1=as.character(chr1)

## ----eval=FALSE----------------------------------------------------------
#  substr(chr1[[12]],25143579,25143579)<-"X"
#  substr(chr1[[12]],25146187,25146187)<-"X"

## ----eval=FALSE----------------------------------------------------------
#  system.time({results<-dtou(chr1,depth=10000,optimizeForSpeed=TRUE)})
#  system.time({results<-dtou(chr1,depth=10000,optimizeForSpeed=FALSE)})

## ----eval=FALSE----------------------------------------------------------
#  dtos<-shrinkToMinimum(results)

## ----echo=FALSE----------------------------------------------------------
#Load the data previously calculated and saved
#load("../notpackage/chr1.rdata")

## ----eval=FALSE----------------------------------------------------------
#  q=sapply(dtos,function(d){quantile(d,1:9*0.1)})

## ----echo=FALSE----------------------------------------------------------
#Load the data previously calculated and saved... the raw data is too large for an R package
fn=system.file("extdata","quantile.rdata",package="dtou")
load(fn) #Loads q

## ----results='asis'------------------------------------------------------
kable(t(q))

## ----eval=FALSE----------------------------------------------------------
#  df=readDNAStringSet("data/Psme_TSA.fasta")
#  ref=as.character(df)
#  df.dtou=dtou(ref,optimizeForSpeed = TRUE)

## ----eval=FALSE----------------------------------------------------------
#  l=sapply(df.dtou,length)
#  max=sapply(df.dtou,max)

## ----echo=FALSE----------------------------------------------------------
fn=system.file("extdata","psme.rdata",package="dtou")
load(fn) #Loads q

## ----fig.width=7---------------------------------------------------------
plot(max/l,pch=".",main="Percentage of transcript covered by longest repeat",ylab="Percentage")

## ----fig.width=7---------------------------------------------------------
plot(sort(max/l,decreasing=TRUE),pch=".",main="Percentage of transcript covered by longest repeat",ylab="Percentage")

## ------------------------------------------------------------------------
sum(max/l==1)

