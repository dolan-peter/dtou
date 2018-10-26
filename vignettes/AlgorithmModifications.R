## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE,fig.width=7-------------------------------------------
library(Biostrings)
library(dtou)
fn=system.file("extdata","ecoli_str_k_12_substr_mg1655.fa",package="dtou")
ecoli<-readDNAStringSet(fn)
str<-as.character(ecoli)
system.time({results<-c_vignetteExample(str,TRUE)})
n<-30
plot(0:n,results[1:(n+1)]/results[2],xlab="k-mer",ylab="Percentage of repetitiveness",type="b",main="Duplication in ecoli")

## ------------------------------------------------------------------------
system.time({
	ecoli.dtou<-dtou(str,rc=TRUE)[[1]]
	r2<-sapply(1:100,function(n){sum(ecoli.dtou>=n)})
})
all(results[-1]/2==r2)

