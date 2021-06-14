## -----------------------------------------------------------------------------
library(coloc)
data(coloc_test_data)
attach(coloc_test_data)

## -----------------------------------------------------------------------------
plot_dataset(D1)
my.res <- finemap.abf(dataset=D1)
my.res[21:30,]

## -----------------------------------------------------------------------------
my.res <- coloc.abf(dataset1=D1,
                    dataset2=D2)
print(my.res) 

## -----------------------------------------------------------------------------
subset(my.res$results,SNP.PP.H4>0.01)

## -----------------------------------------------------------------------------
o <- order(my.res$results$SNP.PP.H4,decreasing=TRUE)
cs <- cumsum(my.res$results$SNP.PP.H4[o])
w <- which(cs > 0.95)[1]
my.res$results[o,][1:w,]$snp

