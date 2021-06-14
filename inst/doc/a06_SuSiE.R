## -----------------------------------------------------------------------------
library(coloc)
data(coloc_test_data)
attach(coloc_test_data) ## datasets D1, D2, D3 and D4

## ----sens0, fig.width=8,fig.height=6------------------------------------------
par(mfrow=c(2,1))
plot_dataset(D3, main="Dataset D3")
plot_dataset(D4, main="Dataset D4")

## ----sens1, fig.width=8,fig.height=6------------------------------------------
my.res <- coloc.abf(dataset1=D3, dataset2=D4)
class(my.res)
## print.coloc_abf
my.res
sensitivity(my.res,"H4 > 0.9")

## -----------------------------------------------------------------------------
check_dataset(D3,req="LD")
check_dataset(D4,req="LD")

## -----------------------------------------------------------------------------
  S3=runsusie(D3)
  summary(S3)
  S4=runsusie(D4)
  summary(S4)

## -----------------------------------------------------------------------------
if(requireNamespace("susieR",quietly=TRUE)) {
  susie.res=coloc.susie(S3,S4)
  print(susie.res$summary)
}

## ----sens, fig.width=8,fig.height=6-------------------------------------------
if(requireNamespace("susieR",quietly=TRUE)) {
  sensitivity(susie.res,"H4 > 0.9",row=1,dataset1=D3,dataset2=D4)
  sensitivity(susie.res,"H4 > 0.9",row=2,dataset1=D3,dataset2=D4)
}

