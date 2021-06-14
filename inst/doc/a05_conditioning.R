## -----------------------------------------------------------------------------
library(coloc)
data(coloc_test_data)
attach(coloc_test_data) # contains D3, D4 that we will use in this vignette

## ----sens0, fig.width=8,fig.height=6------------------------------------------
library(coloc)
my.res <- coloc.abf(dataset1=D3, dataset2=D4)
class(my.res)
## print.coloc_abf
my.res
sensitivity(my.res,"H4 > 0.9")

## -----------------------------------------------------------------------------
finemap.signals(D3,method="cond")
finemap.signals(D4,method="cond")

## -----------------------------------------------------------------------------
finemap.signals(D3,method="cond",pthr=1e-20) ## too small
finemap.signals(D4,method="cond",pthr=0.1) ## too big

## -----------------------------------------------------------------------------
res <- coloc.signals(D3,D4,method="cond",p12=1e-6,pthr=1e-6)
res

## ----sens, fig.width=8,fig.height=6-------------------------------------------
sensitivity(res,"H4 > 0.9",row=1)
sensitivity(res,"H4 > 0.9",row=2)

## ---- fig.width=8,fig.height=6------------------------------------------------
finemap.signals(D3,method="mask")
finemap.signals(D4,method="mask")
resm=coloc.signals(D3,D4,method="mask",p12=1e-6,pthr=1e-6,r2thr=0.01)
resm
sensitivity(resm,"H4 > 0.9",row=1)
sensitivity(resm,"H4 > 0.9",row=2)

