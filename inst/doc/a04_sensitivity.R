## ----prep---------------------------------------------------------------------
library(coloc)
data(coloc_test_data)
attach(coloc_test_data)
my.res <- coloc.abf(dataset1=D1,
                    dataset2=D2,
                    p12=1e-6)
my.res

## ----sens, fig.width=8,fig.height=6-------------------------------------------
sensitivity(my.res,rule="H4 > 0.5") 

## ----echo=FALSE,results="hide"------------------------------------------------
D1a=D1;
D1a$varbeta=D1$varbeta * sqrt(4)
D1a$N=D1a$N/4
D2a=D2;
D2a$varbeta=D2$varbeta * sqrt(2)
D2a$N=D2a$N/2

## ----sens2, fig.width=8,fig.height=6------------------------------------------
my.res <- coloc.abf(dataset1=D1a,
                    dataset2=D2a,
                    p12=1e-6)
my.res
sensitivity(my.res,rule="H4 > 0.5") 

## ----sens3, fig.width=8,fig.height=6------------------------------------------
sensitivity(my.res,rule="H4 > 3*H3 & H0 < 0.1") 

