## -----------------------------------------------------------------------------
library(coloc)
data(coloc_test_data)
attach(coloc_test_data)
minimum_data=D1[c("beta","varbeta","snp","position","type","sdY")]
str(minimum_data)
check_dataset(minimum_data,warn.minp=1e-10)

## -----------------------------------------------------------------------------
plot_dataset(minimum_data)

## -----------------------------------------------------------------------------
minimum_ccdata=D1[c("beta","varbeta","snp","position")]
minimum_ccdata$type="cc"
str(minimum_ccdata)
check_dataset(minimum_ccdata)

## -----------------------------------------------------------------------------
nosdY_data=D1[c("beta","varbeta","snp","position","type","N","MAF")]
str(nosdY_data)
check_dataset(nosdY_data)

## -----------------------------------------------------------------------------
nobeta_data=D1[c("MAF","snp","position","type","sdY","N")]
nobeta_data$pvalues=pnorm(-abs(D1$beta/sqrt(D1$varbeta)))*2
str(nobeta_data)
check_dataset(nobeta_data)

nobeta_ccdata=D1[c("MAF","snp","position","N")]
nobeta_ccdata$pvalues=pnorm(-abs(D1$beta/sqrt(D1$varbeta)))*2
nobeta_ccdata$type="cc"
nobeta_ccdata$s=0.5
str(nobeta_ccdata)
check_dataset(nobeta_ccdata)

## -----------------------------------------------------------------------------
str(D1$LD)

## ----echo=FALSE---------------------------------------------------------------
goodD=D1
rsign=sample(c(1,-1),length(goodD$beta),replace=TRUE)
goodD$beta=goodD$beta * rsign
goodD$LD=goodD$LD * outer(rsign,rsign,"*")
badD=goodD
rsign=sample(c(1,-1),length(goodD$beta),replace=TRUE)
badD$beta=rsign * goodD$beta

## -----------------------------------------------------------------------------
par(mfrow=c(1,2))
check_alignment(goodD)
text(40,200,"good",col="red",cex=4)
check_alignment(badD)
text(0,200,"bad",col="red",cex=4)

