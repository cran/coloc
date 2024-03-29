library(coloc)
data(coloc_test_data)
attach(coloc_test_data)

B1=C1=D1
C1[c("type","s")]=list(type="cc",s=0.5)
B1[c("type","s")]=list(type="cc",s=2)

test_that("missing or bad required elements throws error", {
  expect_null(check_dataset(D1, req = names(D1)))
  expect_error(check_dataset(D1, req = "test"))
  ## type="cc" and s
  expect_null(check_dataset(C1))
  expect_error(check_dataset(B1))
  expect_that(check_dataset(list(), ""), throws_error())
  expect_that(check_dataset(list(beta=1,p=2,type="blah"), ""), throws_error())
})

test_that("LD matrix must have dimnames", {
  expect_null(check_ld(D3, D3$LD))
  ld_no_dimnames <- D3$LD
  attr(ld_no_dimnames, "dimnames") <- NULL
  expect_error(check_ld(D3, ld_no_dimnames))
})

test_that("issue 79", {
  d1=list(snp=letters[1:5],
          position=1:5,
          N=200000,
          MAF=runif(5)/2,
          beta=rnorm(5),
          varbeta=rep(0.01,5),
          type="cc")
  d2=list(snp=letters[1:5],
          position=1:5,
          beta=rnorm(5),
          varbeta=rep(0.01,5),
          type="quant",
          sdY=10)
  expect_error(coloc.abf(d1,d2), NA)
})
