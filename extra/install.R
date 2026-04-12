# This is an auxiliary code for install, formatting and testing the OnlineSurr package

# devtools::document()

devtools::load_all('~/Workbench/kDGLM')
# devtools::install_github('silvaneojunior/kDGLM')
library(ggplot2)
library(dplyr)
library(tidyr)
library(rlang)


devtools::load_all('.')

sd=0.01
N=100
T=20
data=data.frame(id=rep(1:N,each=T),time=rep(1:T,N),
                Y=rnorm(N*T,0,sd),S=rnorm(N*T),X=rep(rnorm(N),each=T),
                G=rep(sample(1:2,N,replace=TRUE),each=T) %>% as.factor)

out=fit.surr(Y~X,surrogate=~pol(X=~S:G,order=2),treat=G,id=id,
             data=data)

devtools::load_all('.')
summary(out)
plot(out)
plot(out,type='CPTE')
plot(out,type='delta')


usethis::use_gpl3_license()
styler::style_pkg()
devtools::document()
devtools::build_manual()
devtools::check(args = c('--as-cran'))
rmarkdown::render('vignettes/vignette.Rmd',output_file = 'vignette.html')
pkgdown::clean_site()
pkgdown::clean_cache()
pkgdown::build_site_github_pages()

pkgdown::build_home()
devtools::build_manual()
tools::Rd2pdf(".", output = "manual.pdf", clean = FALSE)

devtools::install_local('.',force=TRUE)


ggplot(dat)+
  geom_point(aes(x=time,y=y,col=trt,group=id))+
  theme_bw()


f=function(x,lag){
  sapply(1:lag,function(lag){lagged(x,lag)})
}

lm(Y~s(S,K=10),data=data) %>% summary





f <- Y ~ s(S, K = 10):G
update_s_in_formula(f,data=data)


rd
rmarkdown::render("vignettes/onlinesurr-vignette.Rmd")
