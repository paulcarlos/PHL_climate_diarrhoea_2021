# open libraries
library(dlnm)
library(splines)

# load data
dat <- readRDS("sample_data_ncr.rds")

# fixed model specifications
varfun <- "ns"  #natural cubic B-spline
vardf <- 2  #degrees of freedom for variable
lagnk <- 3  #number of knots for lag
yr <- c(4,4)  #number of years 2014-2017
tlag <- 21  #maximum lag days
dfseas <- 4  #degrees of freedom for seasonality and trend control

# formula for model
fmla1 <- list(as.formula(deaths~cb+ns(day,df=dfseas*yr[1])+dow+holiday),
              as.formula(admissions~cb+ns(day,df=dfseas*yr[2])+dow+holiday))

# time lengths
ts1 <- list(seq(as.Date("2006-01-01"),as.Date("2017-12-31"),"day"),
            seq(as.Date("2014-01-01"),as.Date("2017-12-31"),"day"))

# shortened data to 2014-2017 
dat1 <- dat[dat$date %in% ts1[[2]],]

# modelling 2m temperature from ERA5-Land for mortality
argvar <- list(fun=varfun,knots=equalknots(dat1$t2m,df=vardf))
arglag <- list(fun=varfun,knots=logknots(tlag,nk=lagnk))
cb <- crossbasis(dat1$t2m,lag=tlag,argvar=argvar,arglag=arglag)
model <- glm(fmla1[[1]],dat1,family=quasipoisson)
cpred <- crosspred(cb,model=model,cen=quantile(dat1$t2m,0.5),by=0.1)

# simple plot
plot(cpred,"overall",ylim=c(0,5), xlab="2-meter temperature", ylab="relative risks", 
     ci.arg=list(density=30,col="darkgray"))
