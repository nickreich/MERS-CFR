## Analyze MERS data for CFR
## Nicholas Reich
## May 2014

require(RCurl)
require(coarseDataTools)
require(lme4)

## read in Rambaud's data from GitHub
url <- getURL("https://raw.githubusercontent.com/rambaut/MERS-Cases/gh-pages/data/cases.csv",ssl.verifypeer=FALSE)
dat <- read.csv(text=url, colClasses=c(onset="Date", hospitalized="Date", reported="Date", death="Date"), na.strings=c("?", "(Other)", ""))

## impute date onset from date hospitalized or date reported
time_to_hosp <- with(dat, mean(hospitalized-onset, na.rm=TRUE))

## set up data structures for easier analysis
dat$onsetWeek <- floor(as.numeric(dat$onset)/7) ## crude
dat$D <- !is.na(dat$death)
dat$N <- 1

## aggregate data by onset and gender
## columns: grp, new.times, R, D, N


## subset to only cases with onset
dat_sub <- subset(dat, !is.na(onset))

dat_for_analysis <- aggregate(dat_sub[,c("D", "N")], by=list(new.times=dat_sub$onsetWeek, grpFac=dat_sub$gender), FUN=sum)
dat_for_analysis$R <- dat_for_analysis$N - dat_for_analysis$D
dat_for_analysis$grp <- as.numeric(dat_for_analysis$grpFac)

## calculate naive relative CFR
tab <- with(dat_sub, table(gender, died))

## naive relCFR M/F
(tab[2,2]/sum(tab[2,])) / (tab[1,2]/sum(tab[1,]))

## reporting rate adjusted CFR M/F
mdl <- glm(D ~ factor(new.times) + grpFac, data=dat_for_analysis, offset=log(N), family=poisson)
rev(exp(coef(mdl)))[1]

re_mdl <- glmer(D ~ (1|new.times) + grpFac, data=dat_for_analysis, offset=log(N), family=poisson)
exp(fixef(re_mdl))

require(hglm)
mdl1 <- hglm(fixed = D ~ grpFac, random = ~1 | new.times, family = poisson(link = log), rand.family=gaussian(link=identity), data = dat_for_analysis)
exp(mdl1$fixef)

mdl2 <- hglm(fixed = D ~ grpFac, random = ~1 | new.times, family = poisson(link = log), rand.family=Gamma(link=log), data = dat_for_analysis)
exp(mdl2$fixef)

## assume a survival distribution
assumed.nu = c(0, 0.3, 0.4, 0.3) ## WRONG!

## set reporting starting values
## needs length=n-1 where n=length(unique(full.dat[,"new.times"]))
alpha.start <- rep(0, length(unique(dat_for_analysis[,"new.times"]))-1)


## estimate CFR without lag
mdl <- glm(D ~ factor(new.times) + grpFac, data=dat_for_analysis, offset=log(N), family=poisson)
exp(confint(mdl))

## estimate CFR with lag
cfr.ests <- EMforCFR(assumed.nu = assumed.nu, alpha.start.values = alpha.start,
                    full.data = dat_for_analysis, verb = FALSE, 
                    SEM.var = FALSE, ## for now
                    max.iter = 500, tol = 1e-05)
## throwing an error: maybe need to have no missing values between min and max of new.times?
