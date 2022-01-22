#Power simulation
library(lme4)
library(lmerTest)

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Research/PHD_Bio_CSU/CH1_YEWA_telomeres/Analysis")

YEWA_TL_GV_csv=read.csv("YEWA_TL_GV.csv")

#Y: telomere length
#X1: sex
#X2: vulnerabiltiy

#how many sites?
nsite <- 10
#how many obs per site?
site_reps <- 50
#site numbers
site_no <- rep(1:nsite, each = site_reps)
#total n
n <- nsite*site_reps
#simulate sex
X1 <- sample(c(0, 1), n, replace = T)
#simulate vulnerability
X2 <- runif(n, 0.002, 0.018)

#site-to-site variation
sigma_site <- 0.2

#noise variation
sigma <- 0.005

#slope for Y-X2 relationship
beta1 <- -0.5
#intercept for sex = 0
beta0 <- 1.5

#gender effect (diff in intercept between sex's)
#create a sequence of effect sizes
ge <- seq(0, 1, length = 0.01)

#save results
power <- rep(NA, length(ge))

#how many simulations for each effect size?
nsim <- 100


#loop over effect sizes (could replace this with observations per site)
for(i in 1:length(ge)){
	ge_val <- ge[i]
	#store p-values
	pvals <- rep(NA, nsim)
	#simulate random site effect
	site_effs <- rnorm(nsite, mean = 0, sigma_site)
	site_effs <- rep(site_effs, each = site_reps)
	
	for(j in 1:nsim){
		#simulate data at each site
		#specify mean at each site, for each gender
		mu <- site_effs + beta0 + ge_val*X1 + beta1*X2
		#simulate observations
		Y <- rnorm(n, mean = mu, sd = sigma)
		#create a data set
		simdata <- data.frame("site" = as.factor(site_no), Y, X1, X2 )
		#estimate the mixed effects model
		lmerOut <- lmer(Y ~ (1|site) + X1 + X2, data = simdata)
		#get the p-value for sex, save it
		pvals[j] <- summary(lmerOut)$coef[2, 3]
	}
	#compute power
	power[i] <- sum(pvals < 0.05)/nsim
}

#plot power vs effect size

plot(ge, power, type = "b", xlab = "effect size", ylab = "power")


#plot the data from the last simulation
plot(simdata$X2, simdata$Y, col = simdata$X1+1)
