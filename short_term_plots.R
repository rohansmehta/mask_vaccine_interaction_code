library(ggplot2)
library(viridis)
library(deSolve)
library(reshape2)

dres <- read.csv('/Users/rohan/short_res_maxtime_500_smallgrid.csv')
dres$X <- NULL
p <- ggplot(dres[order(-dres$class),], aes(x = maxI , y = totalI)) + 
  geom_point(aes(color=factor(class))) + 
  scale_color_viridis_d(end=0.8, name='Category', labels = c('v=0, z=0', 'v=0, z>0', 'v>0')) +
  theme_bw(base_size = 16) + 
  xlab('Max Infected Frequency') + ylab('Total Infection Load') + 
  geom_hline(yintercept = 1/0.14, linetype = 'dashed') + 
  scale_y_continuous(breaks = c(0, 2, 4, 6, 1/0.14), labels = c(0, 2, 4, 6, "1/g")) + 
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8), limits = c(0, 0.8))
ggsave('short_overall.pdf', p)

ggplot(dres, aes(x = maxI , y = totalI)) + 
  geom_point(aes(color = x)) + scale_color_viridis() + facet_grid(npeaks~class)
ss <- subset(dres, class == 2 & npeaks == 2)
dres[4772,]$npeaks <- 1
dres1 <- dres
dres1$npeaks <- factor(dres1$npeaks, levels = c(0, 1, 2), labels = c('No Peaks', 'One Peak', "Two Peaks"))
dres1$class <- factor(dres1$class, levels = c(1, 2, 3), labels = c('v=0, z=0', 'v=0, z>0', "v>0"))
ggplot(dres1, aes(x = maxI , y = totalI)) + 
  geom_point(alpha=0.5) + scale_color_viridis() + 
  facet_grid(npeaks~class) + theme_bw(base_size=16) + xlab('Maximum Infected Frequency') + 
  ylab('Total Infection Load')


# OK let's take a look at each feature.
# What happens if you start by going straight down:
# 1/(g*1000) -> 0.007142857 
ss <- subset(dres, npeaks == 0)
newnpeaks <- unlist(apply(dres, 2, function(x){return(dim(findpeaks(x))[1])}))
# Note that all of these have maxI == 1/1000. The totalI is minimally 1/(0.14*1000),
# which is when the decay rate is just g. Let's see what parameters affect this,
# and also we can probably plot all of these trajectories
ggplot(ss, aes(x = m, y = totalI)) + geom_point(aes(color = x)) + scale_color_viridis() + 
  scale_y_log10() + facet_wrap(~x)
# Hmmmmmm this effect feels random. Is it the runtime?


# ires <- c()
ires1 <- c()
peaksvec <- c()
for(r1 in ss1){
  row <- ss[r1,1:7]
  parameters <- c(x=row$x, m = row$m, n = row$n, s = row$s, c = row$c, z = row$z, 
                  v = row$v)
  state <- c(1-3/1000 , 1/1000, 1/1000, 0, 1/1000, 0, 0, 0, 0, 0, 0)
  names(state) <- c('SPM', 'SAM', 'SPN', 'SAN',
                    'IPM', 'IAM', 'IPN', 'IAN',
                    'RPM', 'RAM', 'RPN')
  times <- seq(0, maxtime, by = dt)
  
  out <- ode(y = state, times = times, func = fullsystem_no_vital, parms = parameters, 
             method = 'lsoda')
  dout <- data.frame(out)
  dout$RAN <- 1-dout$SPM-dout$SAM-dout$SPN-dout$SAN-dout$IPM-dout$IAM-dout$IPN-
    dout$IAN-dout$RPM-dout$RAM-dout$RPN
  dout <- subset(dout, !is.nan(dout$SPM))
  dout$I <- dout$IPM+dout$IAM+dout$IPN+dout$IAN
  peaks <- findpeaks(dout$I)
  if(is.null(dim(peaks)[1])){npeaks <- 0} else{
    npeaks <- dim(peaks)[1]
  }
  peaksvec <- c(peaksvec, npeaks)
  ires1 <- cbind(ires1, dout$I)
}

ires1 <- cbind(ires1, times)
dires1 <- data.frame(ires1)
names(dires)[1853] <- 'time'
mdires <- melt(dires, id='time')
ggplot(mdires, aes(x = time, y = value, group = variable)) + geom_line()
# OK that's weird let's check
newnpeaks1 <- unlist(apply(dires1, 2, function(x){
  if(is.null(dim(findpeaks(x)))){
    return(0)
  } else{
    return(dim(findpeaks(x))[1])
  }
  }))
# why does findpeaks find one peak here but not in the original loop?
# let's try again?
ss1 <- which(newnpeaks == 1)


maxval <- apply(ires, 2, max)
# aha!
ss1 <- ires[,which(maxval > 0.001)]

row <- ss[r1,1:7]
parameters <- c(x=row$x, m = row$m, n = row$n, s = row$s, c = row$c, z = row$z, 
                v = row$v)
state <- c(1-3/1000, 1/1000, 1/1000, 0, 1/1000, 0, 0, 0, 0, 0, 0)
names(state) <- c('SPM', 'SAM', 'SPN', 'SAN',
                  'IPM', 'IAM', 'IPN', 'IAN',
                  'RPM', 'RAM', 'RPN')
times <- seq(0, maxtime, by = dt)

out <- ode(y = state, times = times, func = fullsystem_no_vital, parms = parameters, 
           method = 'lsoda')
dout <- data.frame(out)
dout$RAN <- 1-dout$SPM-dout$SAM-dout$SPN-dout$SAN-dout$IPM-dout$IAM-dout$IPN-
  dout$IAN-dout$RPM-dout$RAM-dout$RPN
dout <- subset(dout, !is.nan(dout$SPM))
dout$I <- dout$IPM+dout$IAM+dout$IPN+dout$IAN
md <- melt(dout, id = 'time')
ggplot(subset(md, variable == 'I'), aes(x = time, y = value, group = variable, color = variable)) + 
  geom_line() + scale_color_viridis_d() 


peaks <- findpeaks(dout$I)



# Short-term dynamics
xvec <- seq(-1, 1, 0.2)
zvec <- seq(0, 1, 0.2)
vvec <- seq(0, 1, 0.2)
svec <- seq(0, 1, 0.2)
mvec <- svec


params <- expand.grid(xvec, zvec, vvec, svec, mvec)
names(params) <- c('x', 'z', 'v', 's', 'm')
params$c <- 1-params$s
params$n <- 1-params$m
g <- 0.14
R0 <- 15.28
r <- g*R0
b <- 0
fullsystem_no_vital <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    PM <- SPM+IPM+RPM
    PN <- SPN+IPN+RPN
    P <- PM+PN
    A <- 1 - P
    AM <- SAM+IAM+RAM
    M <- PM + AM
    N <- 1 - M
    AN <- N - PN
    I <- IPM + IPN + IAM + IAN
    RAN <- AN-SAN-IAN
    dspm <- (1+x)*(s*P*SAM+m*M*SPN)-SPM*(r*z*I+(1-x)*(c*A+n*N)+v)
    dsam <- (1-x)*(c*A*SPM+m*M*SAN)-SAM*(r*z*I+(1+x)*(s*P+n*N))
    dspn <- (1-x)*(s*P*SAN+n*N*SPM)-SPN*(r*I+(1+x)*(c*A+m*M)+v)
    dsan <- (1+x)*(c*A*SPN+n*N*SAM)-SAN*(r*I+(1-x)*(s*P+m*M))
    dipm <- r*z*I*SPM+(1+x)*(s*P*IAM+m*M*IPN)-IPM*(g+(1-x)*(c*A+n*N))
    diam <- r*z*I*SAM+(1-x)*(c*A*IPM+m*M*IAN)-IAM*(g+(1+x)*(s*P+n*N))
    dipn <- r*I*SPN+(1-x)*(s*P*IAN+n*N*IPM)-IPN*(g+(1+x)*(c*A+m*M))
    dian <- r*I*SAN+(1+x)*(c*A*IPN+n*N*IAM)-IAN*(g+(1-x)*(s*P+m*M))
    drpm <- v*SPM+g*IPM+(1+x)*(s*P*RAM+m*M*RPN)-RPM*((1-x)*(c*A+n*N))
    dram <- g*IAM+(1-x)*(c*A*RPM+m*M*RAN)-RAM*((1+x)*(s*P+n*N))
    drpn <- v*SPN+g*IPN+(1-x)*(s*P*RAN+n*N*RPM)-RPN*((1+x)*(c*A+m*M))
    list(c(dspm, dsam, dspn, dsan, dipm, diam, dipn, dian, drpm, dram, drpn))
  })
}
maxtime <- 500
dt <- 0.1

