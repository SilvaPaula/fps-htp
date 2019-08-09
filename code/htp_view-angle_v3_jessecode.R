library(asreml)
library(plotrix)
library(psych)

plot = read.csv("~/Desktop/plot_18BYD2.csv", header=TRUE, check.names=FALSE)
head(plot)



##ortho = read.csv("~/Downloads/2018_Ash_BYD2_plot-level_NDVI_om.csv", header=TRUE, check.names=FALSE)
ortho =  read.csv("~/Desktop/20180423_Ash_BYD2_NDVI_om.csv", header=TRUE, check.names=FALSE)
colnames(ortho) = tolower(colnames(ortho))
head(ortho)
ortho = merge(plot, ortho, by='plot_id')
nrow(ortho)

ortho$column = as.factor(ortho$column)
ortho$range = as.factor(ortho$range)
ortho$rep = as.factor(ortho$rep)




##photo = read.csv("~/Downloads/20180516_NDVI_op.csv", header=TRUE, check.names=FALSE)
photo = read.csv("~/Desktop/20180423_Ash_BYD2_NDVI_or.csv", header=TRUE, check.names=FALSE)
colnames(photo) = tolower(colnames(photo))

nrow(photo)
photo = merge(plot, photo, by='plot_id', suffixes = c("", ".y"))
nrow(photo)
head(photo)

photo$column = as.factor(photo$column)
photo$range = as.factor(photo$range)
photo$rep = as.factor(photo$rep)

photo = cbind(photo, img=as.factor(substring(photo$original, 21,24)))
##photo$img = as.factor(photo$img)
photo = photo[order(photo$img), ]
head(photo)

photo[photo$img=='0033',]


plot(photo$mean, photo$mode)
abline(a=0, b=1, col='red')
points(photo$mean[photo$plot_name=='ARTxEVEREST_DH115'], photo$mode[photo$plot_name=='ARTxEVEREST_DH115'], lwd=3, col='red')
points(photo$mean[photo$plot_name=='ARTxEVEREST_DH116'], photo$mode[photo$plot_name=='ARTxEVEREST_DH116'], lwd=3, col='green')

# boxplot(photo$mean/photo$mode ~ photo$treatment)
# hist(photo$mode/photo$mean, breaks=50)
# plot(photo$column.x, photo$mean/photo$mode)

# plot(photo$mean, photo$mode)
# points(photo$mean[photo$treatment=='Treated'], photo$mode[photo$treatment=='Treated'], lwd=3, col='red')
# points(photo$mean[photo$treatment=='Untreated'], photo$mode[photo$treatment=='Untreated'], lwd=3, col='green')

# plot(photo$mean, photo$mode)
# points(photo$mean[photo$rep==4], photo$mode[photo$rep==4], lwd=3, col='green')
# points(photo$mean[photo$purpose=='Entry'], photo$mode[photo$purpose=='Entry'], lwd=3, col='red')

# sub = photo$mean/photo$mode<0.94 & photo$mean>0.8 & photo$mode>0.92
# points(photo$mean[sub], photo$mode[sub], lwd=3, col='green')

# head(photo[sub,], n=100)

# hist(photo$v_deg)
# hist(photo$v_deg[sub], add=TRUE, col='grey')

# boxplot(photo$v_deg ~ sub)

# par(mfrow=c(1,2))
# polar.plot(lengths=photo$h_dist, polar.pos=photo$h_deg, rp.type='s')
# polar.plot(lengths=photo$h_dist[sub], polar.pos=photo$h_deg[sub], rp.type='s')
# clr = heat.colors(100)[photo$mean*100]

# plot(photo$range.x, photo$column.x)


# boxplot(photo$)

plot(photo$h_deg, photo$mean)
##plot(photo$v_deg, photo$mean)

points(photo$h_deg[photo$plot_name=='ARTxEVEREST_DH115'], photo$mean[photo$plot_name=='ARTxEVEREST_DH115'], lwd=3, col='red')
points(photo$h_deg[photo$plot_name=='ARTxEVEREST_DH116'], photo$mean[photo$plot_name=='ARTxEVEREST_DH116'], lwd=3, col='green')

# rm(l)
# l = nls(Mode ~ mu + (abs(H_Degree)*b)^e, data=ndvi) ##, start=list(mu=0.6, H_Deg=0) ## 
# ##predict(l)
# cor(ndvi$Mode, predict(l))


# plot(photo$mean, photo$mode)

# polar.plot(lengths=photo$h_dist, polar.pos=photo$h_deg, rp.type='s')
# clr = heat.colors(100)[photo$mean*100]

# head(photo)
# photo$h_degree = abs(photo$h_degree)
# plot(photo$h_deg, photo$mean)

# hist(photo$v_deg)

##photo$column = as.factor(photo$column)
##photo$range = as.factor(photo$range)

reps=2 ## number of reps
trts=2 ## number of treatments

getVariance = function(asr, comp){
	var = summary(asr)$varcomp
	idx = which(rownames(var)==comp)
	v = var$component[idx]
	print(paste("variance component", v))
	return(v)
}



### base models using orthomosic with single value per plot

head(ortho)

asr = asreml(fixed = mean ~ 1 , random =  ~id(plot_name) + ~id(rep) + ~id(rep):(column) + ~id(rep):(range) , data=ortho, maxit=20, subset= treatment=='Untreated') 
summary(asr)$varcomp


Vg = getVariance(asr, 'plot_name')
Ve = getVariance(asr, 'units!R')

Vg/(Vg + Ve/reps)

asr.pv = predict(asr, 'plot_name')
asr.pv$avsed




# with treatment effect

asr = asreml(fixed = mean ~ 1 + treatment, random = ~id(plot_name) + ~id(treatment):(plot_name) + ~id(treatment):(rep) + ~id(treatment):(rep):(column) + ~id(treatment):(rep):(range) , data=ortho, maxit=20, start.values=TRUE) 

 
iv = asr$vparameters.table
iv$Value[2]=1
iv$Value[3]=1
iv$Constraint[2]='U'
iv$Constraint[3]='U'

asr = asreml(fixed = mean ~ 1 + treatment, random = ~id(plot_name) + ~id(treatment):(plot_name) + ~id(treatment):(rep) + ~id(treatment):(rep):(column) + ~id(treatment):(rep):(range) , data=ortho, maxit=20, G.param = iv) 

summary(asr)$varcomp

Vg = getVariance(asr, 'plot_name')
Ve = getVariance(asr, 'units!R')
Vd = getVariance(asr, 'treatment:plot_name')

Vg/(Vg + Vd/trts + Ve/(reps*trts))


v = summary(asr)$varcomp$comp
v[3]/(v[3] + v[6]/trts + v[7]/(reps*trts))

var = summary(asr)$varcomp


comp = 'plot_name'
getVariance(var, 'plot_name')


### 



### using mean value per plot 

plot.mean = aggregate(photo$mean, list(photo$plot_id), mean)
colnames(plot.mean) = c('plot_id', 'mean')
plot.mean = merge(plot, plot.mean)

plot.mean$column = as.factor(plot.mean$column)
plot.mean$range = as.factor(plot.mean$range)
plot.mean$rep = as.factor(plot.mean$rep)

head(plot.mean)

asr = asreml(fixed= mean ~ 1 , random = ~id(plot_name) + ~id(rep) + ~id(rep):(column) + ~id(rep):(range), data=plot.mean, maxit=60, subset= treatment=='Untreated')
summary(asr)$varcomp

Vg = getVariance(asr, 'plot_name')
Ve = getVariance(asr, 'units!R')

Vg/(Vg + Ve/reps)


v = summary(asr)$varcomp$comp
v[2]/(v[2]+v[3]/reps)



# averaged subsamples with treatment effect

asr = asreml(fixed = mean ~ 1 , random =  ~id(treatment) + ~id(plot_name) + ~id(treatment):(plot_name) + ~id(treatment):(rep) + ~id(treatment):(rep):(column) + ~id(treatment):(rep):(range) , data=plot.mean, maxit=20) 
summary(asr)$varcomp

v = summary(asr)$varcomp$comp
v[3]/(v[3] + v[6]/trts + v[7]/(reps*trts))

Vg = getVariance(asr, 'plot_name')
Ve = getVariance(asr, 'units!R')
Vd = getVariance(asr, 'treatment:plot_name')

Vg/(Vg + Vd/trts + Ve/(reps*trts))



v = summary(asr)$varcomp$comp
v[3]/(v[3] + v[6]/trts + v[7]/(reps*trts))




### model including within plot variance (subsampling)

# examing variance within and among plots
plot.var = aggregate(photo$mean, list(photo$plot_id), var)

hist(plot.var$x, breaks=50)
mean(plot.var$x) # within plot variance


var(plot.mean$mean) # among plot variance


### modeling within plot variance (subset of data)

head(photo)

asr = asreml(fixed = mean ~ 1, random =  ~id(plot_name) + ~id(rep) + ~id(rep):(column) + ~id(rep):(range) + ~id(rep):(range):(column):(plot_name) , data=photo, maxit=20, subset= treatment=='Untreated') 
summary(asr)$varcomp

hm = harmonic.mean(summary(photo$plot_id, maxsum=1000))
min.p = min(summary(photo$plot_id, maxsum=1000))
hist(summary(photo$plot_id, maxsum=1000))

Vg = getVariance(asr, 'plot_name')
Ve = getVariance(asr, 'units!R')
Vd = getVariance(asr, 'rep:range:column:plot_name')

Vg/(Vg + Vd/(reps) + Ve/(reps*hm))


Vg/(Vg + Vd/(reps) + Ve/(reps*min.p))


asr.pv = predict(asr, 'plot_name')
asr.pv$avsed


## ortho photo with treatment effect

asr = asreml(fixed = mean ~ 1 + treatment , random = ~id(plot_name) + ~id(treatment):(plot_name) + ~id(treatment):(rep) + ~id(treatment):(rep):(column) + ~id(treatment):(rep):(range) + ~id(treatment):(rep):(range):(column):(plot_name) , data=photo, maxit=60, workspace='2gb') 
summary(asr)$varcomp

Vg = getVariance(asr, 'plot_name')
Ve = getVariance(asr, 'units!R')
Vt = getVariance(asr, 'treatment:plot_name')
Vd = getVariance(asr, 'treatment:rep:range:column:plot_name')

Vg/(Vg + Vt/trts + Vd/(reps*trts) + Ve/(reps*trts*hm))



boxplot(photo$mean ~ photo$rep)



## including view angle


photo$h_degree = abs(photo$h_degree)
hist(photo$h_degree)
plot(photo$h_degree, photo$mean)


asr = asreml(fixed = mean ~ 1 + h_degree + treatment, random =  ~id(plot_name) + ~id(treatment):(plot_name) + ~id(rep):(treatment) + ~id(treatment):(rep):(column) + ~id(treatment):(rep):(range) + ~id(treatment):(rep):(range):(column):(plot_name) , data=photo, maxit=60, workspace='2gb') 
summary(asr)$varcomp

Vg = getVariance(asr, 'plot_name')
Ve = getVariance(asr, 'units!R')
Vt = getVariance(asr, 'treatment:plot_name')
Vd = getVariance(asr, 'treatment:rep:range:column:plot_name')

Vg/(Vg + Vt/trts + Vd/(reps*trts) + Ve/(reps*trts*hm))



## take a look at the residuals

plot(photo$time, photo$mean)


dim(photo)
length(resid(asr))

photo.res = cbind(photo, res=resid(asr))
head(photo.res)


plot(photo.res$mean, photo.res$res)
cor.test(photo.res$mean, photo.res$res)


# look at just one plot
test = '18BYD20053'
test = '18BYD20123'
sub = photo.res$plot_id==test
plot(photo.res$mean[sub], photo.res$res[sub], main=test)


# look at one entry
test = 'ARTxEVEREST_DH072'
sub = photo.res$plot_name==test
plot(photo.res$mean[sub], photo.res$res[sub], main=test)


plot(photo$time, photo$res)

## residuals compare to time stamp

plot(photo.res$time_stamp, photo.res$res)
cor.test(photo.res$time_stamp, photo.res$res)


# look at just one plot
test = '18BYD20053'
test = '18BYD20123'
sub = photo.res$plot_id==test
plot(photo.res$time[sub], photo.res$res[sub], main=test)


# look at one entry
test = 'ARTxEVEREST_DH072'
test = 'ARTxEVEREST_DH122'
sub = photo.res$plot_name==test
plot(photo.res$time[sub], photo.res$res[sub], main=test)







asr = asreml(fixed = mean ~ 1 , random =  ~idv(treatment) + ~id(plot_name) + ~id(treatment):(plot_name) + ~id(treatment):(rep) + ~id(treatment):(rep):(column) + ~id(treatment):(rep):(range) + ~id(treatment):(rep):(range):(column):(plot_name) , data=photo, maxit=20, start.values=TRUE) 
summary(asr)$varcomp

iv = asr$vparameters.table

iv$Value[1]=1

asr = asreml(fixed = mean ~ 1 , random =  ~idv(treatment) + ~id(plot_name) + ~id(treatment):(plot_name) + ~id(treatment):(rep) + ~id(treatment):(rep):(column) + ~id(treatment):(rep):(range) + ~id(treatment):(rep):(range):(column):(plot_name) , data=photo, maxit=20, G.param = iv) 



##########################################################################

# plots per image

hist(summary(photo$img, maxsum=1000), main="Plots per Image", xlab="number of plots within image")

### model with photo effect
head(photo)
photo$img

plot(photo$img, photo$mean)
plot(photo$time_stamp, photo$mean)

asr = asreml(fixed= mean ~ 1 + h_degree, random =  ~id(plot_name) + ~id(rep) + ~id(rep):(column) + ~id(rep):(range) + ~id(rep):(range):(column):(plot_name) + ~id(img), data=photo, maxit=60, subset= treatment=='Untreated') #, start.values=TRUE )
summary(asr)$varcomp

iv = asr$vparameters.table
#iv$Value[2]=1
#iv$Value[3]=1
iv$Constraint[2]='U'
iv$Constraint[6]='U'
iv

asr = asreml(fixed= mean ~ 1 + h_degree + treatment + img , random =  ~id(plot_name) + ~id(rep) + ~id(rep):(column) + ~id(rep):(range) + ~id(rep):(range):(column):(plot_name) + ~id(img), data=photo, maxit=60, subset= treatment=='Untreated', G.param=iv )
summary(asr)$varcomp

Vg = getVariance(asr, 'plot_name')
Ve = getVariance(asr, 'units!R')
Vi = getVariance(asr, 'img')
Vd = getVariance(asr, 'rep:range:column:plot_name')

Vg/(Vg + Vd/(reps) + Vi/(hm) + Ve/(reps*hm))

Vg/(Vg + Vd/(reps) + Ve/(reps*hm))


hm.img = harmonic.mean(summary(photo$img, maxsum=1000))
Vg/(Vg + Vi/(hm) + Ve/(reps*hm))


asr = asreml(fixed= mean ~ 1 + h_degree + treatment, random =  ~id(plot_name) + ~id(rep) + ~id(rep):(column) + ~id(rep):(range) + ~id(rep):(range):(column):(plot_name) , data=photo, maxit=60, subset= treatment=='Untreated')
summary(asr)$varcomp


Vg = getVariance(asr, 'plot_name')
Ve = getVariance(asr, 'units!R')
Vd = getVariance(asr, 'rep:range:column:plot_name')

Vg/(Vg + Vd/(reps) + Ve/(reps*hm))

asr.pv = predict(asr, 'plot_name')
asr.pv$avsed






asr = asreml(fixed= mean ~ 1 + h_degree, random =  ~id(plot_name) + ~id(rep) + ~id(img) + ~id(plot_name):(img), data=photo, maxit=60, subset= treatment=='Untreated', start.values=TRUE )
summary(asr)$varcomp

iv = asr$vparameters.table
iv$Constraint[]='U'
iv

asr = asreml(fixed= mean ~ 1 + h_degree + rep, random =  ~id(plot_name) + ~id(rep):(img) + ~id(rep):(img):(plot_name) + ~id(range) + ~id(column), data=photo, maxit=60, subset= treatment=='Untreated' ) #, G.param=iv 
summary(asr)$varcomp

asr = asreml(fixed= mean ~ 1 + h_degree, random =  ~id(plot_name) + ~id(rep) + ~ar1(img), data=photo, maxit=60, subset= treatment=='Untreated' ) #, G.param=iv 
summary(asr)$varcomp


Vg = getVariance(asr, 'plot_name')
Ve = getVariance(asr, 'units!R')
Vi = getVariance(asr, 'rep:img')



hi = harmonic.mean(summary(photo$img, maxsum=1000))
Vg/(Vg + Vi/reps + Ve/(reps*hi))

Vg/(Vg + Ve/(reps*hi))


