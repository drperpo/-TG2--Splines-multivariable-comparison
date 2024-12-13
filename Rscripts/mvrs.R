library(splines)
# step 1 fit linear model

fit_linear <- glm(Outcome ~ Pregnancies + Glucose + BloodPressure + SkinThickness +
                Insulin + BMI+ 
                DiabetesPedigreeFunction + Age ,
              family = "binomial", data=dat1)

# rank variables in terms of p-values
res_table <- summary(fit_linear)$coefficients[-1,]
rank_vars <- rank(abs(res_table[,3]))
vars_order <- rev(row.names(res_table)[order(rank_vars)])



it <- 5

outs <-  data.frame(Var=character(),
                 df=character(), 
                 stringsAsFactors=FALSE) 

for(it in 1:length(vars_order)){  

d_f <- 5
vars_test <- vars_order
vars_test[it] <- paste("ns(",vars_order[it],",df=d_f)",sep="")
spline_frm <- as.formula(paste("Outcome",paste(vars_test, collapse=" + "), sep="~"))
fit_spline <- glm(spline_frm,  family = "binomial", data=dat1)





### linear model
frm <- as.formula(paste("Outcome", paste(vars_order, collapse=" + "), sep="~"))
fit_linear <- glm(frm,  family = "binomial", data=dat1)

### null  model 
vars_null <- vars_order[-it]
null_frm <- as.formula(paste("Outcome",paste(vars_null, collapse=" + "), sep="~"))
fit_NULL <- glm(null_frm,family = "binomial", data=dat1)



cri <- fit_NULL$deviance-fit_spline$deviance> qchisq(.8, df=d_f+1)
if(!cri){out <- data.frame( Var=vars_order[it], df=0)}
if(cri){
# for smaller p-value fit spline

cri <- TRUE


while(cri & d_f>=1){
  ### splines model
  vars_test <- vars_order
  vars_test[it] <- paste("ns(",vars_order[it],",df=d_f)",sep="")
  spline_frm <- as.formula(paste("Outcome",paste(vars_test, collapse=" + "), sep="~"))
  fit_spline <- glm(spline_frm,  family = "binomial", data=dat1)


# Spline vs NULL
cri <- fit_linear$deviance-fit_spline$deviance< qchisq(.8, df=d_f+1) 
cat(d_f,"   ", cri, "\n")
if(cri){d_f <- d_f-1}
 if (d_f==1){
   break}
}



out <-  data.frame( Var=vars_order[it], df=d_f)

}


outs <- rbind(outs, out)
}
outs

fit_mvss <- glm(Outcome ~ Pregnancies + Glucose  +
                     ns(BMI, df=5) + 
                    ns(DiabetesPedigreeFunction, df=2)+
                ns(Age,df=5),
                  family = "binomial", data=dat1)


###### plots
plotdata <- plot(fit_ps_mgcv, pages = 1)
tplot <- termplot(fit_mvss, pages=1, plot=FALSE)


cols <- c(
"#e41a1c",
"#377eb8",
"#4daf4a",
"#984ea3",
"#ff7f00",
"#ffff33",
"#a65628",
"#f781bf")

###  Glucose
glucose <- plotdata[[2]]$x
plot(glucose, 3.53287*(glucose/100)-3.53287*(50/100),'l',
     ylab="Partial Predictor", xlab="Glucose")
lines(glucose, 0.03377   *glucose-0.03377   *50, col="#e41a1c")
lines(glucose, plotdata[[2]]$fit+2.466689, col="#377eb8")
lines(glucose, plotdata[[2]]$fit+2.466689+plotdata[[2]]$se,col="#377eb8", lty=2)
lines(glucose, plotdata[[2]]$fit+2.466689-plotdata[[2]]$se, col="#377eb8",lty=2)



###  BMI

BMI <-  tplot$BMI[,1]
plot(BMI, -0.16807*(BMI/100)^(-2)-0.16807*(min(BMI)/100)^(-2),'l',ylim=c(-12,0), ylab="Partial Predictor", xlab="BMI")
lines(BMI, tplot$BMI[,2]-6.638834, col=2)
lines(plotdata[[6]]$x, plotdata[[6]]$fit-6.960142,col=3)
lines(plotdata[[6]]$x, plotdata[[6]]$se+plotdata[[6]]$fit-6.960142,col=3,lty=2)
lines(plotdata[[6]]$x,- plotdata[[6]]$se+ plotdata[[6]]$fit-6.960142,col=3,lty=2)



###  Diabetes
diabetes <- tplot$DiabetesPedigreeFunction[,1]
plot(diabetes, 0.82126* diabetes-0.82126* 0.078, 'l',ylim=c(-2,2.5), ylab="Partial Predictor",
     xlab="Diabetes")
lines(diabetes, tplot$DiabetesPedigreeFunction[,2]-tplot$DiabetesPedigreeFunction[1,2], col=2)
lines(plotdata[[7]]$x, plotdata[[7]]$fit-plotdata[[7]]$fit[1],col=3)
lines(plotdata[[7]]$x, plotdata[[7]]$se+plotdata[[7]]$fit-plotdata[[7]]$fit[1],col=3,lty=2)
lines(plotdata[[7]]$x, -plotdata[[7]]$se+plotdata[[7]]$fit-plotdata[[7]]$fit[1],col=3,lty=2)



### Age
Age <- tplot$Age[,1]
plot(Age, -0.06563*I((Age/100)^-2)+0.06563*I((21/100)^-2),'l',ylim=c(-2,2.2),ylab="Partial Predictor",
     xlab="Age")
lines(Age, tplot$Age[,2]-tplot$Age[1,2], col=2)

lines(plotdata[[8]]$x, plotdata[[8]]$fit-plotdata[[8]]$fit[1], col=3)
lines(plotdata[[8]]$x, plotdata[[8]]$se+plotdata[[8]]$fit-plotdata[[8]]$fit[1], col=3,lty=2)
lines(plotdata[[8]]$x, -plotdata[[8]]$se+plotdata[[8]]$fit-plotdata[[8]]$fit[1], col=3,lty=2)








#vdist <- hdist <- 0.2
ylim <- c(-4,4)
layout(matrix(1:4, 2, 2, byrow=TRUE),widths=c(10,10),heights=c(10,10))
par(mar= c(4, 4, 3, 0.2))

###  Glugose
glucose <- plotdata[[2]]$x
plot(glucose, 3.53287*(glucose/100)-3.53287*(50/100),'l',
     ylab="Partial Predictor", xlab="Glucose", axes=FALSE)
lines(glucose, 0.03377   *glucose-0.03377   *50, col="#e41a1c")
lines(glucose, plotdata[[2]]$fit+2.466689, col="#377eb8")
lines(glucose, plotdata[[2]]$fit+2.466689+plotdata[[2]]$se,col="#377eb8", lty=2)
lines(glucose, plotdata[[2]]$fit+2.466689-plotdata[[2]]$se, col="#377eb8",lty=2)
box(); axis(1);axis(2)


###  BMI


par(mar= c(4, .5, 3, 4))

BMI <-  tplot$BMI[,1]
plot(BMI, -0.16807*(BMI/100)^(-2)-0.16807*(min(BMI)/100)^(-2),'l',
    axes=FALSE, ylim=c(-12,0), ylab="Partial Predictor", xlab="BMI")
lines(BMI, tplot$BMI[,2]-6.638834, col="#e41a1c")
lines(plotdata[[6]]$x, plotdata[[6]]$fit-6.960142,col=3)
lines(plotdata[[6]]$x, plotdata[[6]]$se+plotdata[[6]]$fit-6.960142, col="#377eb8",lty=2)
lines(plotdata[[6]]$x,- plotdata[[6]]$se+ plotdata[[6]]$fit-6.960142, col="#377eb8",lty=2)
box();axis(1)



par(mar= c(4, 4, 1.5, 0.2))

###  Diabetes
diabetes <- tplot$DiabetesPedigreeFunction[,1]
plot(diabetes, 0.82126* diabetes-0.82126* 0.078, 'l',ylim=c(-2,2.5), ylab="Partial Predictor",
     xlab="Diabetes", axes=FALSE)
lines(diabetes, tplot$DiabetesPedigreeFunction[,2]-
        tplot$DiabetesPedigreeFunction[1,2], col="#e41a1c")
lines(plotdata[[7]]$x, plotdata[[7]]$fit-plotdata[[7]]$fit[1],col="#377eb8")
lines(plotdata[[7]]$x, plotdata[[7]]$se+plotdata[[7]]$fit-plotdata[[7]]$fit[1],
      col="#377eb8",lty=2)
lines(plotdata[[7]]$x, -plotdata[[7]]$se+plotdata[[7]]$fit-plotdata[[7]]$fit[1],
      col="#377eb8",lty=2)
box(); axis(1);axis(2)



### Age
par(mar= c(4, .5, 1.5, 4))

Age <- tplot$Age[,1]
plot(Age, -0.06563*I((Age/100)^-2)+0.06563*I((21/100)^-2),'l',ylim=c(-2,2.2),ylab="Partial Predictor",
     xlab="Age",axes=FALSE)
lines(Age, tplot$Age[,2]-tplot$Age[1,2], col="#e41a1c")
lines(plotdata[[8]]$x, plotdata[[8]]$fit-plotdata[[8]]$fit[1], col="#377eb8")
lines(plotdata[[8]]$x, plotdata[[8]]$se+plotdata[[8]]$fit-plotdata[[8]]$fit[1],
      col="#377eb8",lty=2)
lines(plotdata[[8]]$x, -plotdata[[8]]$se+plotdata[[8]]$fit-plotdata[[8]]$fit[1], 
      col="#377eb8",lty=2)
box(); axis(1)


####### compare mgcv

mg1 <- plot(fit_tp_mgcv_noS, pages=1)
mg2 <- plot(fit_ts_mgcv, pages=1)
mg3 <- plot(fit_tp_mgcv, pages=1)
mg4 <- plot(fit_ps_mgcv, pages=1)
mg5 <- plot(fit_cs_mgcv, pages=1)



# Glucose
glucose <- mg1[[2]]$x
plot(glucose, mg1[[2]]$fit,'l', ylim=c(-4,3),
     ylab="Partial Predictor", xlab="Glucose")
lines(glucose, mg2[[2]]$fit, col=cols[1])
lines(glucose, mg3[[2]]$fit, col=cols[2])
lines(glucose, mg4[[2]]$fit, col=cols[3])
lines(glucose, mg5[[2]]$fit, col=cols[4])



# BMI
 mg1[[6]]$xlab
plot( mg1[[6]]$x, mg1[[6]]$fit,'l', 
     ylab="Partial Predictor", xlab= mg1[[6]]$xlab)
lines(mg1[[6]]$x, mg2[[6]]$fit, col=cols[1])
lines(mg1[[6]]$x, mg3[[6]]$fit, col=cols[2])
lines(mg1[[6]]$x, mg4[[6]]$fit, col=cols[3])
lines(mg1[[6]]$x, mg5[[6]]$fit, col=cols[4])




# Diabetes

mg1[[7]]$xlab
plot( mg1[[7]]$x, mg1[[7]]$fit,'l', 
      ylab="Partial Predictor", xlab= mg1[[7]]$xlab)
lines(mg2[[7]]$x, mg2[[7]]$fit, col=cols[1])
lines(mg1[[7]]$x, mg3[[7]]$fit, col=cols[2])
lines(mg1[[7]]$x, mg4[[7]]$fit, col=cols[3])
lines(mg1[[7]]$x, mg5[[7]]$fit, col=cols[4])

# Age

mg1[[8]]$xlab
plot( mg1[[8]]$x, mg1[[8]]$fit,'l', 
      ylab="Partial Predictor", xlab= mg1[[8]]$xlab)
lines(mg2[[8]]$x, mg2[[8]]$fit, col=cols[1])
lines(mg1[[8]]$x, mg3[[8]]$fit, col=cols[2])
lines(mg1[[8]]$x, mg4[[8]]$fit, col=cols[3])
lines(mg1[[8]]$x, mg5[[8]]$fit, col=cols[4])
