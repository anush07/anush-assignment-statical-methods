#########################################

##Anush Khadka Final Code for R assessment


############################################
############################################
###importing the given Datasets
## Defining the my X variable as dx and Y variable as ey are provided for this Module assesment.
# Dataset X
dx = as.matrix(read.csv(file ="F:/STW7089CEM Introduction To Statistical Methods For Data Science/Assignment/X.csv", header= F))
colnames(dx) <- c("dX1","dX2","dX3","dX4")  

# Dataset Y
ey = as.matrix(read.csv(file ="F:/STW7089CEM Introduction To Statistical Methods For Data Science/Assignment/Y.csv", header = F))
colnames(ey)<-c("eY")

# Dataset Time

Time = read.csv(file ="F:/STW7089CEM Introduction To Statistical Methods For Data Science/Assignment/Time.csv", header = F, skip = 1)
Time = as.matrix(rbind(0, Time))

########################################################

## library that are used for this all Tasks.
library(ggplot2)
library(rsample)
library(matlib)

#### Task 1 - Preliminary data analysis
## Outlining the values for dx & ey (i.e. X & Y) against for time series plot
dx.ts<-ts(dx,start = c(min(Time),max(Time)),frequency =1)
ey.ts<-ts(ey,start = c(min(Time),max(Time)),frequency =1)

# ploting the Time series graph
plot(dx.ts,main = "Time series plot of X-Signal", xlab = "Time", ylab = "Input signal",col="darkred")
plot(ey.ts,main = "Time series plot of Y-Signal", xlab = "Time", ylab = "Output signal",col="blue")
################################################

### Distribution of each EEG signals
## From X-Dataset, creating a density plot of dX(i.e. X) Signal
den_dx=density(dx)
plot(den_dx,main = "Density plot of whole ax input signal",col="black")
hist(dx,freq = FALSE,main = "Histogram and Density Plot of dx-Signal",xlab="dx-Signal")
lines(den_dx,lwd=2,col="red")
rug(jitter(dx))

## From X-Dataset, creating a density plot of dX1 Signal
den_dx1=density(dx[,"dX1"])
hist(dx[,"dX1"],freq = FALSE,main = "Histogram and Density plot of dX1 signal",xlab = "dX1 Signal")
lines(den_dx1,lwd=2,col="red")
rug(jitter(dx[,"dX1"]))

## From X-Dataset, creating a density plot of dX2 signal
den_dx2=density(dx[,"dX2"])
hist(dx[,"dX2"],freq = FALSE,main = "Histogram and Density plot of dX2 signal",xlab = "dX2 Signal")
lines(den_dx2,lwd=2,col="red")
rug(jitter(dx[,"dX2"]))
## From X-Dataset, creating a density plot of dX3 signal
den_dx3=density(dx[,"dX3"])
hist(dx[,"dX3"],freq = FALSE,main = "Histogram and Density plot of dX3 signal",xlab = "dX3 Signal")
lines(den_dx3,lwd=2,col="red")
rug(jitter(dx[,"dX3"]))
## From X-Dataset, creating a density plot of dX4 signal
den_dx4=density(dx[,"dX4"])
hist(dx[,"dX4"],freq = FALSE,main = "Histogram and Density plot of dX4 signal",xlab = "dX4 Signal")
lines(den_dx4,lwd=2,col="red")
rug(jitter(dx[,"dX4"]))
## From X-Dataset, creating a density plot of ey(i.e. Y) Signal 
den_ey=density(ey)
plot(den_ey,main = "Density plot of ey",xlab = "Output Signal",col="darkblue")
hist(ey,freq = FALSE,main = "Histogram and Density plot of ey",xlab = "Output Signal")
lines(den_ey,lwd=2,col="orange")
rug(jitter(ey))

##############################################################
##Creating the correlation and Scatter Plots
#Plotting dX1 against ey
plot(dx[,"dX1"],ey,main = "Correlation betweeen dX1 and eY signal", xlab = "dX1 signal", ylab = "Output signal",col="darkorange" )
#Plotting dX2 against ey
plot(dx[,"dX2"],ey,main = "Correlation betweeen dX2 and eY signal", xlab = "dX2 signal", ylab = "Output signal",col="darkorange" )
#Plotting dX3 against ey
plot(dx[,"dX3"],ey,main = "Correlation betweeen dX3 and eY signal", xlab = "dX3 signal", ylab = "Output signal",col="darkorange" )
#Plotting dX4 against ey
plot(dx[,"dX4"],ey,main = "Correlation betweeen dX4 and eY signal", xlab = "dX4 signal", ylab = "Output signal",col="darkorange" )

#######################################################################
#######################################################################

###Task 2
##Calculating ones for binding the data
d_ones = matrix(1, length(dx)/4,1)
d_ones
### Task 2.1
## model 1
# Binding data from equation of model 1
dx_Model1<-cbind(d_ones,(dx[,"dX4"]),(dx[,"dX1"]^2),(dx[,"dX1"])^3,(dx[,"dX2"])^4,(dx[,"dX1"])^4)
dx_Model1
#calculationg thetahat of model 1
dModel1_thetahat=solve(t(dx_Model1) %*% dx_Model1) %*% t(dx_Model1) %*% ey
dModel1_thetahat

## model 2
# Binding data from equation of model 1
dx_Model2<-cbind(d_ones,(dx[,"dX4"]),(dx[,"dX1"])^3,(dx[,"dX3"])^4)
dx_Model2
#calculationg thetahat of model 2
dModel2_thetahat=solve(t(dx_Model2) %*% dx_Model2) %*% t(dx_Model2) %*% ey
dModel2_thetahat

## model 3
# Binding data from equation of model 3
dx_Model3<-cbind(d_ones,(dx[,"dX3"])^3,(dx[,"dX3"])^4)
dx_Model3
#calculationg thetahat of model 3
dModel3_thetahat=solve(t(dx_Model3) %*% dx_Model3) %*% t(dx_Model3) %*% ey
dModel3_thetahat

## model 4
# Binding data from equation of model 4
dx_Model4<-cbind(d_ones,dx[,"dX2"],dx[,"dX1"]^3,dx[,"dX3"]^4)
dx_Model4
# calculationg thetahat of model 4
dModel4_thetahat=solve(t(dx_Model4) %*% dx_Model4) %*% t(dx_Model4) %*% ey
dModel4_thetahat

## model 5
# Binding data from equation of model 5 
dx_Model5<-cbind(d_ones,dx[,"dX4"],dx[,"dX1"]^2,dx[,"dX1"]^3,dx[,"dX3"]^4)
dx_Model5
# calculationg thetahat of model 5
dModel5_thetahat=solve(t(dx_Model5) %*% dx_Model5) %*% t(dx_Model5) %*% ey
dModel5_thetahat

### Task 2.2
#Calculating Y-hat and RSS Model 1
ey_hat_m1 = dx_Model1 %*% dModel1_thetahat
ey_hat_m1
# RSS Calculating
RSS_Model_1=sum((ey-ey_hat_m1)^2)
RSS_Model_1
#Calculating Y-hat and RSS of model 2
ey_hat_m2 = dx_Model2 %*% dModel2_thetahat
ey_hat_m2
# RSS Calculating
RSS_Model_2=sum((ey-ey_hat_m2)^2)
RSS_Model_2

#Calculating Y-hat and RSS of model 3
ey_hat_m3 = dx_Model3 %*% dModel3_thetahat
ey_hat_m3
# RSS Calculating
RSS_Model_3=sum((ey-ey_hat_m3)^2)
RSS_Model_3

#Calculating Y-hat and RSS of model 4
ey_hat_m4 = dx_Model4 %*% dModel4_thetahat
ey_hat_m4
# RSS Calculating
RSS_Model_4=sum((ey-ey_hat_m4)^2)
RSS_Model_4

#Calculating Y-hat and RSS of model 5
ey_hat_m5 = dx_Model5 %*% dModel5_thetahat
ey_hat_m5
# RSS Calculating
RSS_Model_5=sum((ey-ey_hat_m5)^2)
RSS_Model_5

### Task 2.3
##Likelihood and Variance
N=length(ey)

#Calculating the Variance of Model 1
Variance_Model1=RSS_Model_1/(N-1)
Variance_Model1
#Calculating the log-likelihood of Model 1
Likelihood_Model_1=-(N/2)*(log(2*pi))-(N/2)*(log(Variance_Model1))-(1/(2*Variance_Model1))*RSS_Model_1
Likelihood_Model_1

#Calculating the Variance of Model 2
Variance_Model2=RSS_Model_2/(N-1)
Variance_Model2
#Calculating the log-likelihood of Model 2
Likelihood_Model_2=-(N/2)*(log(2*pi))-(N/2)*(log(Variance_Model2))-(1/(2*Variance_Model2))*RSS_Model_2
Likelihood_Model_2

#Calculating the Variance of Model 3
Variance_Model3=RSS_Model_3/(N-1)
Variance_Model3
#Calculating the log-likelihood of Model 3
Likelihood_Model_3= -(N/2)*(log(2*pi))-(N/2)*(log(Variance_Model3))-(1/(2*Variance_Model3))*RSS_Model_3
Likelihood_Model_3

#Calculating the Variance of Model 4
Variance_Model4=RSS_Model_4/(N-1)
Variance_Model4
#Calculating the log-likelihood of Model 4
Likelihood_Model_4= -(N/2)*(log(2*pi))-(N/2)*(log(Variance_Model4))-(1/(2*Variance_Model4))*RSS_Model_4
Likelihood_Model_4

#Calculating the Variance of Model 5
Variance_Model5=RSS_Model_5/(N-1)
Variance_Model5
#Calculating the log-likelihood of Model 5
Likelihood_Model_5= -(N/2)*(log(2*pi))-(N/2)*(log(Variance_Model5))-(1/(2*Variance_Model5))*RSS_Model_5
Likelihood_Model_5

### Task 2.4
## calculating AIC and BIC of each Model
# Model 1
de_K_Model1<-length(dModel1_thetahat)
de_K_Model1
de_AIC_Model1=2*de_K_Model1-2*Likelihood_Model_1
de_AIC_Model1
de_BIC_Model1=de_K_Model1*log(N)-2*Likelihood_Model_1
de_BIC_Model1

# MOdel 2
de_K_Model2<-length(dModel2_thetahat)
de_K_Model2
de_AIC_Model2=2*de_K_Model2-2*Likelihood_Model_2
de_AIC_Model2
de_BIC_Model2=de_K_Model2*log(N)-2*Likelihood_Model_2
de_BIC_Model2

# MOdel 3
de_K_Model3<-length(dModel3_thetahat)
de_K_Model3
de_AIC_Model3=2*de_K_Model3-2*Likelihood_Model_3
de_AIC_Model3
de_BIC_Model3=de_K_Model3*log(N)-2*Likelihood_Model_3
de_BIC_Model3

# Model 4
de_K_Model4<-length(dModel1_thetahat)
de_K_Model4
de_AIC_Model4=2*de_K_Model4-2*Likelihood_Model_4
de_AIC_Model4
de_BIC_Model4=de_K_Model4*log(N)-2*Likelihood_Model_4
de_BIC_Model4

# Model 5
de_K_Model5<-length(dModel5_thetahat)
de_K_Model5
de_AIC_Model5=2*de_K_Model5-2*Likelihood_Model_5
de_AIC_Model5
de_BIC_Model5=de_K_Model5*log(N)-2*Likelihood_Model_5
de_BIC_Model5

### Task 2.5
# Error of Model 1
de_Model1_error <- ey-ey_hat_m1
qqnorm(de_Model1_error, col = "blue",main = "QQ plot for Model-1")
qqline(de_Model1_error, col = "darkred",lwd=1)

# Error of Model 2
de_Model2_error <- ey-ey_hat_m2 
qqnorm(de_Model2_error, col = "blue",main = "QQ plot for Model-2")
qqline(de_Model2_error, col = "darkred")

# Error of Model 3
de_Model3_error <- ey-ey_hat_m3
qqnorm(de_Model3_error, col = "blue",main = "QQ plot for Model-3")
qqline(de_Model3_error, col = "darkred")

# Error of Model 4
de_Model4_error <- ey-ey_hat_m4
qqnorm(de_Model4_error, col = "blue",main = "QQ plot for Model-4")
qqline(de_Model4_error, col = "darkred")

# Error of Model 5
de_Model5_error <- ey-ey_hat_m5
qqnorm(de_Model5_error, col = "blue",main = "QQ plot for Model-5")
qqline(de_Model5_error, col = "darkred")

### Task 2.7
#install package tidymodels
install.packages("tidymodels")
library(tidymodels)

#Split-y
de_split_y<-initial_split(data = as.data.frame(ey),prop= .7)
ey_training_set<-training(de_split_y)
ey_training_set

ey_testing_set<-as.matrix(testing(de_split_y))
ey_testing_set
ey_training_data<-as.matrix(ey_training_set)
ey_training_data

#split- X
de_split_x<-initial_split(data = as.data.frame(dx),prop=.7)
dx_training_set<-training(de_split_x)
dx_training_set
dx_testing_set<-as.matrix(testing(de_split_x))
dx_testing_set
dx_testing_data<-as.matrix(dx_testing_set)
dx_testing_data
dx_training_data<-as.matrix(dx_training_set)
dx_training_data

dx_training_one=matrix(1 , length(dx_training_set$dX2),1)
dx_training_model<-cbind(dx_training_one,dx_training_set[,"dX2"],(dx_training_set[,"dX1"])^3,(dx_training_set[
  ,"dX3"])^4)
de_training_thetahat=solve(t(dx_training_model) %*% dx_training_model) %*% t(dx_training_model) %*% ey_training_data

de_training_thetahat

# MOdel Out/Prediction
ey_testing_hat = dx_testing_data %*% de_training_thetahat
ey_testing_hat
de_RSS_testing=sum((ey_testing_set-ey_testing_hat)^2)
de_RSS_testing

t.test(ey_training_data, mu=500, alternative="two.sided", conf.level=0.95)

de_C_I1=-0.2049779
de_C_I2=0.4383936
de_p2 <- plot(density(ey_training_data), col="darkorange", lwd=2, main="Distribution of Training Data")
abline(v=de_C_I1,col="darkred", lty=2)
abline(v=de_C_I2,col="red", lty=2)
abline(v=0,col="darkgreen", lty=2)
de_thetaHat_training =solve(t(dx_training_data) %*% dx_training_data) %*% t(dx_training_data) %*%
  ey_training_data
de_thetaHat_training
length(de_thetaHat_training)
de_dis_test=density(ey_training_data)
plot((de_dis_test))
plot(de_dis_test,main = "Density plot of b-Y Signal")

### Calculating Confidential interval
##(95%) Confidential interval
de_z=1.96 
de_error=((ey_testing_set-ey_testing_hat))
de_error
de_n_len=length(ey_testing_hat)
de_n_len
C_I_1= de_z*sqrt((de_error*(1-de_error))/de_n_len)
C_I_1
C_I_2= de_z*sqrt((de_error*(1+de_error))/de_n_len)
C_I_2

##
de_plot_data = data.frame(
  de_x_Val = dx_Model2,
  de_y_Val = ey,
  de_SD = sqrt(Variance_Model2)  
  )
  
#de_SD = sqrt(Variance_Model2)
#de_SD

de_plot<-ggplot(de_plot_data) +
  geom_bar( aes(x=de_x_Val.2, y=ey), stat="identity", fill="lightgreen", alpha=0.7) +
  geom_errorbar( aes(x=de_x_Val.2, ymin=ey-de_SD, ymax=ey+de_SD), width=0.2, colour="black", alpha=0.5, linewidth=1)

print (de_plot)

####Task 3
de_arr_1=0
de_arr_2=0
de_f_value=0
de_s_value=0
dModel2_thetahat
#values from thetahat
de_thetebias <- 0.483065688 
de_thetaone <- 0.143578928 
de_thetatwo <- 0.010038614 
de_thetathree <- 0.001912836
Epison <- RSS_Model_2 * 2 
num <- 100 #number of iteration
##Calculating Y-hat for performing rejection ABC
counter <- 0
for (i in 1:num) {
  range1 <- runif(1,-0.448299550,0.448299550) 
  range1
  range2 <- runif(1,-0.038109255,0.038109255)
  de_New_thetahat <- matrix(c(range1,range2,de_thetatwo,de_thetathree))
  de_New_Y_Hat <- dx_Model2 %*% de_New_thetahat
  de_new_RSS <- sum((ey - de_New_Y_Hat)^2)
  de_new_RSS
  if (de_new_RSS > Epison){
    de_arr_1[i] <- range1
    de_arr_2[i] <- range2
    counter = counter+1
    de_f_value <- matrix(de_arr_1)
    de_s_value <- matrix(de_arr_2)
  }
}

hist(de_f_value)

hist(de_s_value)

###ploting the graph
plot(de_f_value,de_s_value, col = c("red", "blue"), main = "Joint and Marginal Posterior Distribution")