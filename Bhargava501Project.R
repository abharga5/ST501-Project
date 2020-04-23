#Author: Ashwin Bhargava
#Due Date: 7/25/2019
#Purpose: Use R to simulate random data, approximate quantities,
#and create graphs. Helps understand the convergence concepts
#of Chapter 5.

#Problem 1)

#Parts A, B, and C
#See PDF

#Parts D and E)

#Generate N=50 Values of L and K for every value of n from 1 to
#250. So consider n=1. We want to generate N=50 datasets with n=1.
#For each data set we want to find L and K. Now, for n=2 we want
#to generate 50 datasets, finding L and K for each one. All the way
#until we get to n=250.

n<-seq(1:250)
N<-50

mu<-0
b<-5

#Create Inverse Laplace CDF Function
InvLaplace<-function(mu,b,u){
  
  InvEq<-mu - b*sign(u-.5)*log(1-2*abs(u-.5))
  return(InvEq)
}

#L --> 50

L_matrix<-matrix(data=0, nrow=N, ncol=250)

#function to calculate L

L_val=function(ui){
  
  n=length(ui)
  L=(sum((ui^2))/n)
  return(L)
}

#K -->~ 7.07

K_matrix<-matrix(data=0, nrow=N, ncol=250)

#function to calculate K (from L)

K_val=function(ui){
  
  K=sqrt(L_val(ui))
  return(K)
}

for (j in 1:250){
  # Loop through data sets
  for (i in 1:N){
    u<-runif(n[j],min=0,max=1) #generates random values
    
    L_matrix[i, j] = L_val(InvLaplace(mu=0, b=5, u))
    K_matrix[i, j] = K_val(InvLaplace(mu=0, b=5, u))
  }
}

#plot L estimates

plot(n, L_matrix[1,], main="L Value Estimates",
     ylab="Estimate of L",ylim=c(0,110))

#add lines

abline(a=50,b=0,col="red")
abline(a=50-20,b=0,lty=2,col="blue")
abline(a=50+20,b=0,lty=2,col="blue")

#add legend

legend("bottomright",legend=c("Theoretical Convergence",
      "Epsilon Bounds"),col=c("red","blue"),lty=1.2,cex=0.8)

#plot K estimates

plot(n, K_matrix[1,], main="K Value Estimates",
     ylab="Estimate of K",ylim=c(0,12))

#add lines

abline(a=sqrt(50),b=0,col="red")
abline(a=sqrt(50)-3,b=0,lty=2,col="blue")
abline(a=sqrt(50)+3,b=0,lty=2,col="blue")

#add legend

legend("bottomright",legend=c("Theoretical Convergence",
      "Epsilon Bounds"),col=c("red","blue"),lty=1.2,cex=0.8)

# Q 1E) Part v.

#Graphically, we are seeing the theoretical result of convergence
#in probability. We can observe this by looking at the spread
#around n = 0 versus the spread around n = 250. By LLN,
#as n increases, the random variables converge to the expected
#mean, denoted with the red line. We would expect the spread to
#continue shrinking if n were increased even more.


#Problem 2)

#Take the N values of L you have for each of n=10,50,100,250 and
#create a histogram to inspect the shape of the sampling
#distribution of L for those values of n.

#histograms of chosen n values (columns used from L_matrix)

hist(L_matrix[,10],main="Distribution of the Sample Means for L (n=10)",
    xlab="Sample Means (n=10)", ylab="Frequency",ylim=c(0,25))
hist(L_matrix[,50],main="Distribution of the Sample Means for L (n=50)",
     xlab="Sample Means (n=50)", ylab="Frequency",ylim=c(0,25))
hist(L_matrix[,100],main="Distribution of the Sample Means for L (n=100)",
     xlab="Sample Means (n=100)", ylab="Frequency",ylim=c(0,25))
hist(L_matrix[,250],main="Distribution of the Sample Means for L (n=250)",
xlab="Sample Means (n=250)", ylab="Frequency",ylim=c(0,25))

#Based on the histograms, there is a weak to no trend towards convergence
#to a normal distribution for L (or K). n=250 samples makes for a poor
#fit in a normal distribution. We will need many more samples in order to 
#see an actual "bell-shape" begin to form.

#Problem 3)
#See PDF

#Problem 4)

#Standardize datasets with result from 3)


hist((L_matrix[,10]-2*b^2)/((2*b^2*sqrt(5))/sqrt(10)),
     main="Standardized Distribution of L (n=10)",
     xlab="Sample Means (n=10)", ylab="Frequency",ylim=c(0,25))
hist((L_matrix[,50]-2*b^2)/((2*b^2*sqrt(5))/sqrt(50)),
     main="Standardized Distribution of L (n=50)",
     xlab="Sample Means (n=50)", ylab="Frequency",ylim=c(0,25))
hist((L_matrix[,100]-2*b^2)/((2*b^2*sqrt(5))/sqrt(100)),
     main="Standardized Distribution of L (n=100)",
     xlab="Sample Means (n=100)", ylab="Frequency",ylim=c(0,25))
hist((L_matrix[,250]-2*b^2)/((2*b^2*sqrt(5))/sqrt(250)),
     main="Standardized Distribution of L (n=250)",
     xlab="Sample Means (n=250)", ylab="Frequency",ylim=c(0,25))

#Problem 5)

#n=1000, 10000 with N = 10,000 data sets

n<-c(1000, 10000)
N<-10000

mu<-0
b<-5

#For lack of a better matrix name..

L_big_matrix<-matrix(data=0, nrow=N, ncol=2)

for (j in 1:length(n)){
  # Loop through data sets
  for (i in 1:N){
    u<-runif(n[j],min=0,max=1)
    L_big_matrix[i, j] = L_val(InvLaplace(mu=0, b=5, u))
  }
}

#Problem 6) 

#Two histograms are Sampling Distribution of Sample Mean

hist(L_big_matrix[,1],main="Distribution of the Sample Means for L (n=1000)",
     xlab="Sample Means (n=1000)", ylab="Frequency",ylim=c(0,2500))
hist(L_big_matrix[,2],main="Distribution of the Sample Means for L (n=10000)",
     xlab="Sample Means (n=10000)", ylab="Frequency",ylim=c(0,2500))

#Two histograms are standardized distributions

hist((L_big_matrix[,1]-2*b^2)/((2*b^2*sqrt(5))/sqrt(1000)),
     main="Standardized Distribution of L (n=1000)",
     xlab="Sample Means (n=1000)", ylab="Frequency",ylim=c(0,2500))
      
hist((L_big_matrix[,2]-2*b^2)/((2*b^2*sqrt(5))/sqrt(10000)),
     main="Standardized Distribution of L (n=10000)",
     xlab="Sample Means (n=10000)", ylab="Frequency",ylim=c(0,2500))

#The larger n values tend to give a fairly 'normal'
#histogram compared to the previously small n values. We
#can safely say that n = 30 does not work for this
#distribution. Laplace distribution in particular has
#larger tails (higher kurtosis) which leads to more
#outliers than a normal distribution. So in order to apply
#the CLT to a Laplace distribution, we need many more
#samples than 30.
