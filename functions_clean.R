### One Sample Functions ####
##ss is a vector of sample sizes
##xs is a vector of sample variances
##k is the total number of studies (length of ss or xs)
##bhat is a value for the MMLE of ahat


alpha_hat <- Vectorize(function(ss,xs,k,bhat){
  dvec <- 2*bhat+xs*(ss-1)
  sum((ss-1)/dvec)/((k/bhat)-2*sum(1/dvec))
},vectorize.args="bhat")

##root_eq is the equation that we need to find the root of to get 
##the MLE of alpha and beta

root_eq <-Vectorize(function(bhat,ss,xs,k){
  alpha<-alpha_hat(ss=ss,xs=xs,k=k,bhat=bhat)
  k*(log(bhat)-digamma(alpha))+
    sum(digamma(alpha+(ss-1)/2)-log(bhat+(xs*(ss-1)/2)))
},vectorize.args="bhat")

##Check parameters ####

Find_alphabeta<-function(ss,xs,k){
  bhat<-tryCatch(uniroot(root_eq,interval=c(0.01,100000),ss=ss,xs=xs,k=k)$root,error=function(e) 0)
  ahat<-alpha_hat(ss,xs,k,bhat)
  c(ahat,bhat)
}
## Note that this function outputs (0,0) when the boundary case occurs. We 
## should then accept alpha/beta where beta is very large as the only value
## for the new sigma draw. 

##Power simulate (version 1)

get_ss<-function(a,b,mu,start_n=2,M=10000,alpha=.05,pwr=.8){
  lst <- 1/rgamma(M,a,b)
  n <- start_n
  repeat{
    df <- n-1
    m <- mean(pt(qt(1-alpha,df),df,mu*sqrt(n/lst),lower.tail=F))
    if(m >= pwr ){break} else{n<-n+1}
  }
  list(Sample_Size=n,Power=m)
}

##Power simulate (Version 2: find ballpark first)

get_power<-function(lst,mu,n,alpha=.05){mean(pt(qt(1-alpha,n-1),n-1,mu*sqrt(n/lst),lower.tail=F))}


get_ss2<-function(a,b,mu,alpha=.05,pwr=.8, bp=1000, M=1000000){
  ##Find the ballpark
  lst <- 1/rgamma(bp,a,b)
  n <- 2
  repeat{
    m <- get_power(lst,mu,n, alpha=alpha)
    if(m >= pwr ){break} else{n<-n+1}
  }
  n_bp<-n
  lst <- 1/rgamma(M,a,b)
  m <- get_power(lst,mu,n, alpha=alpha)
  if(m < pwr){
    n <- n+1
    repeat{
      m <- get_power(lst,mu,n,alpha)
      #print(m)
      if(m >= pwr ){break} else{n<-n+1}
    }
    
  } else{
    if(n==2){m <- get_power(lst,mu,n,alpha)
    return(c(n_bp,n,m))}
    n <- n-1
    while(n>1){
      m <- get_power(lst,mu,n,alpha)
      #print(m)
      if(m < pwr ){break} else{n<-n-1}
    }
    n <- n + 1
  }
  m <- get_power(lst,mu,n,alpha=alpha)
  return(c(n_bp,n,m))
}

get_ss3 <- function(a,b,mu,start_n=2,alpha=.05,pwr=.8, M=1000){
  k <- seq(M)/M - .5/M
  lst <- 1/qgamma(1-k,a,b)
  n <- start_n
  repeat{
    df <- n-1
    m <- mean(pt(qt(1-alpha,df),df,mu*sqrt(n/lst),lower.tail=F))
    if(m >= pwr ){break} else{n<-n+1}
  }
  list(Sample_Size=n,Power=m)
}


alpha_hat.multi <- Vectorize(function(ss1,ss2,xs1,xs2,k,bhat){
  dvec <- 2*bhat+xs1*(ss1-1) + xs2*(ss2-1) 
  sum((ss1+ss2-2)/dvec)/((k/bhat)-2*sum(1/dvec))
},vectorize.args="bhat")

root_eq.multi <-Vectorize(function(bhat,ss1,ss2,xs1,xs2,k){
  alpha<-alpha_hat.multi(ss1=ss1,ss2=ss2,xs1=xs1,xs2=xs2,k=k,bhat=bhat)
  k*(log(bhat)-digamma(alpha))+
    sum(digamma(alpha+(ss1+ss2-2)/2)-log(bhat+((xs1*(ss1-1)+xs2*(ss2-1))/2)))
},vectorize.args="bhat")

Find_alphabeta.multi<-function(ss1,ss2,xs1,xs2,k){
  bhat<-tryCatch(uniroot(root_eq.multi,interval=c(0.01,100000),ss1=ss1,ss2=ss2,xs1=xs1,xs2=xs2,k=k)$root,error=function(e) 0)
  ahat<-alpha_hat.multi(ss1=ss1,ss2=ss2,xs1=xs1,xs2=xs2,k=k,bhat=bhat)
  c(ahat,bhat)
}


##Two_sample function ####

get_ss.two_same <- function(a,b,delta,start_n=2,M=10000,alpha=.05,pwr=.8,ratio=1){
  lst <- 1/rgamma(M,a,b)
  n <- start_n
  repeat{
    n2 <- ceiling(ratio*n)
    df <- n+n2-2
    m <- mean(pt(qt(1-alpha,df),df,delta*sqrt(1/(lst*(1/n+1/n2))),lower.tail=F))
    if(m >= pwr ){break} else{n<-n+1}
  }
  list(Sample_Size=c(n,n2),Power=m)
}

get_power.two_same<-function(lst,delta,n,n2,alpha=.05){mean(pt(qt(1-alpha,n-1),n-1,delta*sqrt(1/(lst*(1/n+1/n2))),lower.tail=F))}

get_ss2.two_same<-function(a,b,delta,alpha=.05,pwr=.8, bp=1000, M=1000000,ratio=1){
  ##Find the ballpark
  lst <- 1/rgamma(bp,a,b)
  n <- 2
  repeat{
    n2 <- ceiling(n*ratio)
    m <- get_power.two_same(lst,delta,n,n2,alpha=alpha)
    if(m >= pwr ){break} else{n<-n+1}
  }
  n_bp<-n
  lst <- 1/rgamma(M,a,b)
  m <- get_power.two_same(lst,delta,n,n2,alpha=alpha)
  if(m < pwr){
    n <- n+1
    repeat{
      n2 <- ceiling(n*ratio)
      m <- get_power.two_same(lst,delta,n,n2,alpha=alpha)
      #print(m)
      if(m >= pwr ){break} else{n<-n+1}
    }
    
  } else{
    n <- n-1
    while(n>1){
      n2 <- ceiling(n*ratio)
      m <- get_power.two_same(lst,delta,n,n2,alpha=alpha)
      #print(m)
      if(m < pwr ){break} else{n<-n-1}
    }
    n <- n + 1
    n2 <- ceiling(n*ratio)
  }
  m <- get_power.two_same(lst,delta,n,n2,alpha=alpha)
  return(c(n_bp,n,n2,m))
}

get_ss3.two_same <- function(a,b,delta,start_n=2,alpha=.05,pwr=.8, M=1000,ratio=1){
  k <- seq(M)/M - .5/M
  lst <- 1/qgamma(1-k,a,b)
  n <- start_n
  repeat{
    n2 <- ceiling(n*ratio)
    df <- n+n2-2
    m <- mean(pt(qt(1-alpha,df),df,delta*sqrt(1/(lst*(1/n+1/n2))),lower.tail=F))
    if(m >= pwr ){break} else{n<-n+1}
  }
  list(Sample_Size=c(n,n2),Power=m)
}
