# Bayesian-Meta-SSD
Code containing all functions and algorithms used in the paper "A Meta-analysis based Hierarchical Variance Model for Powering One and Two-sample t-tests"

The file "functions_clean.r" contains several functions pertaining to this paper, however I am highlighting those that are important and can be used for analysis here:

Find_alphabeta(ss,xs,k):  this function outputs MMLE values for the shape and scale parameters of the inverse gamma distribution followed by theta. ss is a vector
                          of sample sizes from related studies, xs is a vector of sample variances from related studies, k is the number of studies. Note that ss and 
                          xs must be of length k
                          
get_ss2(a,b,mu,alpha=.05,pwr=.8, bp=1000, M=1000000): This is the SRS algorithm for detecting the sample size. a/b are the shape and scale parameters for the
                                                                the IG distribution, mu is the unstandardized effect size, alpha is the Type I error rate, pwr is 
                                                                targeted power, bp is the number of draws to determine a "ballpark" estimate of the sample size and M is 
                                                                the number of draws used to determine a more accurate sample size estimate. Note that bp is meant to be 
                                                                a smaller number to increase efficiency, while M should be large to maintain accuracy.
                                                                
get_ss3(a,b,mu,start_n=2,alpha=.05,pwr=.8, M=1000): This is the discretized algorithm for detecting the sample size. a/b are the shape and scale parameters
                                                                for the IG distribution, mu is the unstandardized effect size, start_n is the starting value guess for 
                                                                the sample size (must be a low guess), alpha is the Type I error rate, pwr is targeted power, and M is
                                                                the number of descritized points in the distribution. As mentioned in the paper, this algorithm is very
                                                                accurate even for small M (100-1000) and is generally very fast. We find it unneccesary to have start_n
                                                                above 2.
                                                                
Two sample versions of these functions (Find_alphabeta.multi, get_ss2.multi, get_ss3.multi) exist as well with relatively intuitive extensions. These functions are 
meant to be used when it can be assumed that sample variances from different groups from a specific study have come from the same "true" variance.
                                                               
                          
                          
