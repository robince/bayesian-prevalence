# Computation of HDPIs

# If the package "nleqslv" is not installed on your system, then
# uncomment the following line and install this package.

#   install.packages("nleqslv", dependencies =TRUE)


   
#    Load the nonlinear equation solver
   
 library(nleqslv)
         
 #    Set the working directory

 setwd("/users/jk/Desktop/PrevStuff")
         
 #   Load the prev_hdpi function from the current workng directory

 source("bayesprev.R")
   
#   Set the parameters: Example 1

 
k = 7; n = 50; p = 0.96
 
 int1 <-  hdpi(p, k, n)  
 int1


 #   Set the parameters: Example 2

 
 k = 11; n = 50; a = 0.0; p = 0.96
 
 int2 <-  hdpi(p, k, n)
      
 int2


 #   Set the parameters: Example 3

 
 k = 32; n = 50; p = 0.96
 
 int3 <-  hdpi(p, k, n)
      
 int3


 #   Set the parameters: Example 4

 
 k = 35; n = 50; p = 0.96
 
 int4 <-  hdpi(p, k, n)
      
 int4

 #   Calculation of lower bounds

 gc1 <-  bound(0.05, 7, 50)
 gc2 <-  bound(0.05, 11, 50)
 gc3 <-  bound(0.05, 32, 50)
 gc4 <-   bound(0.05, 35, 50)

 c(gc1, gc2, gc3, gc4)



 #   Calculation of MAP estimates

 map1 <-  map(7, 50)
 map2 <-  map(11, 50)
 map3 <-  map(32, 50)
 map4 <-  map(35, 50)
    
 c(map1, map2, map3, map4)
    
#  Plot posterior

xvals <- seq(0, 1, .01)
pvals <- posterior(xvals, 32, 50)
plot(xvals, pvals, type ="l", xlab ="x", ylab ="Posterior density")


                 
