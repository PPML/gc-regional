##################################
### calculate population sizes ###
##################################
#takes total popuation size and calculates the number of individuals in each sex, age, subpopulation, and sexual activity group
#i=subpopulation, j=sex, k=sexual activity group, l=age group

pop_calc <- function(n.total){

  # retrieve environment variables
  n.total <- gc_env$n.total
  p.msm <- gc_env$p.msm
  p.s.1 <- gc_env$p.s.1
  p.s.2 <- gc_env$p.s.2
  p.s.3 <- gc_env$p.s.3
  p.sa.m.y <- gc_env$p.sa.m.y
  p.sa.m.o <- gc_env$p.sa.m.o
  p.sa.f.y <- gc_env$p.sa.f.y
  p.sa.f.o<- gc_env$p.sa.f.o
  aging.rate.m <- gc_env$aging.rate.m
  aging.rate.f <- gc_env$aging.rate.f
  p.low <- gc_env$p.low
  p.low.msm <- gc_env$p.low.msm
  p.sa.m <- gc_env$p.sa.m
  p.sa.f <- gc_env$p.sa.f
  i <- gc_env$i
  j <- gc_env$j
  k <- gc_env$k
  l <- gc_env$l
  index <- gc_env$index
  age.cat<- gc_env$age.cat

  p.s <-c(p.s.1,1-(p.s.1+p.s.3), p.s.3) # proportion of population in subpopulations 1, 2, 3
  p.i<- array(0,dim=c(j,k,l,i)) #array of population distribution by subpopulation
  p.i[,,,1] <-p.s[1] #proportion of pop in i=1 (black)
  p.i[,,,3] <-p.s[3] #proportion of pop in i=3 (Hispanic)
  p.i[,,,2] <-p.s[2] #proportion of pop in i=2 (other)
  p.i[1,,,4] <-p.msm #proportion of pop in i=4 (MSM)
  p.l<-array(dim=c(j,k,l,i)) #age distribution
  p.l[,,1,]<-age.cat[1]/sum(age.cat) #proportion of population in l=1 (young age cat)
  p.l[,,2,]<-age.cat[2]/sum(age.cat) #proportion of population in l=2 (old age cat)
  p.j<-array(dim=c(j,k,l,i)) # pop dist by sex (assume equal size populations for male and female)
  p.j[1,,,1:3]<-0.5*(1-p.msm) #male heterosexual population - adjust proportion of heterosexual male population to account for MSM
  p.j[2,,,1:3]<-0.5 #female population
  p.j[1,,,4]<-0.5 #MSM population
  p.j[2,,,4]<-0 # no females in MSM subpop (4)
  p.j.1.1.1 <-rep(p.low,2) #proportion of M and F in AC=1, for i=1, l=1
  p.j.2.2.1 <-rep(p.low,2) #proportion of M and F in AC=1, for i=1 , l=2
  p.j.1.1.2 <-rep(p.low,2) #proportion of M and F in AC=1 for i=2, l=1
  p.j.2.2.2 <-rep(p.low,2) #proportion of M and F in AC=1 for i=2, l=2
  p.j.1.1.3 <-rep(p.low,2) #proportion of M and F in AC=1 for i=3, l=1
  p.j.2.2.3 <-rep(p.low,2) #proportion of M and F in AC=1 for i=3, l=2
  p.j.1.1.4 <-c(p.low.msm, 0) #proportion of M in AC=1 for i=4 (MSM), l=1
  p.j.2.2.4 <-c(p.low.msm,0)  #proportion of M in AC=1 for i=4, l=2
  p.k <- array(c(p.j.1.1.1,  #distribution of population by low and high activity status
                 1-p.j.1.1.1,
                 p.j.2.2.1,
                 1-p.j.2.2.1,
                 p.j.1.1.2,
                 1-p.j.1.1.2,
                 p.j.2.2.2,
                 1-p.j.2.2.2,
                 p.j.1.1.3,
                 1-p.j.1.1.3,
                 p.j.2.2.3,
                 1-p.j.2.2.3,
                 p.j.1.1.4,
                 1-p.j.1.1.4,
                 p.j.2.2.4,
                 1-p.j.2.2.4),
               dim=c(j,k,l,i)) #for each i,l, j*k matrix of subpop dist'n
  p.k[2,,,4]<-0  #no females in subpop 4 (MSM)
  gc_assign(p.k) #used by aging_fun

  p.sa.array<-array(dim=c(j,k,l,i))  #proportion of population sexually active
  p.sa.array[1,,1,]<-rep(p.sa.m.y, each=2)
  p.sa.array[1,,2,]<-rep(p.sa.m.o, each=2)
  p.sa.array[2,,1,]<-rep(p.sa.f.y, each=2)
  p.sa.array[2,,2,]<-rep(p.sa.f.o, each=2)
  gc_assign(p.sa.array) #used by aging_fun

  p.dist.global<-p.i*p.j*p.k*p.l #population distribution across subpopulations
  n.dist<-n.total*p.dist.global #population size
  gc_assign(n.dist)
  n.dist.sa<-n.dist*p.sa.array #population size for sexually active
  gc_assign(n.dist.sa)
  p.dist <- sweep(p.dist.global,c(3,4),apply(p.dist.global,c(3,4),sum),"/") #population distribution within subpopulation and age group
  gc_assign(p.dist)
  p.s.dist <- array(0,dim=c(rep((i),2)))   #relative sizes of different subpopulations - assume that mixing outside of subpopulation is proportionate to number of individuals in each subpopulation, EXCLUDING MSM HERE
  p.s.dist[1,2:3] <-p.s[2:3]/sum(p.s[2:3])
  p.s.dist[2,c(1,3)] <- p.s[c(1,3)]/sum(p.s[c(1,3)])
  p.s.dist[3,1:2] <- p.s[1:2]/sum(p.s[1:2])
  p.s.dist[4,1:3] <- p.s[1:3]/sum(p.s[1:3])
  gc_assign(p.s.dist)
  n.s.dist <- apply(n.dist,c(1,4),sum)  # population sizes by sex (j) and subpop (i)
  gc_assign(n.s.dist)
  n.s.dist.sa <- apply(n.dist,c(1,4),sum)  # population sizes by sex (j) and subpop (i) for sexually active population
  gc_assign(n.s.dist.sa)
  n.s.pop <- c(n.s.dist[1,1:3],n.s.dist[2,1:3]) # pop sizes, minus MSM
  gc_assign(n.s.pop)
  n.s.pop.sa <- c(n.s.dist.sa[1,1:3],n.s.dist.sa[2,1:3]) # pop sizes, minus MSM, for sexually active population
  gc_assign(n.s.pop.sa)
  n.i<-c(n.dist[1,,,],n.dist[2,,,]) #(pop size, by sex, subpop, AC, age) M, subpop=1: age=AC=1, age=1 AC=2,age=2 AC=1, age=2 AC=2; M subpop=2: age=AC=1,etc. Then F, subpop1...
  gc_assign(n.i)
  n.sa<-c(n.dist.sa[1,,,],n.dist.sa[2,,,]) # size of sexually active population
  gc_assign(n.sa)
  n.nsa<- n.i - n.sa # size of not sexually active population
  gc_assign(n.nsa)
} #calculate population distribution
