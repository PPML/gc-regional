#############################################
### functions to calculate mixing matrix ####
#############################################

# minimun rates of partner change
# generate an array with minimum rates for each sex and age group (MSM separate
# from other males), used for calculating contact matrix
# @param c.min.params vector of min partnership rates, 1=M age 1, 2=M age 2, 3=F
#   age 1, 4=F age 2, 5= MSM age 1, 6=MSM age 2
c_min_fun <- function(c.min.params) {
  j <- gc_env$j
  k <- gc_env$k
  l <- gc_env$l
  i <- gc_env$i

  #mean rate  of partner change, M and F, age cat 1
  c.min.1<-matrix(c(c.min.params[1], c.min.params[3]),nrow=j,ncol=k)

  #mean rate of partner change, M and F, age cat 2
  c.min.2<-matrix(c(c.min.params[2], c.min.params[4]),nrow=j,ncol=k)

  #mean rate of partner change, MSM, age cat 1
  c.min.3<-matrix(c.min.params[5],nrow=j,ncol=k)

  #mean rate of partner change, MSM, age cat 2
  c.min.4<-matrix(c.min.params[6],nrow=j,ncol=k)

  #array: 1=sex, 2=AC, 3=age, 4=subpop
  c.min <- array(c(rep(c(c.min.1, c.min.2),(i-1)),c.min.3, c.min.4),dim=c(k,j,l,i))

  c.min
}

#generate arrays of age assortativity params to use in contact matrix calculations
#p.all = vector of age assortative mixing params, for each subpop, 2 ages per subpop
pi_fun <- function(pi.all) {
  j <- gc_env$j
  l <- gc_env$l
  i <- gc_env$i
  pi.1.m<-diag(pi.all[1:2],l) #age assortativity, black M
  pi.2.m<-diag(pi.all[3:4],l) #age assortativity, other M
  pi.3.m<-diag(pi.all[5:6],l) #age assortativity, Hispanic M
  pi.4.m<-diag(pi.all[7:7],l) #age assortativity, MSM
  pi.1.f<-diag(pi.all[2:1],l) #age assortativity, black F
  pi.2.f<-diag(pi.all[4:3],l) #age assortativity, other F
  pi.3.f<-diag(pi.all[6:5],l) #age assortativity, Hispanic F
  pi.4.f<-diag(c(0,0),l)  #no females in subpop 4
  pi.m <- c(pi.all,pi.all[7])
  pi.f <- c(pi.all[2:1],pi.all[4:3],pi.all[6:5],0,0)
  pi.all.m<-array(c(pi.1.m,pi.2.m,pi.3.m,pi.4.m),dim=c(j,l,i))
  pi.all.f<-array(c(pi.1.f,pi.2.f,pi.3.f,pi.4.f),dim=c(j,l,i))
  inv.pi.all.m<-aperm(array(t(matrix(rep((1-pi.m),2),ncol=j)),dim=c(l,j,i)),c(2,1,3))
  inv.pi.all.f<-aperm(array(t(matrix(rep((1-pi.f),2),ncol=j)),dim=c(l,j,i)),c(2,1,3))
  list(pi.all.m, pi.all.f, inv.pi.all.m, inv.pi.all.f)
}

####################################
### calculates the mixing matrix ###
####################################
# this function uses the approach of Garnett et al. to balance partnerships
# it starts by calculating the desired number of partners for each subgroup
# and balances them to ensure the total number of partnerships in males = total number in females
# for within subpop partnerships, assume omega = 0.5 (average between males and females)
# across subpops, smaller subpopulation decides
mixing <- function(epsilon,pi.all,theta, c.min, rp.all){

  n.dist.sa <- gc_env$n.dist.sa
  p.s.dist <- gc_env$p.s.dist
  omega <- gc_env$omega
  omega.t <- gc_env$omega.t
  i <- gc_env$i
  j <- gc_env$j
  l <- gc_env$l
  k <- gc_env$k

  #### define arrays and grids for mixing calculations ###
  #var1 = AC of opposite sex, var2= AC of individual, var3= age of opposite sex,
  #var4=age of individual, var5=sex, var6=subpopulation  --> used as index for
  #some calculations
  grid<- as.matrix(expand.grid(1:k,1:k,1:l, 1:l,1:j,1:i))

  # x7=subpop(i) of contact; x6=subpop(i); x5=sex(j);x4=age of contact; x3=age; x2=AC of contact ; x1=AC
  grid.btwn<- as.matrix(expand.grid(1:k,1:k,1:l,1:l,1:j,1:i,1:i))
  colnames(grid)<-c("k","k'","l","l'","j","i")
  colnames(grid.btwn)<-c("k","k'","l","l'","j","i","i'")

  # contact distribution, probablity of contact between individual of activity
  # class k with opposite sex person in AC k, by sex (j) and subpop (i)
  p.c.ac <- array(dim=c(k,k,j,i))

  # contact distribution, probablity of mixing between age groups , by sex (j) and subpop (i)
  p.c.age <- array(dim=c(l,l,j,i))

  # check for balance between F AC=k, male AC=k, F age=l, M age =l, subpop=i
  rb <- array(dim=c(k,k,l,l,i))

  # check if partnerships are balanced between subpopulations
  rb.btwn <- array(0,dim=c(k,k,l,l,i,i))

  #re-checking after balancing within subpop, should all=1
  rb.bal <- array(dim=c(k,k,l,l,i))

  #re-checking between subpop partnerships after balancing, should all=1
  rb.bal.btwn <- array(0,dim=c(k,k,l,l,i,i))

  ### calculate partership acquistion rates and total partnerships ###

  #relative rates of partner acquistion by activity group, subpop 1, age group
  #1, M and F  (from nsfg lifetime partners -nsfg_fitted.distributions.xlsx)
  rp.1 <- matrix(c(rp.all[1:4]),ncol=2, byrow=T)
  rp.2 <- matrix(c(rp.all[5:8]),ncol=2, byrow=T)  #subpop 1, age cat 2
  rp.3 <- matrix(c(rp.all[9:12]),ncol=2, byrow=T)  #subpop 2, age cat 1
  rp.4 <- matrix(c(rp.all[13:16]),ncol=2, byrow=T)   #subpop 2, age cat 2
  rp.5 <- matrix(c(rp.all[17:20]),ncol=2, byrow=T) #subpop 3, age cat 1
  rp.6 <- matrix(c(rp.all[21:24]),ncol=2, byrow=T)  #subpop 3, age cat 2
  rp.7 <- matrix(c(rp.all[25:28]),ncol=2, byrow=T)  #subpop 4, age cat 1
  rp.8 <- matrix(c(rp.all[29:32]),ncol=2, byrow=T)  #subpop 4,age cat2

  # 1=sex, 2=AC, 3=age, 4=subpop
  rp <- array(as.numeric(c(rp.1,rp.2,rp.3,rp.4,rp.5,rp.6,rp.7,rp.8)),dim=c(j,k,l,i))

  c.r <- c.min*rp #average number of contacts by sex(j),AC(k),age(l),subpop(i)
  c.r[2,,,4]<-0  #set population of females in MSM subpop to 0

  #################################################################
  #################################################################

  # total partnerships formed, by subpop(i) and sex (j) - only counting sexually active
  p.tot <- apply (c.r*n.dist.sa, MARGIN=c(4,1), sum)

  #total number of partnerships formed by females, replicated by number of AC, by subpop
  p.tot.f <- matrix(rep(p.tot[,2],k), nrow=k, byrow=T)

  #total number of partnerships formed by males, replicated by number of AC, by subpop
  p.tot.m <- matrix(rep(p.tot[,1],k), nrow=k, byrow=T)

  ### Contact probabilities for mixing by AC ###

  # proportion of contacts allocated to assortative mixing, based on subpop (by AC of M (row) and AC of F (column))
  p.assort.ac <- array(sapply(1:i,function(i) diag(rep(epsilon[i],k)),simplify='array'),dim=c(j,k,i))

  #proportionate contacts = 1-epsilon, setting up proportion of contacts allocated to proportionate mixing, based on AC
  p.prop.ac <- array(sapply(1:(i*l), function(i) matrix(1-epsilon[i],j,k),simplify="array"),dim=c(j,k,i))

  #proportion of partnerships available by males in different strata
  q<-array(rep(as.vector(apply(c.r[1,,,]*n.dist.sa[1,,,],c(1,3),sum)/p.tot.m),each=2),dim=c(j,k,i))

  p.prop.f<-q*p.prop.ac  #proportionate contacts, taking into account distribution of available partnerhsips

  #proportion of partnerships available by females in different strata
  q1<-array(rep(as.vector(apply(c.r[2,,,]*n.dist.sa[2,,,],c(1,3),sum)/p.tot.f),each=2),dim=c(j,k,i))


  q1[,,4]<-q[,,4]  #for MSM, calculate distribution of male partners, same calculation as used for females
  p.prop.m<-q1*p.prop.ac #proportionate contacts, taking into account distribution of available partnerships

  p.m.ac <- p.prop.m + p.assort.ac  #distribution of contacts for males, by AC and subpop
  p.f.ac <- p.prop.f + p.assort.ac # distribution of contacts for females, by AC and subpop
  p.c.ac[,,1,1:i]<-p.m.ac[,,1:i] #(k,k,j,l,i) #p contact between person of AC k and sex j and AC k' in subpop i
  p.c.ac[,,2,1:i]<-p.f.ac[,,1:i]

  ### Contact probabilities for mixing by age group ###
  pi.all.fun<-pi_fun(pi.all)
  p.assort.age.m <- pi.all.fun[[1]] #proportion of contacts allocated for mixing with same age group, M
  p.assort.age.f<-pi.all.fun[[2]] #proportion of contacts allocated for mixing with same age group, F
  p.prop.age.m <- pi.all.fun[[3]] #proportionate mixing, M
  p.prop.age.f <- pi.all.fun[[4]] #proportionate mixing, F

  #proportionate distribution of partnerships
  q<-array(rep(as.vector(apply(c.r[1,,,]*n.dist.sa[1,,,],c(2,3),sum)/p.tot.m),each=2),dim=c(j,l,i))

  p.prop.f.age<-q*p.prop.age.f #proportion of contacts for proportionate mixing, F

  #proportionate distribution of partnerships
  q1<-array(rep(as.vector(apply(c.r[2,,,]*n.dist.sa[2,,,],c(2,3),sum)/p.tot.f),each=2),dim=c(j,l,i))
  q1[,,4]<-q[,,4]  #for MSM, calculate partnership distribution among males
  p.prop.m.age<-q1*p.prop.age.m #proportion of contacts for proportionate mixing, M

  #distribution of contacts for males, by AC and subpop
  p.m.age <- p.prop.m.age + p.assort.age.m

  # distribution of contacts for females, by AC and subpop
  p.f.age <- p.prop.f.age + p.assort.age.f

  #(l,l',j,i) #p contact between person of age l and age l' and sex j  in subpop i
  p.c.age[,,1,1:i]<-p.m.age[,,1:i]
  p.c.age[,,2,1:i]<-p.f.age[,,1:i]

  ### Contact probabilities for mixing by AC and age ###
  #contact matrix for p.c.age*p.c.ac, for k,k'(ac of partner),l,l'(age of
  #partner),j,i; p.c should add up to 1 for a person of AC k, and age l (i.e.,
  #summed over k',l')
  p.c<-array(dim=c(k,k,l,l,j,i))
  for (x in 1:i){    #for given subpopulation
    for (s in 1:j) { #for given sex
      p.c[,,1,1,s,x]<-p.c.ac[,,s,x]*p.c.age[1,1,s,x]
      p.c[,,2,1,s,x]<-p.c.ac[,,s,x]*p.c.age[2,1,s,x]
      p.c[,,1,2,s,x]<-p.c.ac[,,s,x]*p.c.age[1,2,s,x]
      p.c[,,2,2,s,x]<-p.c.ac[,,s,x]*p.c.age[2,2,s,x]
    }
  }


  ### Total number of partnerships within and across subpopulations ###
  p.total <-
    array(aaply(grid, 1, function(x)
      #total partnerships within each stratum x6:subpop(i); x5=sex(j);x4=age of
      #contact; x3=age; x2=AC of contact ; x1=AC
      theta[x[6], x[5]] * p.c[x[1], x[2], x[3], x[4], x[5], x[6]] * c.r[x[5],
      x[1], x[3], x[6]] * n.dist.sa[x[5], x[1], x[3], x[6]]),
      dim = c(k, k, l, l, j, i))
  p.total.btwn <-
    array(aaply(grid, 1, function(x) {
      #total partnerships with other subpops (will be divided between other 2
      #subpops) x6:subpop(i); x5=sex(j);x4=age of contact; x3=age; x2=AC of
      #contact ; x1=AC
    if(x[6]!=4) (
      if(x[5]==1)
        (1-theta[x[6],x[5]]) * p.c[x[1],x[2],x[3],x[4],x[5],x[6]] *
        c.r[x[5],x[1],x[3],x[6]] * n.dist.sa[x[5],x[1],x[3],x[6]]
      else
        #minus MSM partnerships for F
        (1-theta[x[6],x[5]]) *theta[4,1] * p.c[x[1],x[2],x[3],x[4],x[5],x[6]] *
        c.r[x[5],x[1],x[3],x[6]] * n.dist.sa[x[5],x[1],x[3],x[6]] ) 
    else (0) }),
    dim=c(k,k,l,l,j,i))

  p.total.msm <- array(aaply(grid,1, function(x){   
      ### Total partnerships between MSM and other subpopulations (i=1-3 is desired
      ### partnerships for F with MSM, i=4 is desired partnerships with F for MSM)
    if(theta[4,1]<1)(
      if(x[5]==2)   #assign proportion of btwn partnerships in females to be with MSM subpop
        (1-theta[x[6],x[5]]) * (1-theta[4,1]) *
        p.c[x[1],x[2],x[3],x[4],x[5],x[6]] * c.r[x[5],x[1],x[3],x[6]] *
        n.dist.sa[x[5],x[1],x[3],x[6]] #total partnerships with other subpops
        # (will be divided between other 2 subpops) x6:subpop(i);
        # x5=sex(j);x4=age of contact; x3=age; x2=AC of contact ; x1=AC
      else if(x[6]==4)
        (1-theta[x[6],x[5]]) * p.c[x[1],x[2],x[3],x[4],x[5],x[6]] *
        c.r[x[5],x[1],x[3],x[6]] * n.dist.sa[x[5],x[1],x[3],x[6]]
      else(0))
    else (0)}),
    dim=c(k,k,l,l,j,i))

  #round(apply((p.total+ p.total.btwn+p.total.msm), c(6,5),sum))  #check that total partnerships = p.tot

    #ratio of partnerships for given AC (F, M), age (F,M), and subpop
  rb[,,1,1,] <- p.total[,,1,1,2,]/aperm(p.total[,,1,1,1,],c(2,1,3)) 
  rb[,,1,2,] <- p.total[,,1,2,2,]/aperm(p.total[,,2,1,1,],c(2,1,3))
  rb[,,2,1,] <- p.total[,,2,1,2,]/aperm(p.total[,,1,2,1,],c(2,1,3))
  rb[,,2,2,] <- p.total[,,2,2,2,]/aperm(p.total[,,2,2,1,],c(2,1,3))

    #for MSM: ratio of partnerships for given AC (M, M), age (M,M), and subpop
  rb[,,1,1,4] <-  p.total[,,1,1,1,4]/aperm(p.total[,,1,1,1,4],c(2,1))  
  rb[,,1,2,4] <-  p.total[,,1,2,1,4]/aperm(p.total[,,1,2,1,4],c(2,1))
  rb[,,2,1,4] <-  p.total[,,2,1,1,4]/aperm(p.total[,,2,1,1,4],c(2,1))
  rb[,,2,2,4] <-  p.total[,,2,2,1,4]/aperm(p.total[,,2,2,1,4],c(2,1))

  rb[is.na(rb)] <-0 #to prevent errors if epsilon =1
  
    # partnerships for given AC(F,M) and age(F,M) for F of subpop 1 (w M of subpop 2)
  rb.btwn[,,1,1,1,2] <- if (theta[1,2]<1) (p.s.dist[1,2]*p.total.btwn[,,1,1,2,1])/(p.s.dist[2,1]*t(p.total.btwn[,,1,1,1,2])) else (0) 
  rb.btwn[,,1,2,1,2] <- if (theta[1,2]<1) (p.s.dist[1,2]*p.total.btwn[,,1,2,2,1])/(p.s.dist[2,1]*t(p.total.btwn[,,2,1,1,2])) else (0)
  rb.btwn[,,2,1,1,2] <- if (theta[1,2]<1) (p.s.dist[1,2]*p.total.btwn[,,2,1,2,1])/(p.s.dist[2,1]*t(p.total.btwn[,,1,2,1,2])) else (0)
  rb.btwn[,,2,2,1,2] <- if (theta[1,2]<1) (p.s.dist[1,2]*p.total.btwn[,,2,2,2,1])/(p.s.dist[2,1]*t(p.total.btwn[,,2,2,1,2])) else (0)

    # partnerships for given AC(F,M) and age(F,M) for F of subpop 1 (w M of subpop 3)
  rb.btwn[,,1,1,1,3] <- if (theta[1,2]<1) (p.s.dist[1,3]*p.total.btwn[,,1,1,2,1])/(p.s.dist[3,1]*t(p.total.btwn[,,1,1,1,3])) else (0) 
  rb.btwn[,,1,2,1,3] <- if (theta[1,2]<1) (p.s.dist[1,3]*p.total.btwn[,,1,2,2,1])/(p.s.dist[3,1]*t(p.total.btwn[,,2,1,1,3])) else (0)
  rb.btwn[,,2,1,1,3] <- if (theta[1,2]<1) (p.s.dist[1,3]*p.total.btwn[,,2,1,2,1])/(p.s.dist[3,1]*t(p.total.btwn[,,1,2,1,3])) else (0)
  rb.btwn[,,2,2,1,3] <- if (theta[1,2]<1) (p.s.dist[1,3]*p.total.btwn[,,2,2,2,1])/(p.s.dist[3,1]*t(p.total.btwn[,,2,2,1,3])) else (0)

    # partnerships for given AC(F,M) and age(F,M) for F of subpop 1 (w M of subpop 3)
  rb.btwn[,,1,1,1,4] <- if (theta[4,1]<1) (1*p.total.msm[,,1,1,2,1])/(p.s.dist[4,1]*t(p.total.msm[,,1,1,1,4])) else (0) 
  rb.btwn[,,1,2,1,4] <- if (theta[4,1]<1) (1*p.total.msm[,,1,2,2,1])/(p.s.dist[4,1]*t(p.total.msm[,,2,1,1,4])) else (0)
  rb.btwn[,,2,1,1,4] <- if (theta[4,1]<1) (1*p.total.msm[,,2,1,2,1])/(p.s.dist[4,1]*t(p.total.msm[,,1,2,1,4])) else (0)
  rb.btwn[,,2,2,1,4] <- if (theta[4,1]<1) (1*p.total.msm[,,2,2,2,1])/(p.s.dist[4,1]*t(p.total.msm[,,2,2,1,4])) else (0)

    # partnerships for given AC(F,M) and age(F,M) for F of subpop 1 (w M of subpop 3)
  rb.btwn[,,1,1,2,4] <- if (theta[4,1]<1) (1*p.total.msm[,,1,1,2,2])/(p.s.dist[4,2]*t(p.total.msm[,,1,1,1,4])) else (0) 
  rb.btwn[,,1,2,2,4] <- if (theta[4,1]<1) (1*p.total.msm[,,1,2,2,2])/(p.s.dist[4,2]*t(p.total.msm[,,2,1,1,4])) else (0)
  rb.btwn[,,2,1,2,4] <- if (theta[4,1]<1) (1*p.total.msm[,,2,1,2,2])/(p.s.dist[4,2]*t(p.total.msm[,,1,2,1,4])) else (0)
  rb.btwn[,,2,2,2,4] <- if (theta[4,1]<1) (1*p.total.msm[,,2,2,2,2])/(p.s.dist[4,2]*t(p.total.msm[,,2,2,1,4])) else (0)

    # partnerships for given AC(F,M) and age(F,M) for F of subpop 1 (w M of subpop 3)
  rb.btwn[,,1,1,3,4] <- if (theta[4,1]<1) (1*p.total.msm[,,1,1,2,3])/(p.s.dist[4,3]*t(p.total.msm[,,1,1,1,4])) else (0) 
  rb.btwn[,,1,2,3,4] <- if (theta[4,1]<1) (1*p.total.msm[,,1,2,2,3])/(p.s.dist[4,3]*t(p.total.msm[,,2,1,1,4])) else (0)
  rb.btwn[,,2,1,3,4] <- if (theta[4,1]<1) (1*p.total.msm[,,2,1,2,3])/(p.s.dist[4,3]*t(p.total.msm[,,1,2,1,4])) else (0)
  rb.btwn[,,2,2,3,4] <- if (theta[4,1]<1) (1*p.total.msm[,,2,2,2,3])/(p.s.dist[4,3]*t(p.total.msm[,,2,2,1,4])) else (0)

    # partnerships for given AC(F,M) and age(F,M) for F of subpop 2 (w M of subpop 1)
  rb.btwn[,,1,1,2,1] <- if (theta[2,2]<1) (p.s.dist[2,1]*p.total.btwn[,,1,1,2,2])/(p.s.dist[1,2]*t(p.total.btwn[,,1,1,1,1])) else (0) 
  rb.btwn[,,1,2,2,1] <- if (theta[2,2]<1) (p.s.dist[2,1]*p.total.btwn[,,1,2,2,2])/(p.s.dist[1,2]*t(p.total.btwn[,,2,1,1,1])) else (0)
  rb.btwn[,,2,1,2,1] <- if (theta[2,2]<1) (p.s.dist[2,1]*p.total.btwn[,,2,1,2,2])/(p.s.dist[1,2]*t(p.total.btwn[,,1,2,1,1])) else (0)
  rb.btwn[,,2,2,2,1] <- if (theta[2,2]<1) (p.s.dist[2,1]*p.total.btwn[,,2,2,2,2])/(p.s.dist[1,2]*t(p.total.btwn[,,2,2,1,1])) else (0)

    # partnerships for given AC(F,M) and age(F,M) for F of subpop 2 (w M of subpop 3)
  rb.btwn[,,1,1,2,3] <- if (theta[2,2]<1) (p.s.dist[2,3]*p.total.btwn[,,1,1,2,2])/(p.s.dist[3,2]*t(p.total.btwn[,,1,1,1,3])) else (0) 
  rb.btwn[,,1,2,2,3] <- if (theta[2,2]<1) (p.s.dist[2,3]*p.total.btwn[,,1,2,2,2])/(p.s.dist[3,2]*t(p.total.btwn[,,2,1,1,3])) else (0)
  rb.btwn[,,2,1,2,3] <- if (theta[2,2]<1) (p.s.dist[2,3]*p.total.btwn[,,2,1,2,2])/(p.s.dist[3,2]*t(p.total.btwn[,,1,2,1,3])) else (0)
  rb.btwn[,,2,2,2,3] <- if (theta[2,2]<1) (p.s.dist[2,3]*p.total.btwn[,,2,2,2,2])/(p.s.dist[3,2]*t(p.total.btwn[,,2,2,1,3])) else (0)

    # partnerships for given AC(F,M) and age(F,M) for F of subpop 3 (w M of subpop 1)
  rb.btwn[,,1,1,3,1] <- if (theta[3,2]<1) (p.s.dist[3,1]*p.total.btwn[,,1,1,2,3])/(p.s.dist[1,3]*t(p.total.btwn[,,1,1,1,1])) else (0) 
  rb.btwn[,,1,2,3,1] <- if (theta[3,2]<1) (p.s.dist[3,1]*p.total.btwn[,,1,2,2,3])/(p.s.dist[1,3]*t(p.total.btwn[,,2,1,1,1])) else (0)
  rb.btwn[,,2,1,3,1] <- if (theta[3,2]<1) (p.s.dist[3,1]*p.total.btwn[,,2,1,2,3])/(p.s.dist[1,3]*t(p.total.btwn[,,1,2,1,1])) else (0)
  rb.btwn[,,2,2,3,1] <- if (theta[3,2]<1) (p.s.dist[3,1]*p.total.btwn[,,2,2,2,3])/(p.s.dist[1,3]*t(p.total.btwn[,,2,2,1,1])) else (0)

    # partnerships for given AC(F,M) and age(F,M) for F of subpop 3 (w M of subpop 2)
  rb.btwn[,,1,1,3,2] <- if (theta[3,2]<1) (p.s.dist[3,2]*p.total.btwn[,,1,1,2,3])/(p.s.dist[2,3]*t(p.total.btwn[,,1,1,1,2])) else (0) 
  rb.btwn[,,1,2,3,2] <- if (theta[3,2]<1) (p.s.dist[3,2]*p.total.btwn[,,1,2,2,3])/(p.s.dist[2,3]*t(p.total.btwn[,,2,1,1,2])) else (0)
  rb.btwn[,,2,1,3,2] <- if (theta[3,2]<1) (p.s.dist[3,2]*p.total.btwn[,,2,1,2,3])/(p.s.dist[2,3]*t(p.total.btwn[,,1,2,1,2])) else (0)
  rb.btwn[,,2,2,3,2] <- if (theta[3,2]<1) (p.s.dist[3,2]*p.total.btwn[,,2,2,2,3])/(p.s.dist[2,3]*t(p.total.btwn[,,2,2,1,2])) else (0)

  rb.btwn[is.na(rb.btwn)] <-0 #to prevent errors if epsilon/theta/pi =1

  c.bal <-array(aaply(grid,1,function(x){  #for person of AC=k, with partner of AC=k', age=l, partner age=l', sex=j, subpop=i (kk'll'ji)
    if(x[6]!=4) ( #if not MSM
      if (x[5]==1) #if non-MSM and M
        c.r[x[5],x[1],x[3],x[6]]*rb[x[2],x[1],x[4],x[3],x[6]]^omega
      else  # if non-MSM and F
        (c.r[x[5],x[1],x[3],x[6]]/rb[x[1],x[2],x[3],x[4],x[6]]^(1-omega))
    )
    else  # if MSM
      (c.r[x[5],x[1],x[3],x[6]]/rb[x[1],x[2],x[3],x[4],x[6]]^(1-omega))
  }), dim=c(k,k,l,l,j,i))

  c.bal[is.na(c.bal)] <-0  #to prevent errors if epsilon =1
  c.bal[is.infinite(c.bal)] <-0

  c.bal.btwn <-array(aaply(grid.btwn,1,function(x) {
    if (theta[1,1]<1)      ##using omega.t to define degree of compromise -- right now smaller subpops determine mixing
      if (x[6]!=x[7]) (
        if (x[5]==1)   ### if sex = 1 (M)
          c.r[x[5],x[1],x[3],x[6]]*rb.btwn[x[2],x[1],x[4],x[3],x[7],x[6]]^omega.t[x[7],x[6]] ### omega.t[subpop F, subpop M]
        else (c.r[x[5],x[1],x[3],x[6]]/rb.btwn[x[1],x[2],x[3],x[4],x[6],x[7]]^(1-omega.t[x[6],x[7]])))
    else (0)
    else (0)
  }), dim=c(k,k,l,l,j,i,i))

  c.bal.btwn[is.na(c.bal.btwn)] <-0  #to prevent errors if epsilon =1
  c.bal.btwn[is.infinite(c.bal.btwn)] <-0


  #total partnerships within each stratum x6:subpop(i); x5=sex(j);x4=age of contact; x3=age; x2=AC of contact ; x1=AC
  p.total.bal <- array(aaply(grid,1, function(x) theta[x[6],x[5]] *
      p.c[x[1],x[2],x[3],x[4],x[5],x[6]] * c.bal[x[1],x[2],x[3],x[4],x[5],x[6]] *
      n.dist.sa[x[5],x[1],x[3],x[6]]), dim=c(k,k,l,l,j,i))
  gc_env$p.total.bal <- p.total.bal


  #total partnerships with other subpops (will be divided between other 2 subpops) x6:subpop(i); x5=sex(j);x4=age of contact; x3=age; x2=AC of contact ; x1=AC
  p.total.bal.btwn <-array(aaply(grid.btwn, 1, function(x) {
    if(x[6]!=4 & x[7]!=4) (
      if(x[5]==1)
        (1-theta[x[6],x[5]]) * p.s.dist[x[6],x[7]] *
        p.c[x[1],x[2],x[3],x[4],x[5],x[6]] *
        c.bal.btwn[x[1],x[2],x[3],x[4],x[5],x[6],x[7]] *
        n.dist.sa[x[5],x[1],x[3],x[6]]
      else
        (1-theta[x[6],x[5]]) *theta[4,1] * p.s.dist[x[6],x[7]] *
        p.c[x[1],x[2],x[3],x[4],x[5],x[6]] *
        c.bal.btwn[x[1],x[2],x[3],x[4],x[5],x[6],x[7]] *
        n.dist.sa[x[5],x[1],x[3],x[6]])
    else (0) }),
    dim=c(k,k,l,l,j,i,i))
  gc_env$p.total.bal.btwn <- p.total.bal.btwn

  p.total.bal.msm <- array(aaply(grid.btwn,1, function(x){   ### Total partnerships between MSM and other subpopulations
    if(theta[4,1]<1)(
      if(x[5]==2 & x[7]==4) (  #assign proportion of btwn partnerships in females to be with MSM subpop
        #total partnerships with other subpops (will be divided between other 2 subpops) x6:subpop(i); x5=sex(j);x4=age of contact; x3=age; x2=AC of contact ; x1=AC
        (1-theta[x[6],x[5]]) * (1-theta[4,1]) * 1 *
        p.c[x[1],x[2],x[3],x[4],x[5],x[6]] *
        c.bal.btwn[x[1],x[2],x[3],x[4],x[5],x[6],x[7]] *
        n.dist.sa[x[5],x[1],x[3],x[6]])
      else if(x[5]==1 & x[6]==4)
        (1-theta[4,1]) *p.s.dist[x[6],x[7]] *
        p.c[x[1],x[2],x[3],x[4],x[5],x[6]] *
        c.bal.btwn[x[1],x[2],x[3],x[4],x[5],x[6],x[7]] *
        n.dist.sa[x[5],x[1],x[3],x[6]]
      else(0))
    else(0)}),
    dim=c(k,k,l,l,j,i,i))
  gc_env$p.total.bal.msm <- p.total.bal.msm

  ### can run code below to check that rb.bal and rb.bal.btwn  =1 (i.e., partners are indeed balanced)  ##
  # rb.bal[,,1,1,] <- p.total.bal[,,1,1,2,]/aperm(p.total.bal[,,1,1,1,],c(2,1,3)) #ratio of partnerships for given AC (F, M), age (F,M), and subpop
  # rb.bal[,,1,2,] <- p.total.bal[,,1,2,2,]/aperm(p.total.bal[,,2,1,1,],c(2,1,3))
  # rb.bal[,,2,1,] <- p.total.bal[,,2,1,2,]/aperm(p.total.bal[,,1,2,1,],c(2,1,3))
  # rb.bal[,,2,2,] <- p.total.bal[,,2,2,2,]/aperm(p.total.bal[,,2,2,1,],c(2,1,3))
  # rb.bal[,,1,1,4] <- p.total.bal[,,1,1,1,4]/aperm(p.total.bal[,,1,1,1,4],c(2,1)) #for MSM: ratio of partnerships for given AC (M, M), age (M,M), and subpop
  # rb.bal[,,1,2,4] <- p.total.bal[,,1,2,1,4]/aperm(p.total.bal[,,1,2,1,4],c(2,1))
  # rb.bal[,,2,1,4] <- p.total.bal[,,2,1,1,4]/aperm(p.total.bal[,,2,1,1,4],c(2,1))
  # rb.bal[,,2,2,4] <- p.total.bal[,,2,2,1,4]/aperm(p.total.bal[,,2,2,1,4],c(2,1))
  #
  # rb.bal[is.na(rb)]
  #
  # rb.bal.btwn[,,1,1,1,2] <- if (theta[2,1]<1) (p.total.bal.btwn[,,1,1,2,1,2])/(t(p.total.bal.btwn[,,1,1,1,2,1])) else (0) # partnerships for given AC(F,M) and age(F,M) for F of subpop 1 (w M of subpop 2)
  # rb.bal.btwn[,,1,2,1,2] <- if (theta[2,1]<1) (p.total.bal.btwn[,,1,2,2,1,2])/(t(p.total.bal.btwn[,,2,1,1,2,1])) else (0)
  # rb.bal.btwn[,,2,1,1,2] <- if (theta[2,1]<1) (p.total.bal.btwn[,,2,1,2,1,2])/(t(p.total.bal.btwn[,,1,2,1,2,1])) else (0)
  # rb.bal.btwn[,,2,2,1,2] <- if (theta[2,1]<1) (p.total.bal.btwn[,,2,2,2,1,2])/(t(p.total.bal.btwn[,,2,2,1,2,1])) else (0)
  #
  # rb.bal.btwn[,,1,1,1,3] <- if (theta[2,1]<1) (p.total.bal.btwn[,,1,1,2,1,3])/(t(p.total.bal.btwn[,,1,1,1,3,1])) else (0) # partnerships for given AC(F,M) and age(F,M) for F of subpop 1 (w M of subpop 3)
  # rb.bal.btwn[,,1,2,1,3] <- if (theta[2,1]<1) (p.total.bal.btwn[,,1,2,2,1,3])/(t(p.total.bal.btwn[,,2,1,1,3,1])) else (0)
  # rb.bal.btwn[,,2,1,1,3] <- if (theta[2,1]<1) (p.total.bal.btwn[,,2,1,2,1,3])/(t(p.total.bal.btwn[,,1,2,1,3,1])) else (0)
  # rb.bal.btwn[,,2,2,1,3] <- if (theta[2,1]<1) (p.total.bal.btwn[,,2,2,2,1,3])/(t(p.total.bal.btwn[,,2,2,1,3,1])) else (0)
  #
  # rb.bal.btwn[,,1,1,2,1] <- if (theta[2,2]<1) (p.total.bal.btwn[,,1,1,2,2,1])/(t(p.total.bal.btwn[,,1,1,1,1,2])) else (0) # partnerships for given AC(F,M) and age(F,M) for F of subpop 2 (w M of subpop 1)
  # rb.bal.btwn[,,1,2,2,1] <- if (theta[2,2]<1) (p.total.bal.btwn[,,1,2,2,2,1])/(t(p.total.bal.btwn[,,2,1,1,1,2])) else (0)
  # rb.bal.btwn[,,2,1,2,1] <- if (theta[2,2]<1) (p.total.bal.btwn[,,2,1,2,2,1])/(t(p.total.bal.btwn[,,1,2,1,1,2])) else (0)
  # rb.bal.btwn[,,2,2,2,1] <- if (theta[2,2]<1) (p.total.bal.btwn[,,2,2,2,2,1])/(t(p.total.bal.btwn[,,2,2,1,1,2])) else (0)
  #
  # rb.bal.btwn[,,1,1,2,3] <- if (theta[2,2]<1) (p.total.bal.btwn[,,1,1,2,2,3])/(t(p.total.bal.btwn[,,1,1,1,3,2])) else (0) # partnerships for given AC(F,M) and age(F,M) for F of subpop 2 (w M of subpop 3)
  # rb.bal.btwn[,,1,2,2,3] <- if (theta[2,2]<1) (p.total.bal.btwn[,,1,2,2,2,3])/(t(p.total.bal.btwn[,,2,1,1,3,2])) else (0)
  # rb.bal.btwn[,,2,1,2,3] <- if (theta[2,2]<1) (p.total.bal.btwn[,,2,1,2,2,3])/(t(p.total.bal.btwn[,,1,2,1,3,2])) else (0)
  # rb.bal.btwn[,,2,2,2,3] <- if (theta[2,2]<1) (p.total.bal.btwn[,,2,2,2,2,3])/(t(p.total.bal.btwn[,,2,2,1,3,2])) else (0)
  #
  # rb.bal.btwn[,,1,1,3,1] <- if (theta[3,2]<1) (p.total.bal.btwn[,,1,1,2,3,1])/(t(p.total.bal.btwn[,,1,1,1,1,3])) else (0) # partnerships for given AC(F,M) and age(F,M) for F of subpop 3 (w M of subpop 1)
  # rb.bal.btwn[,,1,2,3,1] <- if (theta[3,2]<1) (p.total.bal.btwn[,,1,2,2,3,1])/(t(p.total.bal.btwn[,,2,1,1,1,3])) else (0)
  # rb.bal.btwn[,,2,1,3,1] <- if (theta[3,2]<1) (p.total.bal.btwn[,,2,1,2,3,1])/(t(p.total.bal.btwn[,,1,2,1,1,3])) else (0)
  # rb.bal.btwn[,,2,2,3,1] <- if (theta[3,2]<1) (p.total.bal.btwn[,,2,2,2,3,1])/(t(p.total.bal.btwn[,,2,2,1,1,3])) else (0)
  #
  # rb.bal.btwn[,,1,1,3,2] <- if (theta[3,2]<1) (p.total.bal.btwn[,,1,1,2,3,2])/(t(p.total.bal.btwn[,,1,1,1,2,3])) else (0) # partnerships for given AC(F,M) and age(F,M) for F of subpop 3 (w M of subpop 2)
  # rb.bal.btwn[,,1,2,3,2] <- if (theta[3,2]<1) (p.total.bal.btwn[,,1,2,2,3,2])/(t(p.total.bal.btwn[,,2,1,1,2,3])) else (0)
  # rb.bal.btwn[,,2,1,3,2] <- if (theta[3,2]<1) (p.total.bal.btwn[,,2,1,2,3,2])/(t(p.total.bal.btwn[,,1,2,1,2,3])) else (0)
  # rb.bal.btwn[,,2,2,3,2] <- if (theta[3,2]<1) (p.total.bal.btwn[,,2,2,2,3,2])/(t(p.total.bal.btwn[,,2,2,1,2,3])) else (0)
  #
  # rb.bal.btwn[,,1,1,1,4] <- if (theta[4,1]<1) (1*p.total.bal.msm[,,1,1,2,1,4])/(t(p.total.bal.msm[,,1,1,1,4,1])) else (0) # partnerships for given AC(F,M) and age(F,M) for F of subpop 1 (w M of subpop 3)
  # rb.bal.btwn[,,1,2,1,4] <- if (theta[4,1]<1) (1*p.total.bal.msm[,,1,2,2,1,4])/(t(p.total.bal.msm[,,2,1,1,4,1])) else (0)
  # rb.bal.btwn[,,2,1,1,4] <- if (theta[4,1]<1) (1*p.total.bal.msm[,,2,1,2,1,4])/(t(p.total.bal.msm[,,1,2,1,4,1])) else (0)
  # rb.bal.btwn[,,2,2,1,4] <- if (theta[4,1]<1) (1*p.total.bal.msm[,,2,2,2,1,4])/(t(p.total.bal.msm[,,2,2,1,4,1])) else (0)
  #
  # rb.bal.btwn[,,1,1,2,4] <- if (theta[4,1]<1) (1*p.total.bal.msm[,,1,1,2,2,4])/(t(p.total.bal.msm[,,1,1,1,4,2])) else (0) # partnerships for given AC(F,M) and age(F,M) for F of subpop 1 (w M of subpop 3)
  # rb.bal.btwn[,,1,2,2,4] <- if (theta[4,1]<1) (1*p.total.bal.msm[,,1,2,2,2,4])/(t(p.total.bal.msm[,,2,1,1,4,2])) else (0)
  # rb.bal.btwn[,,2,1,2,4] <- if (theta[4,1]<1) (1*p.total.bal.msm[,,2,1,2,2,4])/(t(p.total.bal.msm[,,1,2,1,4,2])) else (0)
  # rb.bal.btwn[,,2,2,2,4] <- if (theta[4,1]<1) (1*p.total.bal.msm[,,2,2,2,2,4])/(t(p.total.bal.msm[,,2,2,1,4,2])) else (0)
  #
  # rb.bal.btwn[,,1,1,3,4] <- if (theta[4,1]<1) (1*p.total.bal.msm[,,1,1,2,3,4])/(t(p.total.bal.msm[,,1,1,1,4,3])) else (0) # partnerships for given AC(F,M) and age(F,M) for F of subpop 1 (w M of subpop 3)
  # rb.bal.btwn[,,1,2,3,4] <- if (theta[4,1]<1) (1*p.total.bal.msm[,,1,2,2,3,4])/(t(p.total.bal.msm[,,2,1,1,4,3])) else (0)
  # rb.bal.btwn[,,2,1,3,4] <- if (theta[4,1]<1) (1*p.total.bal.msm[,,2,1,2,3,4])/(t(p.total.bal.msm[,,1,2,1,4,3])) else (0)
  # rb.bal.btwn[,,2,2,3,4] <- if (theta[4,1]<1) (1*p.total.bal.msm[,,2,2,2,3,4])/(t(p.total.bal.msm[,,2,2,1,4,3])) else (0)

  cm <- array(aaply(grid,1, function(x) theta[x[6],x[5]] * p.c[x[1],x[2],x[3],x[4],x[5],x[6]] * c.bal[x[1],x[2],x[3],x[4],x[5],x[6]]), dim=c(k,k,l,l,j,i))

  cm.btwn <-array(aaply(grid.btwn, 1, function(x) { #total heterosexual partnerships with other subpops (will be divided between other 2 subpops) x6:subpop(i); x5=sex(j);x4=age of contact; x3=age; x2=AC of contact ; x1=AC
    if(theta[2,1]<1)(
      if(x[6]!=4 & x[7]!=4) (
        if(x[5]==1)
          (1-theta[x[6],x[5]]) * p.s.dist[x[6],x[7]] *
          p.c[x[1],x[2],x[3],x[4],x[5],x[6]] *
          c.bal.btwn[x[1],x[2],x[3],x[4],x[5],x[6],x[7]]
        else
          (1-theta[x[6],x[5]]) *theta[4,1] * p.s.dist[x[6],x[7]] *
          p.c[x[1],x[2],x[3],x[4],x[5],x[6]] *
          c.bal.btwn[x[1],x[2],x[3],x[4],x[5],x[6],x[7]] ) #excluding F partnerships with MSM
      else (0))
    else(0)}),
    dim=c(k,k,l,l,j,i,i))

  cm.msm <- array(aaply(grid.btwn,1, function(x){   ### Total partnerships between MSM and other subpopulations
    if(theta[4,1]<1) (
      if(x[5]==2 & x[7]==4) (  #assign proportion of btwn partnerships in females to be with MSM subpop
        # total partnerships with other subpops (will be divided between other 2
        # subpops) x6:subpop(i); x5=sex(j);x4=age of contact; x3=age; x2=AC of
        # contact ; x1=AC
        (1-theta[x[6],x[5]]) * (1-theta[4,1]) * 1 *
        p.c[x[1],x[2],x[3],x[4],x[5],x[6]] *
        c.bal.btwn[x[1],x[2],x[3],x[4],x[5],x[6],x[7]] )
      else if(x[5]==1 & x[6]==4)
        (1-theta[4,1]) *p.s.dist[x[6],x[7]] *
        p.c[x[1],x[2],x[3],x[4],x[5],x[6]] *
        c.bal.btwn[x[1],x[2],x[3],x[4],x[5],x[6],x[7]]
      else(0))
    else(0)
  }),
  dim=c(k,k,l,l,j,i,i))

  cm.all <- cm.btwn + cm.msm
  cm.all[,,,,,1,1] <- cm[,,,,,1]
  cm.all[,,,,,2,2] <- cm[,,,,,2]
  cm.all[,,,,,3,3] <- cm[,,,,,3]
  cm.all[,,,,,4,4] <- cm[,,,,,4]

  cm.low.1 <- c(cm.all[,1,,,1,1,],cm.all[,1,,,2,1,]) #contacts with person of opposite sex in low AC (by sex, AC, age, and subpop) for person in subpop1
  cm.high.1 <- c(cm.all[,2,,,1,1,],cm.all[,2,,,2,1,]) #contacts with person of opposite sex in high AC (by sex, AC, age, and subpop) for person in subpop1
  cm.low.2 <- c(cm.all[,1,,,1,2,],cm.all[,1,,,2,2,]) #contacts with person of opposite sex in low AC (by sex, AC, age, and subpop) for person in subpop2
  cm.high.2 <- c(cm.all[,2,,,1,2,],cm.all[,2,,,2,2,]) #contacts with person of opposite sex  in high AC (by sex, AC, age, and subpop) for person in subpop2
  cm.low.3 <- c(cm.all[,1,,,1,3,],cm.all[,1,,,2,3,]) #contacts with person of opposite sex  in low AC (by sex, AC, age, and subpop) for person in subpop3
  cm.high.3 <- c(cm.all[,2,,,1,3,],cm.all[,2,,,2,3,]) #contacts with person of opposite sex in high AC (by sex, AC, age, and subpop) for person in subpop3
  cm.low.4 <- c(cm.all[,1,,,1,4,]) #contacts with person of same sex  in low AC (by sex, AC, age, and subpop) for person in subpop4
  cm.high.4 <- c(cm.all[,2,,,1,4,]) #contacts with person of same sex in high AC (by sex, AC, age, and subpop) for person in subpop4
  cm.list <- list(cm.low1=cm.low.1, cm.low2=cm.low.2, cm.low3=cm.low.3, cm.low4=cm.low.4, cm.high1=cm.high.1, cm.high2=cm.high.2, cm.high3=cm.high.3, cm.high4=cm.high.4)
  gc_env$cm.list <- cm.list

  part.all <- p.total.bal.btwn + p.total.bal.msm #total partners
  part.all[,,,,,1,1] <- p.total.bal[,,,,,1]
  part.all[,,,,,2,2] <- p.total.bal[,,,,,2]
  part.all[,,,,,3,3] <- p.total.bal[,,,,,3]
  part.all[,,,,,4,4] <- p.total.bal[,,,,,4]
  gc_assign(part.all)
  part.all.m <-apply(part.all[,,,,1,,],c(5:6) ,sum)
  part.all.f <-apply(part.all[,,,,2,,],c(5:6) ,sum)
  part.num.age <- as.vector(apply(part.all, c(3,5), sum))
  part.num.sub <-c(apply(part.all.m[1:3,1:3],1,sum), apply(part.all.f[1:3,1:3],1,sum))
  gc_assign(part.all.m)
  gc_assign(part.all.f)
  gc_assign(part.num.age)
  gc_assign(part.num.sub)
}
