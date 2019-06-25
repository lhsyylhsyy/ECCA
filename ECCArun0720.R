
setwd('C:\\Users\\xun02\\Dropbox\\Shao\\SIR') #work

#load("run0720.RData")

library(MASS)
library(dr)
library(R.matlab)
library(lattice)

#start MatLab
matlab <- Matlab()
isOpen <- open(matlab)


#my angle calculation function
angle <- function(x,y)
{
  x <- as.vector(x)/sqrt(sum(x^2))
  y <- as.vector(y)/sqrt(sum(y^2))
  rho <- acos(round(t(x)%*% y, 10))#added an approximation to avoid errors when x and y are in exact opposite direction
  rho
}

#I guess I'll have to write a function
multiSIR <- function(x,y,
                     # rounddown=0, #any coefficient with absolute value below this number is rounded down to 0
                     niters=10,
                     nd=1, #working dimensions of the space
                     x.initial=NULL, y.initial=NULL
)
{
  x <- as.matrix(x)
  y <- as.matrix(y)
  #step 1, CCA
  #step 1, initial value
  if(!is.null(x.initial))
  {
    vec1x <- x.initial
    if(!is.null(y.initial))
    {
      vec1 <- y.initial
    } else
    {
      tmp <- cancor(x %*% vec1x, y)
      vec1 <- tmp$ycoef[,1]
    }
  } else if (!is.null(y.initial))
  {
    vec1 <- y.initial
    tmp <- cancor(x, y%*% vec1)
    vec1x <- tmp$xcoef[,1]
  } else
  {
    step1 <- cancor(x,y)
    vec1 <- step1$ycoef[,1]
    vec1x <- step1$xcoef[,1]
  }
  #standardize
  vec1y.ini <- vec1/sqrt(t(vec1)%*% vec1)
  vec1x.ini <- vec1x/sqrt(t(vec1x)%*% vec1x)
  vec1y <- vec1y.ini
  vec1x <- vec1x.ini
  Yvec <- list(vec1y)
  Xvec <- list(vec1x)
  #records value of previous iteration, initial value set to 0
  vec0y <- 0*vec1y
  vec0x <- 0*vec1x
  #starting point
  ynew <- y %*% vec1y
  cor1 <- rep(NA, niters)
  i <- 1
  while(i <= niters & sum(vec0y!=vec1y)+sum(vec0x != vec1x) >0)
  {
    vec0y <- vec1y
    vec0x <- vec1x
    #use reduced y dimension to find SIR directions
    step2 <- dr(ynew~x, method='sir')
    #keep significant directions
    #reduce x dimension
    vec1x <- step2$evectors[,1:nd]
    vec1x <- as.matrix(vec1x)
    for(j in 1:nd)
    {
      if((angle(vec1x[,j], as.matrix(vec1x.ini)[,j])>pi/2))
        vec1x[,j] <- -vec1x[,j]
    }
    xnew <- x %*% vec1x
    Xvec[[i+1]] <- vec1x
    #use reduced x dimension to find MP directions
    step3 <- dr(xnew~y, method='sir')
    vec1y <-summary(step3)$evectors[,1:nd]
    vec1y <- as.matrix(vec1y)
    cor1[i] <- step3$evalues[1]
    for(j in 1:nd)
    {
      if((angle(vec1y[,j], as.matrix(vec1y.ini)[,j])>pi/2))
        vec1y[,j] <- -vec1y[,j]
    }
    ynew <- y %*% vec1y
    Yvec[[i+1]] <- vec1y
    i <- i+1
  }
  if(sum(vec0y!=vec1y)+sum(vec0x != vec1x) >0) converge <- FALSE else 
    converge <- TRUE
  Xvec <- Reduce(cbind, Xvec)
  Yvec <- Reduce(cbind, Yvec)
  return(list(Xvec=Xvec, Yvec=Yvec,x=x, y=y,
              x.original =x, y.original=y,
              evalue = cor1, converge=converge))
}

#extended CCA
ECCA <- function(x,y,niters=20,
                 x.original =x, y.original=y,
                 standardize=F,
                 x.initial=NULL, y.initial=NULL)
{

  #remove missing value
  keep <- apply(!is.na(cbind(x,y)), 1, prod)==1
  x <- as.matrix(x[keep,])
  y <- as.matrix(y[keep,])
  if(standardize==T)
  {
    #standardize
    x <- apply(x,2, function(x) (x-mean(x))/sd(x))
    y <- apply(y,2, function(x) (x-mean(x))/sd(x))
  }
  
  if(dim(x)[2]==1 & dim(y)[2] >1)
  {
    result <- result1 <- result2 <- dr(x~y,method='sir')
    Xvec <- 1
    Yvec <- result$evectors[,1]
    cor1 <- result$evalues[1]
  } else if(dim(y)[2]==1 & dim(x)[2] >1)
  {
    result <- result1 <- result2 <-  dr(y~x,method='sir')
    Xvec <- result$evectors[,1]
    Yvec <- 1
    cor1 <- result$evalues[1]
  } else if(dim(x)[2]==1 & dim(y)[2] ==1)
    {
    Xvec <- Yvec <- 1
  cor1 <- cor(x,y)
  result1 <- result2 <- NULL
  }  else
  {
    #step 1, initial value
    if(!is.null(x.initial))
    {
      vec1x <- x.initial
      if(!is.null(y.initial))
      {
        vec1 <- y.initial
      } else
      {
        tmp <- cancor(x %*% vec1x, y)
        vec1 <- tmp$ycoef[,1]
      }
    } else if (!is.null(y.initial))
    {
      vec1 <- y.initial
      tmp <- cancor(x, y%*% vec1)
      vec1x <- tmp$xcoef[,1]
    } else
    {
      step1 <- cancor(x,y)
      vec1 <- step1$ycoef[,1]
      vec1x <- step1$xcoef[,1]
    }
    #standardize
    vec1.ini <- vec1/sqrt(t(vec1)%*% vec1)
    vec1x.ini <- vec1x/sqrt(t(vec1x)%*% vec1x)
    #SIR-MAVE
    converge <- NA
    cor1 <- NA
    vec1 <- vec1.ini
    vec1x <- vec1x.ini
    Yvec <- list(vec1)
    Xvec <- list(vec1x)
    #records value of previous iteration, initial value set to 0
    vec0 <- 0*vec1
    vec0x <- 0*vec1x
    #starting point
    ynew <- y %*% vec1
    i <- 1
    while(i <= niters & sum(vec0!=vec1)+sum(vec0x != vec1x) >0)
    {
      vec0 <- vec1
      vec0x <- vec1x
      #use reduced y dimension to find SIR directions
      step2 <- dr(ynew~x, method='sir')
      #keep significant directions
      #reduce x dimension, only keep 1 dimension, that is how I know it right now
      vec1x <- step2$evectors[,1]
      cor1 <- cbind(cor1,step2$evalues[1])
      #keep signs of each vec2 the same
      if (angle(vec1x, Xvec[[i]])>pi/2)
        vec1x <- -vec1x
      xnew <- x %*% vec1x
      Xvec[[i+1]] <- vec1x
      #step3 use reduced x dimension as y, full y as x, fit MAVE in matlab
      setVariable(matlab,Y=xnew,X=y)
      evaluate(matlab,'[B, cv]= rMAVE(X, Y, 0.5, 1);')
      vec1 <- getVariable(matlab, "B")$B
      # cor2[i] <- sqrt(getVariable(matlab, "cv")$cv)
      #change sign of vec 3 if it's opposite of the original direction
      if (angle(vec1, Yvec[[i]])>pi/2)
        vec1 <- -vec1
      ynew <- y %*% vec1
      Yvec[[i+1]] <- vec1
      i <- i+1
    }
    Xvec <- Reduce(cbind, Xvec)
    Yvec <- Reduce(cbind, Yvec)
    if(sum(vec0!=vec1)+sum(vec0x != vec1x) >0) converge <- FALSE else 
      converge <- TRUE
    
    result1 <- list(method='SIR-MAVE',Xvec=Xvec, Yvec=Yvec, x=x, y=y,
                    x.original =x.original, y.original=y.original, 
                    evalue=cor1, converge=converge)
    #MAVE-SIR
    converge <- NA
    #standardize
    vec1 <- vec1.ini
    vec1x <- vec1x.ini
    Yvec <- list(vec1)
    Xvec <- list(vec1x)
    cor1 <- NA
    #records value of previous iteration, initial value set to 0
    vec0 <- 0*vec1
    vec0x <- 0*vec1x
    #starting point
    ynew <- y %*% vec1
    i <- 1
    while(i <= niters & sum(vec0!=vec1)+sum(vec0x != vec1x) >0)
    {
      vec0 <- vec1
      vec0x <- vec1x
      #use reduced y dimension to find SIR directions
      setVariable(matlab,Y=ynew,X=x)
      evaluate(matlab,'[B, cv]= rMAVE(X, Y, 0.5, 1);')
      vec1x <- getVariable(matlab, "B")$B
      #keep signs of each vec2 the same
      if (angle(vec1x, Xvec[[i]])>pi/2)
        vec1x <- -vec1x
      xnew <- x %*% vec1x
      Xvec[[i+1]] <- vec1x
      #step3 use reduced x dimension as y, full y as x, fit MAVE in matlab
      #use reduced x dimension to find MP directions
      step3 <- dr(xnew~y, method='sir')
      vec1 <-summary(step3)$evectors[,1]
      # vec3[abs(vec3)<rounddown] <- 0
      vec1 <- as.matrix(vec1)
      cor1 <- cbind(cor1, step3$evalues[1])
      # cor2[i] <- sqrt(getVariable(matlab, "cv")$cv)
      #change sign of vec 3 if it's opposite of the original direction
      if (angle(vec1, Yvec[[i]])>pi/2)
        vec1 <- -vec1
      ynew <- y %*% vec1
      Yvec[[i+1]] <- vec1
      i <- i+1
    }
    Xvec <- Reduce(cbind, Xvec)
    Yvec <- Reduce(cbind, Yvec)
    if(sum(vec0!=vec1)+sum(vec0x != vec1x) >0) converge <- FALSE else 
      converge <- TRUE
    result2 <- list(method='MAVE-SIR',Xvec=Xvec, Yvec=Yvec, x=x, y=y,
                    x.original =x.original, y.original=y.original, 
                    evalue=cor1, converge=converge)
  }
  return(list(result1=result1, result2=result2))
}

direction <- function(output,niter=dim(output$Xvec)[2])
{
  vecx <-output$Xvec[,niter]
  vecy <- output$Yvec[,niter]
  evalue <- mean(output$evalue[niter])
  method <- output$method
  return(list(method=method
              ,xdirections=vecx, ydirections=vecy, evalue= evalue,
              x.original=output$x.original,y.original=output$y.original,
              converge=output$converge))
}


#sink('model2.txt')
#model 2 1000 repetitions
time0 <- Sys.time()
set.seed(7202016)
R <- 1000 #number of repititions
n <- 200
refx1 <- c(0,0,0,0,0,1)
refx2 <- c(1,1,0,0,0,0)
dir2.CCA <- array(NA, c(6,2,R))
dir2.SIR <- array(NA, c(6,2,R))
dir2.ECCA <- array(NA, c(6,2,R))
cor2.CCA <- rep(NA,R)
cor2.SIR <- rep(NA,R)
cor2.ECCA <- rep(NA,R)
converge2.SIR <- rep(NA,R)
converge2.ECCA <- rep(NA,R)
ind2.CCA<- rep(NA,R)
ind2.SIR<- rep(NA,R)
ind2.ECCA<- rep(NA,R)
pb <- winProgressBar("progress bar", "0 % done",
                     0, 100)
for(i in 341:R)
{
  #generate data
  x <- matrix(rnorm(5*n, 0,1),n,5)
  y1 <- 0.3*(x[,1]+x[,2])^2 +0.2* rnorm(n,0,1)
  y26 <- matrix(rnorm(5*n, 0,1),n,5)
  y <- cbind(y1, y26)
  x6 <- 0.8*log(abs(y[,3]+y[,4])) + 0.4*rnorm(n, 0, 1)
  x <- cbind(x, x6)
  
  #CCA
  cca <- cancor(x,y)
  #in case the second CCA direction comes out as refx1
  ind.cca <- ifelse(angle(cca$xcoef[,1], refx1) <= angle(cca$xcoef[,1], refx2),
                    1,2)
  ind2.CCA[i] <- ind.cca
  dir2.CCA[,,i] <- cbind(cca$xcoef[,ind.cca],cca$ycoef[,ind.cca])
  cor2.CCA[i] <- cca$cor[ind.cca]
  
  #multi SIR
  sir <- multiSIR(x,y, niters = 100)
  converge2.SIR[i] <- sir$converge
  if(sir$converge==TRUE)
  {
    nn <- dim(sir$Xvec)[2]
    #make sure the first element is positive before comparing angles
    xvec <- sir$Xvec[,nn]* sign(sir$Xvec[1,nn])
    yvec <- sir$Yvec[,nn]* sign(sir$Yvec[1,nn])
    if(angle(xvec,refx1)<= angle(xvec, refx2))
    {
      dir2.SIR[,,i] <- cbind(xvec,yvec)
      ind2.SIR[i] <- 1
    } else ind2.SIR[i] <- 2
    cor2.SIR[i] <- sir$evalue[nn]
  }
  
  #ECCA
  ecca2 <- ECCA(x,y)
  dir1.forward <- direction(ecca2$result1)
  dir1.inverse <- direction(ecca2$result2)
  if(dir1.forward$converge==TRUE)
  {
    if(dir1.forward$evalue>=dir1.inverse$evalue |
       dir1.inverse$converge==FALSE)
      ecca <- dir1.forward else
        ecca <- dir1.inverse
  } else if(dir1.inverse$converge==TRUE)
    ecca <- dir1.inverse else
    {
      if(dir1.forward$evalue>=dir1.inverse$evalue)
        ecca <- dir1.forward else
          ecca <- dir1.inverse
    }
  converge2.ECCA[i] <- ecca$converge
  if(ecca$converge==TRUE)
  {
    xvec <- ecca$xdirections * sign(ecca$xdirections[1])
    yvec <- ecca$ydirections * sign(ecca$ydirections[1])
    if(angle(xvec,refx1)<= angle(xvec, refx2))
    {
      dir2.ECCA[,,i] <- cbind(xvec,yvec)
      ind2.ECCA[i] <- 1
    } else ind2.ECCA[i] <- 2
    cor2.ECCA[i] <- ecca$evalue
  }
  setWinProgressBar(pb,100*i/R, label = sprintf("%f%% done",100*i/R))
  if(i %% 20 ==0)
    save.image("run0720.RData")
}
close(pb)
Sys.time() - time0


#sink('model1.txt')
#model 1 repeat 1000 times
time0 <- Sys.time()
time0
set.seed(7212016)
R <- 1000 #number of repititions
n <- 200
refx1 <- c(1,-1,0,0,0,0)
refx2 <- c(1,1,0,0,0,0)
dir1.CCAx <- matrix(NA,6,R)
dir1.SIRx <- matrix(NA,6,R)
dir1.ECCAx <- matrix(NA,6,R)
dir1.CCAy <- matrix(NA,9,R)
dir1.SIRy <- matrix(NA,9,R)
dir1.ECCAy <- matrix(NA,9,R)
cor1.CCA <- rep(NA,R)
cor1.SIR <- rep(NA,R)
cor1.ECCA <- rep(NA,R)
converge1.SIR <- rep(NA,R)
converge1.ECCA <- rep(NA,R)
ind1.CCA<- rep(NA,R)
ind1.SIR<- rep(NA,R)
ind1.ECCA<- rep(NA,R)
pb <- winProgressBar("progress bar", "0 % done",
                     0, 100)
for(i in 1:R)
{
  #generate data
  x <- matrix(rnorm(6*n,0,1),n,6)
  eps <- mvrnorm(n, c(0,0), matrix(c(2,-1,-1,1),2,2))
  mu1 <- sin(x[,1]+x[,2])
  y12 <- cbind(mu1, 2*mu1) + 0.2*eps
  y3 <- rnorm(n,0,1)
  y34 <- log(x[,1] - x[,2] +5) + 0.2* rnorm(n,0,1)
  y4 <- y34-y3
  y5_9 <- matrix(rnorm(5*n, 0,1),n,5)
  y <- cbind(y12, y3, y4, y5_9)
  
  #CCA
  cca <- cancor(x,y)
  #in case the second CCA direction comes out as refx1
  ind.cca <- ifelse(angle(cca$xcoef[,1], refx1) <= angle(cca$xcoef[,1], refx2),
                    1,2)
  ind1.CCA[i] <- ind.cca
  dir1.CCAx[,i] <- cca$xcoef[,ind.cca]
  dir1.CCAy[,i] <- cca$ycoef[,ind.cca]
  cor1.CCA[i] <- cca$cor[ind.cca]
  
  #multi SIR
  sir <- multiSIR(x,y, niters = 100)
  converge1.SIR[i] <- sir$converge
  if(sir$converge==TRUE)
  {
    nn <- dim(sir$Xvec)[2]
    #make sure the first element is positive before comparing angles
    xvec <- sir$Xvec[,nn]* sign(sir$Xvec[1,nn])
    yvec <- sir$Yvec[,nn]* sign(sir$Yvec[1,nn])
    if(angle(xvec,refx1)<= angle(xvec, refx2))
    {
      dir1.SIRx[,i] <- xvec
      dir1.SIRy[,i] <- yvec
      ind1.SIR[i] <- 1
    } else ind1.SIR[i] <- 2
    cor1.SIR[i] <- sir$evalue[nn]
  }
  
  #ECCA
  ecca2 <- ECCA(x,y)
  dir1.forward <- direction(ecca2$result1)
  dir1.inverse <- direction(ecca2$result2)
  if(dir1.forward$converge==TRUE)
  {
    if(dir1.forward$evalue>=dir1.inverse$evalue |
       dir1.inverse$converge==FALSE)
      ecca <- dir1.forward else
        ecca <- dir1.inverse
  } else if(dir1.inverse$converge==TRUE)
    ecca <- dir1.inverse else
    {
      if(dir1.forward$evalue>=dir1.inverse$evalue)
        ecca <- dir1.forward else
          ecca <- dir1.inverse
    }
  converge1.ECCA[i] <- ecca$converge
  if(ecca$converge==TRUE)
  {
    xvec <- ecca$xdirections * sign(ecca$xdirections[1])
    yvec <- ecca$ydirections * sign(ecca$ydirections[1])
    if(angle(xvec,refx1)<= angle(xvec, refx2))
    {
      dir1.ECCAx[,i] <- xvec
      dir1.ECCAy[,i] <- yvec
      ind1.ECCA[i] <- 1
    } else ind1.ECCA[i] <- 2
    cor1.ECCA[i] <- ecca$evalue
  }
  setWinProgressBar(pb,100*i/R, label = sprintf("%f%% done",100*i/R))
  if(i %% 20 ==0)
  save.image("run0720.RData")
}
close(pb)
Sys.time() - time0


#just the first model 2 direction
#sink('model21.txt')
#model 2 1000 repetitions
time0 <- Sys.time()
set.seed(7212016)
R <- 1000 #number of repititions
n <- 200
refx1 <- c(0,0,0,0,0,1)
refx2 <- c(1,1,0,0,0,0)
dir21.CCA <- array(NA, c(6,2,R))
dir21.SIR <- array(NA, c(6,2,R))
dir21.ECCA <- array(NA, c(6,2,R))
cor21.CCA <- rep(NA,R)
cor21.SIR <- rep(NA,R)
cor21.ECCA <- rep(NA,R)
converge21.SIR <- rep(NA,R)
converge21.ECCA <- rep(NA,R)
ind21.CCA<- rep(NA,R)
ind21.SIR<- rep(NA,R)
ind21.ECCA<- rep(NA,R)
pb <- winProgressBar("progress bar", "0 % done",
                     0, 100)
for(i in 1:R)
{
  #generate data
  x <- matrix(rnorm(5*n, 0,1),n,5)
  y <- matrix(rnorm(6*n, 0,1),n,6)
  x6 <- 0.8*log(abs(y[,3]+y[,4])) + 0.4*rnorm(n, 0, 1)
  x <- cbind(x, x6)
  
  #CCA
  cca <- cancor(x,y)
  #in case the second CCA direction comes out as refx1
  ind.cca <- ifelse(angle(cca$xcoef[,1], refx1) <= angle(cca$xcoef[,1], refx2),
                    1,2)
  ind21.CCA[i] <- ind.cca
  dir21.CCA[,,i] <- cbind(cca$xcoef[,ind.cca],cca$ycoef[,ind.cca])
  cor21.CCA[i] <- cca$cor[ind.cca]
  
  #multi SIR
  sir <- multiSIR(x,y, niters = 100)
  converge21.SIR[i] <- sir$converge
  if(sir$converge==TRUE)
  {
    nn <- dim(sir$Xvec)[2]
    #make sure the first element is positive before comparing angles
    xvec <- sir$Xvec[,nn]* sign(sir$Xvec[1,nn])
    yvec <- sir$Yvec[,nn]* sign(sir$Yvec[1,nn])
    if(angle(xvec,refx1)<= angle(xvec, refx2))
    {
      dir21.SIR[,,i] <- cbind(xvec,yvec)
      ind21.SIR[i] <- 1
    } else ind21.SIR[i] <- 2
    cor21.SIR[i] <- sir$evalue[nn]
  }
  
  #ECCA
  ecca2 <- ECCA(x,y)
  dir1.forward <- direction(ecca2$result1)
  dir1.inverse <- direction(ecca2$result2)
  if(dir1.forward$converge==TRUE)
  {
    if(dir1.forward$evalue>=dir1.inverse$evalue |
       dir1.inverse$converge==FALSE)
      ecca <- dir1.forward else
        ecca <- dir1.inverse
  } else if(dir1.inverse$converge==TRUE)
    ecca <- dir1.inverse else
    {
      if(dir1.forward$evalue>=dir1.inverse$evalue)
        ecca <- dir1.forward else
          ecca <- dir1.inverse
    }
  converge21.ECCA[i] <- ecca$converge
  if(ecca$converge==TRUE)
  {
    xvec <- ecca$xdirections * sign(ecca$xdirections[1])
    yvec <- ecca$ydirections * sign(ecca$ydirections[1])
    if(angle(xvec,refx1)<= angle(xvec, refx2))
    {
      dir21.ECCA[,,i] <- cbind(xvec,yvec)
      ind21.ECCA[i] <- 1
    } else ind21.ECCA[i] <- 2
    cor21.ECCA[i] <- ecca$evalue
  }
  setWinProgressBar(pb,100*i/R, label = sprintf("%f%% done",100*i/R))
  if(i %% 20 ==0)
  save.image("run0720.RData")
}
close(pb)
Sys.time() - time0
