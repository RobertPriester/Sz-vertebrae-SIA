
# Acoompanying functions to calculate SEV, SEVc, SEVb and ellipsoid overlaps for Skinner et al., 2019

#############################################################################################

#############################################################################################

### Estimating SEV and SEVc ###

# EstimateSEVc() estimates the ellipsoid volume of a dataset under the assumption of a multivariate normal distribution
# via the sample covariance matrix estimated via maximum likelihood.
# The argument sample.size.correction defaults to TRUE giving SEVc, to obtain SEV, turn this argument to FALSE
# Input data should be given in a standard format (dataframe, matrix or 2d array) with each row corresponding to an 
# observation and each column a variable
# There is the option give a desired p.interval (proportion of data to encompass), otherwise the default returns the
# standard ellipsoid which in the 3d case covers ~20% of the data

# Note that this is a general function and should work on any n-dimensional dataset
EstimateSEVc <- function(data, sample.size.correction = TRUE, p.interval = NULL) {
  if(is.null(p.interval)) {message("No p.interval provided - defaulting to standard ellipsoid")}
  
  m <- ncol(data)
  n <- nrow(data)
  cov_mat <- cov(data)
  eig <- eigen(cov_mat)
  sqrt(eig$values)
  
  if(is.null(p.interval)){p.interval <- stats::pchisq(1, m)}
  
  if(sample.size.correction==FALSE) {scaling <- sqrt(stats::qchisq(p.interval, m))}
  if(sample.size.correction==TRUE) {scaling <- sqrt(qf(p.interval, m, (n-m))*m*(n-1)/(n-m))}
  
  x <- m/2
  (pi^(x))/gamma(x+1)
  
  ellipsoid.volume <- (scaling^m)*prod(sqrt(eig$values))*(pi^(x))/gamma(x+1)
  return(ellipsoid.volume)
}


#############################################################################################

#############################################################################################

# Bayesian determination of covariance and means of n-dimensional dataset
# pre-requisit for SEVb and overlaps functions

# the fitMVNdirect() function estimates the means and covariance matrix of a 3 dimensional multivariate normal distribution
# to data, without any rescaling or translation of data prior to Bayesian estimation. It is directly modified from
# the fitEllipse() function provided in the SIBER package.
# Sensistivity analysis indicates that this approach exhibits less bias compared to fitMVN(), which appears to over-estimate
# volume at smaller sample sizes
# The default burnin, iterations and thinning are set to those used in the sensitivity analysis
# The priors - R, K and tau.mu - which are the scale matrix and degrees of freedom of the Wishart distribution, and mean precision
# (centred about zero) are defaulted to vague, uniformative priors.

fitMVNdirect <- function (data,
                          nChain = 3,
                          nAdapt = 1000,
                          nBurnIn = 10000,
                          nThin = 25,
                          nIter = 100000,  ## altered!!
                          save.model.txt = FALSE,
                          R = diag(1,3),
                          k = ncol(data),
                          tau.mu = 10^-3) 
{
  modelstring <- "\n
  \n    model {
  \n    # ----------------------------------
  \n    # define the priors
  \n    # ----------------------------------
  \n    
  \n    # this loop defines the priors for the means
  \n    for (i in 1:n.iso) {
  \n      mu[i] ~ dnorm (0, tau.mu)
  \n    }
  \n    
  \n    # prior for the precision matrix
  \n    tau[1:n.iso,1:n.iso] ~ dwish(R[1:n.iso,1:n.iso],k)
  \n    
  \n    # convert to covariance matrix
  \n    Sigma2[1:n.iso, 1:n.iso] <- inverse(tau[1:n.iso, 1:n.iso]) 
  \n    
  \n    # calculate correlation coefficient
  \n    # rho <- Sigma2[1,2]/sqrt(Sigma2[1,1]*Sigma 2[2,2])
  \n    
  \n    #----------------------------------------------------
  \n    # specify the likelihood of the observed data
  \n    #----------------------------------------------------
  \n    
  \n    for(i in 1:n.obs) {                             
  \n      Y[i,1:n.iso] ~ dmnorm(mu[1:n.iso],tau[1:n.iso,1:n.iso])
  \n    }
  \n    
  \n    
  \n    
  \n  }"
  Y = data
  n.obs <- nrow(Y)
  n.iso <- ncol(Y)
  jags.data <- list(Y = Y, n.obs = n.obs, n.iso = n.iso, R = R, 
                    k = k, tau.mu = tau.mu)
  inits <- list(list(mu = stats::rnorm(n.iso, 0, 1)), list(mu = stats::rnorm(n.iso, 
                                                                             0, 1)))
  parameters <- c("mu", "Sigma2")
  model <- rjags::jags.model(textConnection(modelstring), data = jags.data, 
                             n.chains = nChain)
  output <- rjags::coda.samples(model = model, variable.names = c("mu", 
                                                                  "Sigma2"), n.iter = nIter, thin = nThin)
  return(output)
}

#############################################################################################

#############################################################################################

# Estimating SEVb

# Function fit3dEllipsoid() notes
# This function takes the output from fitMVNdirect() for the 3d case only!
# The default setting returns the volumes of the standard 3d ellipsoid.
# One can specify a proportion of data coverage with the p.interval argument (note that the standard
# ellipsoid covers approximately only 20% of the data!).
# Because it uses Bayesian draws, there is no need for sample size correction in this instance.
# If magnitude and vectors of major and minor axes are desired, set vol.only to FALSE

fit3dEllipsoid <- function(fitMVNdirect.output, vol.only = TRUE, p.interval = NULL) {
  if(!is.list(fitMVNdirect.output)){
    stop("Incorrect data input")
  }
  if(is.null(p.interval)) {message("No p.interval provided - defaulting to standard ellipsoid, approximately 20% data coverage")}
  MCMCdata <- as.matrix(fitMVNdirect.output) # Combine the MCMC chains
  if(ncol(MCMCdata) != 12){
    stop("Incorrect data input - none 3d dataset used")
  }
  if(vol.only){Ellipsoid_volume <- matrix(NA, ncol = 1, nrow = nrow(MCMCdata))
  colnames(Ellipsoid_volume) <- c("Volume")
  }else{Ellipsoid_volume <- matrix(NA, ncol = 13, nrow = nrow(MCMCdata))
  colnames(Ellipsoid_volume) <- c("Volume", "a", "b", "c", "VectorA1", "VectorA2", "VectorA3",
                                  "VectorB1", "VectorB2", "VectorB3", "VectorC1", "VectorC2",
                                  "VectorC3")}
  
  # Determine the scaling to be applied to a, b & c depending on whether we 1) correct for sample size and 2)
  # the proportion of data to cover
  
  if(is.null(p.interval)) {scaling <- 1}
  if(!is.null(p.interval)) {scaling <- sqrt(stats::qchisq(p.interval, 3))}
  
  for(i in 1:nrow(MCMCdata)) {
    # Firstly, determine the covariance matrix from sd's and rho's
    cov_matrix <- matrix(MCMCdata[i,1:9], ncol = 3, byrow = TRUE)
    
    # Decompose the covariance matrix into eigenvalues and eigenvectors
    # This decomposition provides a rotational matrix (~eigenvectors) and linear scalar (~eigenvalues)
    # that describes the transformation of "white noise" (i.e. variances of 1 along x,y,z 
    # dimensions with no covariance structure) into the covariance matrix.
    # Because the input is a covariance matrix, the eigenvalues are the maximal variances in
    # n-orthoganal directions, with the eigenvectors as the corresponding unit vectors of these
    # maximal variances in x,y,z (i.e. the "new" axes).
    
    eig <- eigen(cov_matrix)
    a <- sqrt(eig$values[1])*scaling # we take the sqrt because the input is variances and we want sd.
    b <- sqrt(eig$values[2])*scaling
    c <- sqrt(eig$values[3])*scaling
    Vol <- a*b*c*pi*4/3
    
    if(vol.only){
      Ellipsoid_volume[i,1] <- Vol
    }else{
      Ellipsoid_volume[i,1] <- Vol
      Ellipsoid_volume[i,2] <- a
      Ellipsoid_volume[i,3] <- b
      Ellipsoid_volume[i,4] <- c
      Ellipsoid_volume[i,5] <- eig$vectors[1,1]
      Ellipsoid_volume[i,6] <- eig$vectors[2,1]
      Ellipsoid_volume[i,7] <- eig$vectors[3,1]
      Ellipsoid_volume[i,8] <- eig$vectors[1,2]
      Ellipsoid_volume[i,9] <- eig$vectors[2,2]
      Ellipsoid_volume[i,10] <- eig$vectors[3,2]
      Ellipsoid_volume[i,11] <- eig$vectors[1,3]
      Ellipsoid_volume[i,12] <- eig$vectors[2,3]
      Ellipsoid_volume[i,13] <- eig$vectors[3,3]
    }
  }
  return(Ellipsoid_volume)
}

#############################################################################################

#############################################################################################

# Estimating overlap of two ellipsoids

# Estimating the volume overlap means first defining the intersection of two ellipsoids
# Much work has been done on intersecting ellipsoids mostly due to its applications for 
# rendering computer models and simulation, and typically fall into algorithms for testing
# whether two ellipsoids intersect or not (i.e. do two objects collide?) or are based on 
# optimising 'packing' problems (i.e. what is the best way to fill a defined space with ellipsoids,
# with or without overlap)

# I have been unable to find an arithmetic solution to estimate the volume of intersection - 
# most likely because the two issues regarded above do not necessarily require an absolute measure
# of intersection volume (but also because mathematics is not my strong point!). For issues around
# packing, most studies seem to approximate volume by optimally fitting another ellipsoid within
# the intersection to measuring the volume of that... I do not think this is necessarily appropiate
# here, as it will always lead to an underestimation of volume and therefore perceived overlap
# (although I am certain the underestimation won't be large, regardless of intersection shape)

# This leaves a numerical solution (approximation) of the volume which is computationally easier
# to compute. SIBER estimates ellipse overlap by first approximating the two ellipses as polygons
# with a high number of vertices located on the ellipse perimeters. The intersection of these two
# polygons is then defined (as a polygon convex hull of the overlapping vertices), from which area
# can be calculated.

# The 3d equivalent could be applied, which would involve approximating the two ellipsoids as 3d
# meshes, define the intersection of these meshes as a convex hull, from which a volume can be
# estimated. 

library('rgl') # This package contains functions to produce 3d meshes from defined 3 dimensional shapes
library('geometry') # This package contains geometric functions that interface well with rgl mesh objects
# and can be used to define intersections

# the ellipse3d function in 'rgl' takes a covariance matrix, a centre (3 means), and a p.interval.
# So we can take the output from fitMVNdirect() as before and construct the covariance matrix


ellipsoid.Overlap <- function(ellipsoid1, ellipsoid2,
                              p.interval = NULL, subdivide = 4) {
  message("Please note that this function may take a while to complete depending on the number
          of posterior draws - run time can be reduced by lowering the subdivide value but this results
          in poorer ellipsoid approximation.")
  require('rgl')
  require('geometry')
  MCMCdata1 <- as.matrix(ellipsoid1) # Convert posteriors into matrices
  MCMCdata2 <- as.matrix(ellipsoid2)
  n1 <- nrow(MCMCdata1) # Lengths of draws, check if they are equal and limit larger if not
  n2 <- nrow(MCMCdata2)
  if(n1 != n2) {
    message(paste("Different number of posterior draws, using first", ifelse(n1>n2, n2, n1),
                  "only to estimate overlap volumes"))
    ifelse(n1>n2, MCMCdata1 <- MCMCdata1[1:n2,], MCMCdata2 <- MCMCdata2[1:n1,])
  }
  if(is.null(p.interval)) {message("No p.interval given - defaulting to the standard ellipsoid, which encompasses approximately 20% of the data")}
  if(is.null(p.interval)) {p.interval <- pchisq(1,3)}
  # Create dataframe to fill and return
  Overlaps <- as.data.frame(matrix(NA, nrow = nrow(MCMCdata1), ncol = 5))
  colnames(Overlaps) <- c("Overlap.Volume", "Ellipsoid1.Volume", "Ellipsoid2.Volume",
                          "Prop.Overlap.Ellipsoid1", "Prop.Overlap.Ellipsoid2")
  cf1 <- 1
  cf2 <- 1
  
  # setup progress bar
  pb <- txtProgressBar(min = 0, max = nrow(MCMCdata1), style = 3)
  
  for(i in 1:nrow(MCMCdata1)){
    # Calculate means and covariance matrix for ellipsoid 1
    means1 <- MCMCdata1[i, 10:12]
    cov_matrix1 <- matrix(MCMCdata1[i,1:9], ncol = 3, byrow = FALSE)
    # Calculate means and covariance matrix for ellipsoid 2
    means2 <- MCMCdata2[i, 10:12]
    cov_matrix2 <- matrix(MCMCdata2[i,1:9], ncol = 3, byrow = FALSE)
    
    # Create 3d mesh approximations of the two ellipsoids
    ellipsoid1 <- ellipse3d(cov_matrix1, centre = means1, level = p.interval, subdivide = subdivide)
    ellipsoid2 <- ellipse3d(cov_matrix2, centre = means2, level = p.interval, subdivide = subdivide)
    
    # Calculate the convex hull of the intersection of the two ellipsoids
    intersection <- intersectn(t(ellipsoid1[["vb"]])[,1:3], t(ellipsoid2[["vb"]])[,1:3], tol=0,
                               return.chs = TRUE)
    # Store information
    Overlaps[i,1] <- round(intersection[['ch']]$vol, 3)
    Overlaps[i,2] <- round(intersection[['ch1']]$vol, 3)
    Overlaps[i,3] <- round(intersection[['ch2']]$vol, 3)
    Overlaps[i,4] <- round(intersection[['ch']]$vol/intersection[['ch1']]$vol, 4)
    Overlaps[i,5] <- round(intersection[['ch']]$vol/intersection[['ch2']]$vol, 4)
    
    # Update progress bar
    setTxtProgressBar(pb, i)
  }
  
  return(Overlaps)
}

#############################################################################################

#############################################################################################

### Plotting SEVc ###

# The ellipse3d() function provides a means to produce meshes of ellipsoids. Although this function provides a p.interval
# (argument called level), it does not have a sample size correction.
# Here, the function ellipse3d.samplesize() provides a correction for sample size before passing arguments to the ellipse3d()
# function by defining the appropriate t-statistic on the boundary

ellipse3d.samplesize <- function(cov_matrix, means, sample.size = NULL, p.interval=NULL, subdivide = 5) {
  require('rgl')
  if(is.null(sample.size)) {stop("No sample size provided")}
  if(is.null(p.interval)) {message("No p.interval given - defaulting to the standard ellipsoid, which excompasses approximately 20% of the data")}
  if(is.null(p.interval)) {p.interval <- pchisq(1,3)}
  cf1 <- 3*(sample.size - 1)/(sample.size - 3)
  return(ellipse3d(cov_matrix, centre = means, t = sqrt(qf(p.interval, 3, (sample.size-3))*cf1), subdivide = subdivide))
}

#####################################################################################################################
#####################################################################################################################
#####################################################################################################################






