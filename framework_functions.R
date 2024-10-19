library(adiv)
library(Select)
library(fundiversity)
library(RColorBrewer)

# function to generate simulated communities: ##############
comSimulation <- function(trait, ava, it, rich, cwm, cwv, rao, cost, dens, ref,
                          min_p, phi = 1){
  # generates simulated communities and calculate its parameters.
  # ARGUMENTS:
  # trait: data frame or matrix with species traits. Traits as columns and species as rows.
  # ava: vector indicating availability of species (1 or 0)
  # it: number of iterations (communities)
  # rich: range of richness values
  # cwm: character vector with trait names to calculate Community Weighted Mean (CWM). One CWM is calculated for each trait.
  # cwv: vector vector with trait names to calculate Community Weighted Variance. One CWV is calculated for each trait.
  # rao: vector with traits to calculate Rao Quadratic Entropy, or distance matrix (class dist)
  # cost: vector of species cost per individual
  # dens: vector of planting density
  # ref: matrix with proportion of species in the reference sites. NAs not accepted.
  
  # OUTPUTS ###############################################
  simComm <- list(composition = c(), parameters = c())
  refComm <- list(composition = c(), parameters = c())
  
  # CREATE RANDOM COMMUNITIES #############################
  species <- row.names(trait)
  resample <- function(x, ...) x[sample.int(length(x), ...)]
  
  # Find distant species: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  if(!missing(rao)){
    if(inherits(rao, 'character')){
      t2d <- as.matrix(scale(trait[, rao, drop = FALSE]) )
    } else if(inherits(rao, 'dist')){
      t2d <- rao
    }
    if(!missing(cwm)){
      t2c <- as.matrix(trait[,cwm, drop = FALSE])
      constraints <- apply(t2c, 2, function(x){
        cons <- seq(min(x), max(x), length.out = 8)
        return(cons[c(-1,-8)])
      })
      selSpp <- vector(mode="list", length=length(cwm))
      names(selSpp) <- cwm
      for(i in cwm){
        cons_i <- constraints[,i]
        t2c_i <- t2c[,i, drop=FALSE]
        selSpp_i <- c()
        for(j in cons_i){
          names(j) <- colnames(t2c_i)
          invisible(capture.output(selSpp_j <- selectSpecies(t2c_i, j, t2d, phi = phi) ) ) 
          selSpp_i <- cbind(selSpp_i, selSpp_j$prob)
        }
        selSpp[[i]] <- selSpp_i 
      }
      propMatrixSelSpp <- do.call(cbind, selSpp)
      sppMax <- c()
      propMin <- 0.1
      while(length(sppMax) < rich[1]){
        sppMax <- which(propMatrixSelSpp > propMin, arr.ind = TRUE)
        sppMax <- unique(species[sppMax[,1]])
        propMin <- 0.5*propMin
      }
      propMatrixSelSpp <- round(t(propMatrixSelSpp),3)
    } else {
      invisible( capture.output( selSpp <- selectSpecies(t2d = t2d, phi = phi)$prob ))
      propMatrixSelSpp <- selSpp
      sppMax <- c()
      propMin <- 0.1
      while(length(sppMax) < rich[1]){
        sppMax <- species[propMatrixSelSpp > propMin]
        propMin <- 0.5*propMin
      }
    }
    
    # Find distant species that are available: <<<<<<<<<<<<<<
    avaLog <- as.logical(ava)
    speciesAva <- species[avaLog]
    if(inherits(rao, 'character')){
      t2d <- as.matrix(scale(trait[avaLog, rao, drop = FALSE]) )
    } else if(inherits(rao, 'dist')){
      t2d <- as.dist(as.matrix(rao)[avaLog, avaLog])
    }
    if(!missing(cwm)){
      t2c <- as.matrix(trait[avaLog, cwm, drop = FALSE])
      constraints <- apply(t2c, 2, function(x){
        cons <- seq(min(x), max(x), length.out = 8)
        return(cons[c(-1,-8)])
      })
      selSppAva <- vector(mode="list", length=length(cwm))
      names(selSppAva) <- cwm
      for(i in cwm){
        cons_i <- constraints[,i]
        t2c_i <- t2c[,i, drop=FALSE]
        selSppAva_i <- c()
        for(j in cons_i){
          names(j) <- colnames(t2c_i)
          invisible(capture.output(selSppAva_j <- selectSpecies(t2c_i, j, t2d, phi = phi) ))
          selSppAva_i <- cbind(selSppAva_i, selSppAva_j$prob)
        }
        selSppAva[[i]] <- selSppAva_i 
      }
      propMatrixSelSppAva <- do.call(cbind, selSppAva)
      sppMaxAva <- c()
      propMin <- 0.1
      while(length(sppMaxAva) < rich[1]){ 
        sppMaxAva <- which(propMatrixSelSppAva > propMin, arr.ind = TRUE)
        sppMaxAva <- unique(speciesAva[sppMaxAva[,1]])
        propMin <- 0.5*propMin
      }
      propMatrixSelSppAva <- round(t(propMatrixSelSppAva),3)
    } else {
      invisible(capture.output(selSppAva <- selectSpecies(t2d = t2d, phi = phi)$prob ))
      propMatrixSelSppAva <- selSpp
      sppMaxAva <- c() 
      propMin <- 0.1 
      while(length(sppMaxAva) < rich[1]){
        sppMaxAva <- species[propMatrixSelSppAva > propMin]
        propMin <- 0.5*propMin
      }
    }
    
    # number of iterations for simulations: <<<<<<<<<<<<<<<<<
    itMax <- round(0.25*it)
    itMaxAva <- round(0.25*it)
    itAva <- round(0.25*it)
    itAll <- it - itMax - itMaxAva - itAva
    
    # Run simulation with species that maximize rao: <<<<<<<<<
    if(length(sppMax) < rich[2]){
      nsp <- length(sppMax)
      vLen <- nsp
    } else {
      nsp <- rich[2]
      vLen <- length(sppMax)
    }
    propMatrixSelSpp2 <- matrix(rep(0,length(species)*itMax),
                                ncol=length(species), nrow=itMax)
    sppMaxPos <- species %in% sppMax
    for(i in 1:itMax){
      nsp_i <-  resample(rich[1]:nsp, 1)
      ocor = sample( c(rep(1, nsp_i), rep(0, vLen - nsp_i)) )
      abund = rlnorm(vLen)
      abund <- abund * ocor
      prop <- abund/sum(abund)
      propMatrixSelSpp2[i,sppMaxPos] <- prop
    }
    
    # Run simulation with species that maximize rao and are available: <<<<
    if(length(sppMaxAva) < rich[2]){
      nsp <- length(sppMaxAva) 
      vLen <- nsp
    } else {
      nsp <- rich[2]
      vLen <- length(sppMaxAva)
    }
    propMatrixSelSppAva2 <- matrix(rep(0,length(species)*itMaxAva),
                                   ncol=length(species), nrow=itMaxAva)
    sppMaxAvaPos <- species %in% sppMaxAva
    for(i in 1:itMaxAva){
      nsp_i <-  resample(rich[1]:nsp, 1)
      ocor = sample( c(rep(1, nsp_i), rep(0, vLen - nsp_i)) )
      abund = rlnorm(vLen)
      abund <- abund * ocor
      prop <- abund/sum(abund)
      propMatrixSelSppAva2[i,sppMaxAvaPos] <- prop
    }
  }
  else {
    avaLog <- as.logical(ava)
    itAva <- round(0.5*it)
    itAll <- it - itAva
  }
  
  # Run simulation with available species: <<<<<<<<<<<<<<<<
  if(sum(ava) < rich[2]){
    nsp <- sum(ava) 
    vLen <- nsp
  } else {
    nsp <- rich[2]
    vLen <- sum(ava)
  }
  propMatrixAva <- matrix(rep(0,length(species)*itAva),
                          ncol=length(species), nrow=itAva)
  if(rich[1] > nsp){
    stop("Minimum richness is higher than number of available species")
  }
  for(i in 1:itAva){
    nsp_i <-  resample(c(rich[1]:nsp), 1)
    ocor = sample( c(rep(1, nsp_i), rep(0, vLen - nsp_i)) )
    abund = rlnorm(vLen)
    abund <- abund * ocor
    prop <- abund/sum(abund)
    propMatrixAva[i,avaLog] <- prop
  }
  
  # Run simulation with all species: <<<<<<<<<<<<<<<<<<<<<<
  nsp <- length(species)
  propMatrixPool <- matrix(rep(0,length(species)*itAll),
                           ncol=length(species), nrow=itAll)
  for(i in 1:itAll){
    nsp_i <-  resample(c(rich[1]:rich[2]), 1)
    ocor = sample( c(rep(1, nsp_i), rep(0, nsp - nsp_i)) )
    abund = abund = rlnorm(nsp)
    abund <- abund * ocor
    prop <- abund/sum(abund)
    propMatrixPool[i,] <- prop
  }
  
  if(!missing(rao)){
    propMatrix <- rbind(propMatrixSelSpp2,propMatrixSelSppAva2,
                        propMatrixAva, propMatrixPool)
  } else {
    propMatrix <- rbind(propMatrixAva, propMatrixPool)
  }
  
  rownames(propMatrix) <- sprintf("sim%d",seq(1:nrow(propMatrix)))
  colnames(propMatrix) <- species
  
  simComm$composition <- propMatrix
  
  # Exclude rare species:
  if(!missing(min_p)){
    pos <- propMatrix < min_p
    propMatrix[pos] <- 0
    propMatrix <- propMatrix/rowSums(propMatrix)
    pos <- is.na(propMatrix[,1])
    propMatrix <- propMatrix[!pos,]
  }
  
  # CALCULATE PARAMETERS ##################################
  if(!require(adiv)){
    stop("Package adiv not found")
  }
  if(!require(fundiversity)){
    stop("Package fundiversity not found")
  }
  
  out <- NULL
  
  if(!missing(ref)){
    ref <- ref/rowSums(ref)
    x <- rbind(propMatrix, ref)
  } else {x <- propMatrix}
  
  S <- apply(x, 1, FUN = function(a) sum(a > 0))
  out <- cbind(out, richness = S)
  
  UNA <- apply(x, 1, FUN = function(a) sum(a[!as.logical(ava)] > 0) )
  out <- cbind(out, unavailable = UNA)
  
  if(!missing(cwm)){
    if(inherits(cwm, 'character')){
      traitSub <- trait[,cwm, drop=FALSE]
      CWM <- apply(x, 1, FUN=function(p){
        colSums(traitSub*p, na.rm = T)
      })
      if(class(CWM)[1] == 'matrix'){
        CWM <- t(CWM)
      } else {
        CWM <- as.matrix(CWM)
      }
      #colnames(CWM) <- paste0('cwm_', cwm)
      out <- cbind(out, CWM)
    } else {message('** CWM skipped. **')}
  }
  
  if(!missing(cwv)){
    if(inherits(cwv, 'character')){
      traitSub <- trait[,cwv, drop=FALSE]
      CWM <- apply(x, 1, FUN=function(p){
        colSums(traitSub*p, na.rm = T)
      })
      if(class(CWM)[1] == 'matrix'){
        CWM <- t(CWM)
      } else {
        CWM <- as.matrix(CWM)
      }
      
      CWV <- NULL
      for(i in 1:nrow(CWM)){ #para cada linha de CWMs
        CWM_i <- CWM[i,] #linha i
        traitMod <- NULL #traits modificados pelo CWM i
        for(j in 1:ncol(traitSub)){
          traitMod_j <- traitSub[,j] - CWM_i[j] #subtrai pelo CWM
          traitMod_j <- traitMod_j^2 #eleva ao quadrado
          traitMod <- cbind(traitMod, traitMod_j) #adiciona na tabela
        }
        #CWV: traits modificados pelo CWM i vezes proporcao i
        CWV_i <- colSums(traitMod*x[i,])
        CWV <- rbind(CWV, CWV_i)
      }
      #colnames(CWV) <- paste0('cwv_', cwv)
      colnames(CWV) <- cwv
      CWV <- apply(CWV, 2, FUN = function(x){x/max(x)}) #padroniza pelo maximo
      out <- cbind(out, CWV)
    } else {message('** CWV skipped. **')}
  }
  
  if(!missing(rao)){
    if(inherits(rao, 'character')){
      traitSub <- scale(trait[,rao, drop=FALSE] )
      RAO <- fd_raoq(traitSub, x)$Q
    } else if(inherits(rao, 'dist')){
      RAO <- fd_raoq(sp_com = x, dist_matrix = rao)$Q
    } else{
      message('** RAO skipped')
    }
    RAO <- RAO/max(RAO)
    out <- cbind(out, rao = RAO)
    colnames(out)[colnames(out) == 'rao'] <- rao
  }
  
  if(!missing(cost)){
    COST <- apply(x, 1, FUN=function(p){
      COST_i <- sum(p*cost*dens, na.rm = TRUE)
      return(COST_i)
    })
    out <- cbind(out, cost = COST)
  }
  
  pos <- duplicated(colnames(out)) | duplicated(colnames(out), fromLast=TRUE)
  dup_col <- colnames(out)[pos]
  correct_col <- paste0(dup_col, 1:length(dup_col)) #add index to duplicated
  colnames(out)[pos] <- correct_col
  
  pos <- 1:nrow(propMatrix)
  out1 <- as.data.frame(out[pos,])
  simComm$parameters <- out1
  
  if(!missing(ref)){
    pos <- ( nrow(propMatrix)+1 ):( nrow(propMatrix)+nrow(ref) )
    out2 <- as.data.frame(out[pos,,drop = FALSE])
    refComm$composition <- ref
    refComm$parameters <- out2
  }
  
  if(missing(ref)){
    o <- list(sim_communities = simComm)
  } else {
    o <- list(sim_communities = simComm,
              ref_communities = refComm)
  }
  
  return(o)
}

# function to generate simulated communities to add on ongoing restoration: ####
addSimulation <- function(trait, ava, und, it, rich, cwm, cwv, rao, cost, dens,
                          stan, ref, rest, max_add, min_p, phi = 1){
  # generates simulated communities and calculate its parameters.
  # ARGUMENTS:
  # trait: data frame or matrix with species traits. Traits as columns and species as rows.
  # ava: vector indicating availability of species (1 or 0)
  # und: vector indicating undesired species (1 or 0)
  # it: number of iterations (communities) for each area ongoing restoration
  # rich: range of richness values
  # cwm: vector with traits to calculate Community Weighted Mean (CWM). One CWM is calculated for each trait.
  # rao: vector with traits to calculate Rao Quadratic Entropy, or distance matrix (class dist)
  # cost: vector of species cost per individual
  # dens: vector of planting density
  # stan: which parameters should be standardized by the maximum? Must provide trait names
  # ref: matrix with proportion of species in the reference sites. NAs not accepted.
  # rest: matrix with proportion of species in the restored sites. NAs not accepted.
  # max_add: maximum proportion that can be added to restored areas.
  # min_p: minimum species proportion to consider that species occurs in the area
  
  # CREATE RANDOM COMMUNITIES #############################
  species <- row.names(trait)
  resample <- function(x, ...) x[sample.int(length(x), ...)]
  
  # Find distant species: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  if(!missing(rao)){
    if(inherits(rao, 'character')){
      t2d <- as.matrix(scale(trait[, rao, drop = FALSE]) )
    } else if(inherits(rao, 'dist')){
      t2d <- rao
    }
    if(!missing(cwm)){
      t2c <- as.matrix(trait[,cwm, drop = FALSE])
      constraints <- apply(t2c, 2, function(x){
        cons <- seq(min(x), max(x), length.out = 8)
        return(cons[c(-1,-8)])
      })
      selSpp <- vector(mode="list", length=length(cwm))
      names(selSpp) <- cwm
      for(i in cwm){
        cons_i <- constraints[,i]
        t2c_i <- t2c[,i, drop=FALSE]
        selSpp_i <- c()
        for(j in cons_i){
          names(j) <- colnames(t2c_i)
          invisible(capture.output(selSpp_j <- selectSpecies(t2c_i, j,
                                                             t2d, phi = phi) ) ) 
          selSpp_i <- cbind(selSpp_i, selSpp_j$prob)
        }
        selSpp[[i]] <- selSpp_i 
      }
      propMatrixSelSpp <- do.call(cbind, selSpp)
      sppMax <- c()
      propMin <- 0.1
      while(length(sppMax) < rich[1]){
        sppMax <- which(propMatrixSelSpp > propMin, arr.ind = TRUE)
        sppMax <- unique(species[sppMax[,1]])
        propMin <- 0.5*propMin
      }
      propMatrixSelSpp <- round(t(propMatrixSelSpp),3)
    } else {
      invisible( capture.output( selSpp <- selectSpecies(t2d = t2d, phi = phi)$prob ))
      propMatrixSelSpp <- selSpp
      sppMax <- c()
      propMin <- 0.1
      while(length(sppMax) < rich[1]){
        sppMax <- species[propMatrixSelSpp > propMin]
        propMin <- 0.5*propMin
      }
    }
    
    # Find distant species that are available: <<<<<<<<<<<<<<
    avaLog <- as.logical(ava)
    speciesAva <- species[avaLog]
    if(inherits(rao, 'character')){
      t2d <- as.matrix(scale(trait[avaLog, rao, drop = FALSE]) )
    } else if(inherits(rao, 'dist')){
      t2d <- as.dist(as.matrix(rao)[avaLog, avaLog])
    }
    if(!missing(cwm)){
      t2c <- as.matrix(trait[avaLog, cwm, drop = FALSE])
      constraints <- apply(t2c, 2, function(x){
        cons <- seq(min(x), max(x), length.out = 8)
        return(cons[c(-1,-8)])
      })
      selSppAva <- vector(mode="list", length=length(cwm))
      names(selSppAva) <- cwm
      for(i in cwm){
        cons_i <- constraints[,i]
        t2c_i <- t2c[,i, drop=FALSE]
        selSppAva_i <- c()
        for(j in cons_i){
          names(j) <- colnames(t2c_i)
          invisible(capture.output(selSppAva_j <- selectSpecies(t2c_i, j, t2d, phi = phi) ))
          selSppAva_i <- cbind(selSppAva_i, selSppAva_j$prob)
        }
        selSppAva[[i]] <- selSppAva_i 
      }
      propMatrixSelSppAva <- do.call(cbind, selSppAva)
      sppMaxAva <- c()
      propMin <- 0.1
      while(length(sppMaxAva) < rich[1]){ 
        sppMaxAva <- which(propMatrixSelSppAva > propMin, arr.ind = TRUE)
        sppMaxAva <- unique(speciesAva[sppMaxAva[,1]])
        propMin <- 0.5*propMin
      }
      propMatrixSelSppAva <- round(t(propMatrixSelSppAva),3)
    } else {
      invisible(capture.output(selSppAva <- selectSpecies(t2d = t2d, phi = phi)$prob ))
      propMatrixSelSppAva <- selSpp
      sppMaxAva <- c() 
      propMin <- 0.1 
      while(length(sppMaxAva) < rich[1]){
        sppMaxAva <- species[propMatrixSelSppAva > propMin]
        propMin <- 0.5*propMin
      }
    }
    
    # number of iterations for simulations: <<<<<<<<<<<<<<<<<
    itMax <- round(0.25*it)
    itMaxAva <- round(0.25*it)
    itAva <- round(0.25*it)
    itAll <- it - itMax - itMaxAva - itAva
    
    # Run simulation with species that maximize rao: <<<<<<<<<
    if(length(sppMax) < rich[2]){
      nsp <- length(sppMax)
      vLen <- nsp
    } else {
      nsp <- rich[2]
      vLen <- length(sppMax)
    }
    propMatrixSelSpp2 <- matrix(rep(0,length(species)*itMax),
                                ncol=length(species), nrow=itMax)
    sppMaxPos <- species %in% sppMax
    for(i in 1:itMax){
      nsp_i <-  resample(rich[1]:nsp, 1)
      ocor = sample( c(rep(1, nsp_i), rep(0, vLen - nsp_i)) )
      abund = rlnorm(vLen)
      abund <- abund * ocor
      prop <- abund/sum(abund)
      propMatrixSelSpp2[i,sppMaxPos] <- prop
    }
    
    # Run simulation with species that maximize rao and are available: <<<<
    if(length(sppMaxAva) < rich[2]){
      nsp <- length(sppMaxAva) 
      vLen <- nsp
    } else {
      nsp <- rich[2]
      vLen <- length(sppMaxAva)
    }
    propMatrixSelSppAva2 <- matrix(rep(0,length(species)*itMaxAva),
                                   ncol=length(species), nrow=itMaxAva)
    sppMaxAvaPos <- species %in% sppMaxAva
    for(i in 1:itMaxAva){
      nsp_i <-  resample(rich[1]:nsp, 1)
      ocor = sample( c(rep(1, nsp_i), rep(0, vLen - nsp_i)) )
      abund = rlnorm(vLen)
      abund <- abund * ocor
      prop <- abund/sum(abund)
      propMatrixSelSppAva2[i,sppMaxAvaPos] <- prop
    }
  } else {
    avaLog <- as.logical(ava)
    itAva <- round(0.5*it)
    itAll <- it - itAva
  }
  
  # Run simulation with available species: <<<<<<<<<<<<<<<<
  if(sum(ava) < rich[2]){
    nsp <- sum(ava) 
    vLen <- nsp
  } else {
    nsp <- rich[2]
    vLen <- sum(ava)
  }
  propMatrixAva <- matrix(rep(0,length(species)*itAva),
                          ncol=length(species), nrow=itAva)
  if(rich[1] > nsp){
    stop("Minimum richness is higher than number of available species")
  }
  for(i in 1:itAva){
    nsp_i <-  resample(c(rich[1]:nsp), 1)
    ocor = sample( c(rep(1, nsp_i), rep(0, vLen - nsp_i)) )
    abund = rlnorm(vLen)
    abund <- abund * ocor
    prop <- abund/sum(abund)
    propMatrixAva[i,avaLog] <- prop
  }
  
  # Run simulation with all species: <<<<<<<<<<<<<<<<<<<<<<
  nsp <- length(species)
  propMatrixPool <- matrix(rep(0,length(species)*itAll),
                           ncol=length(species), nrow=itAll)
  for(i in 1:itAll){
    nsp_i <-  resample(c(rich[1]:rich[2]), 1)
    ocor = sample( c(rep(1, nsp_i), rep(0, nsp - nsp_i)) )
    abund = abund = rlnorm(nsp)
    abund <- abund * ocor
    prop <- abund/sum(abund)
    propMatrixPool[i,] <- prop
  }
  
  # Bind all matrices: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  if(!missing(rao)){
    propMatrix <- rbind(propMatrixSelSpp2,propMatrixSelSppAva2,
                        propMatrixAva, propMatrixPool)
  } else {
    propMatrix <- rbind(propMatrixAva, propMatrixPool)
  }
  
  rownames(propMatrix) <- sprintf("sim%d",seq(1:nrow(propMatrix)))
  colnames(propMatrix) <- species
  
  # Set prop = 0 to undesired species: <<<<<<<<<<<<<<<<<<<<
  if(!missing(und)){
    pos <- as.logical(und)
    propMatrix[,pos] <- 0
    pos <- rowSums(propMatrix) == 0 #se ==0 simulacao contem apenas und
    propMatrix <- propMatrix[!pos,] #elimina linhas sem spp
    propMatrix <- propMatrix/rowSums(propMatrix)
  }
  
  row.names(propMatrix) <- sprintf("sim%d",seq(1:nrow(propMatrix)))
  
  # TRANSFORM PROPORTIONS AND SUM TO REST #################
  propMatrixAdd <- propMatrix * max_add #transforma matriz
  #Set prop = 0 to rare species: <<<<<<<<<<<<<<<<<<<<<<<<<
  if(!missing(min_p)){
    pos <- propMatrixAdd < min_p
    propMatrixAdd[pos] <- 0
    propMatrixAdd <- (propMatrixAdd/rowSums(propMatrixAdd)) * max_add
    pos <- is.na(propMatrixAdd[,1])
    propMatrixAdd <- propMatrixAdd[!pos,]
  }
  propMatrixList <- apply(rest, 1, FUN=function(x){ #para cada comun restaurada
    propMatrix_x <- apply(propMatrixAdd, 1, FUN = function(y){ #para cada comun simulada
      x_y <- x + y #restaurada + simulada
      return(x_y)
    })
    propMatrix_x <- t(propMatrix_x)
    propMatrix_x <- propMatrix_x/rowSums(propMatrix_x)
    return(propMatrix_x)
  }, simplify = FALSE)
  propMatrixTab <- do.call(rbind, propMatrixList)
  
  # CALCULATE PARAMETERS ##################################
  if(!require(adiv)){
    stop("Package adiv not found")
  }
  if(!require(fundiversity)){
    stop("Package fundiversity not found")
  }
  
  out <- NULL
  
  #rest <- rest/rowSums(rest)
  x <- rbind(propMatrixTab, rest)
  
  if(!missing(ref)){
    ref <- ref/rowSums(ref)
    x <- rbind(x, ref)
  } else {x <- x}
  
  UNA <- apply(x, 1, FUN = function(a) sum(a[!as.logical(ava)] > 0) )
  out <- cbind(out, unavailable = UNA)
  
  S <- apply(x, 1, FUN = function(a) sum(a > 0))
  out <- cbind(out, richness = S)
  
  if(!missing(cwm)){
    if(inherits(cwm, 'character')){
      traitSub <- trait[,cwm, drop=FALSE]
      CWM <- apply(x, 1, FUN=function(p){
        colSums(traitSub*p, na.rm = T)
      })
      if(class(CWM)[1] == 'matrix'){
        CWM <- t(CWM)
      } else {
        CWM <- as.matrix(CWM)
      }
      #colnames(CWM) <- paste0('cwm_', cwm)
      out <- cbind(out, CWM)
    } else {message('** CWM skipped. **')}
  }
  
  if(!missing(cwv)){
    if(inherits(cwv, 'character')){
      traitSub <- trait[,cwv, drop=FALSE]
      CWM <- apply(x, 1, FUN=function(p){
        colSums(traitSub*p, na.rm = T)
      })
      if(class(CWM)[1] == 'matrix'){
        CWM <- t(CWM)
      } else {
        CWM <- as.matrix(CWM)
      }
      
      CWV <- NULL
      for(i in 1:nrow(CWM)){ #para cada linha de CWMs
        CWM_i <- CWM[i,] #linha i
        traitMod <- NULL #traits modificados pelo CWM i
        for(j in 1:ncol(traitSub)){
          traitMod_j <- traitSub[,j] - CWM_i[j] #subtrai pelo CWM
          traitMod_j <- traitMod_j^2 #eleva ao quadrado
          traitMod <- cbind(traitMod, traitMod_j) #adiciona na tabela
        }
        #CWV: traits modificados pelo CWM i vezes proporcao i
        CWV_i <- colSums(traitMod*x[i,])
        CWV <- rbind(CWV, CWV_i)
      }
      #colnames(CWV) <- paste0('cwv_', cwv)
      colnames(CWV) <- cwv
      #CWV <- apply(CWV, 2, FUN = function(x){x/max(x)}) #padroniza pelo maximo #$$$$$$$
      out <- cbind(out, CWV)
    } else {message('** CWV skipped. **')}
  }
  
  if(!missing(rao)){
    if(inherits(rao, 'character')){
      traitSub <- scale(trait[,rao, drop=FALSE] )
      RAO <- fd_raoq(traitSub, x)$Q
    } else if(inherits(rao, 'dist')){
      RAO <- fd_raoq(sp_com = x, dist_matrix = rao)$Q
    } else{
      message('** RAO skipped')
    }
    #RAO <- RAO/max(RAO) #$$$$$$$$$$$$
    out <- cbind(out, rao = RAO)
    colnames(out)[colnames(out) == 'rao'] <- rao
  }
  
  if(!missing(cost)){
    COST <- apply(x, 1, FUN=function(p){
      COST_i <- sum(p*cost*dens, na.rm = TRUE)
      return(COST_i)
    })
    out <- cbind(out, cost = COST)
  }
  
  if(!missing(stan)){
    out[,stan] <- out[,stan,drop=F]/max(out[,stan,drop=F])
  } #$$$$$$$$$$$$$$$$$$$
  
  if(!missing(ref)){
    f <- c(rep(row.names(rest),each = nrow(propMatrixList$X1) ), #posicao comunidades restauradas modificadas
           rep('restored', nrow(rest)), #posicao comunidades restauradas
           rep('reference', nrow(ref)) #posicao comunidades ref
    )
  } else {
    f <- c(rep(row.names(rest),each = nrow(propMatrixList$X1)  ), #posicao comunidades restauradas modificadas
           rep('restored', nrow(rest)), #posicao comunidades restauradas
    )
  } #https://stackoverflow.com/questions/13060000/how-can-i-separate-a-matrix-into-smaller-ones-in-r
  
  pos <- duplicated(colnames(out)) | duplicated(colnames(out), fromLast=TRUE)
  dup_col <- colnames(out)[pos]
  correct_col <- paste0(dup_col, 1:length(dup_col)) #add index to duplicated
  colnames(out)[pos] <- correct_col
  
  outList <- lapply( split( out, f ), matrix, ncol=ncol(out))
  for(i in 1:length(outList)){
    colnames(outList[[i]]) <- colnames(out)
    row.names(outList[[i]]) <- sprintf("sim%d",seq(1:nrow(outList[[i]])))
  }
  
  if(!missing(ref)){
    propMatrixList2 <- c(propMatrixList, list(restored = rest), list(reference = ref) )
  } else {
    propMatrixList2 <- c(propMatrixList, list(restored = rest) )
  }
  
  ordem <- match(names(propMatrixList2), names(outList))
  outList <- outList[ordem]
  all(names(outList) == names(propMatrixList2))
  
  out2 <- list(compositions = propMatrixList2, parameters = outList)
  
  return(out2)
}

# function to calculate multifunctionality: #############
multCalculation <- function(x, th, bel){
  #x: list of table with parameters of different areas
  #th: named vector of thresholds. Names must match x columns
  #bel: character vector indicating functions that must be below threshold. Must match x columns.
  funs <- names(th) #functions
  if(!missing(bel)){
    th[bel] <- -th[bel]
    bel2 <- bel #para usar no apply
    } #reflect functions
  x2_nrow <- NULL
  e <- environment() #para buscar bel2, abaixo
  x2 <- lapply(x, FUN = function(y){ #for each table:
    if(exists('bel2', envir = e, inherits = FALSE)){
      y[,bel2] <- -y[,bel2]
      } #reflect functions. Nota sobre exists: procura no environment definido em 'envir', no caso, environment da funcao multCalculation. Porem, considera que 'bel' existe mesmo se nao tiver valor atribuido, por que aparece no environment como 'missing argument'. Por isso, preciso salvar bel em bel2, acima. missing() ou hasArg() nao funcionam pq buscam apenas no environment da funcao anonima FUN. 
    x3 <- apply(y, 1, FUN = function(y2){ #for each line
      tests <- as.numeric(y2[funs] > th) #is function above threshold?
      mult <- sum(tests) #sum of functions above threshold
      return(c(tests, mult)) #output
    })
    x3 <- t(x3) #transpose
    colnames(x3) <- c(funs, 'multifunctionality') #rename
    return(x3)
    x2_nrow <- c(x2_nrow, nrow(x3))
  })
  
  return(x2)
  
}

# function to calculate functional dissimilarity: #########
disCalculation <- function(sim, trait, dis, ref){
  # Calculate functional dissimilarity between simulated communities and reference sites
  # ARGUMENTS:
  # sim: data frame or matrix with species proportions in simulated communities. Species as columns and simulated communities as rows.
  # trait: data frame or matrix with species traits. Traits as columns and species as rows.
  # dis: matrix of trait distances
  # ref: data frame or matrix with species proportions in reference sites. Species as columns and simulated communities as rows.
  
  if(inherits(trait, 'data.frame') | inherits(trait, 'matrix')){
    dis <- dist(scale(trait))
    message('*** Calculating functional dissimilarity. This may take a while. ***')
  } else if(inherits(trait, 'dist')){
    dis <- trait
    message('*** Calculating functional dissimilarity. This may take a while. ***')
  }
  i = 0
  DISSIM <- apply(ref, 1, FUN=function(r){
    i <<-  i+1 
    message(paste('##### Reference site number: ',i, '#####') )
    j = 1
    pb <- txtProgressBar(min = 0, max = nrow(sim), style = 3)
    DISSIM_i <- apply(sim, 1, FUN=function(p){ 
      setTxtProgressBar(pb, j)
      j <<- j +1
      comm <- rbind(r, p)
      DISSIM_j <- discomQE(comm, dis, formula = "QE")
      return(DISSIM_j)
    })
    close(pb)
    return(DISSIM_i)
  })
  DISSIM <- apply(DISSIM, 1, mean)
  DISSIM <- DISSIM/max(DISSIM)
  
  return(DISSIM)
  
}

# function to select communities: #########################
comSelection <- function(param, comp, tests){
  #select simulated communities based on tests provided and add information to x
  #ARGUMENTS:
  #param: data.frame with parameters of simulated communities
  #comp: data.frame with compositions of simulated communities
  #tests: list with tests to be performed
  
  xPar <- param
  completeString <- paste0('xPar', '$', tests)
  testsEval <- sapply(completeString, function(a) eval(parse(text=a)))
  pos <- apply(testsEval, 1, all) 
  selPar <- xPar[pos,] 
  selCom <- comp[pos,]
  
  #number of selected communities:
  nSel <- apply(testsEval, 2, sum)
  names(nSel) <- tests
  nSel <- c(nSel, all = sum(pos))
  
  #thresholds:
  testsSplit <- strsplit(tests, ' ')
  trsh <- as.numeric(sapply(testsSplit, '[', 3))
  names(trsh) <- sapply(testsSplit, '[', 1)
  
  outSel <- list(parameters = selPar,
                 composition = selCom,
                 N = nSel,
                 thresholds = trsh)
  
  return(outSel)
  
}

# function to visualize results: ##########################
viewResults <- function(x, y, sim,
                        cols = c('grey', 'black', brewer.pal(3, 'Paired')),
                        xlab, ylab,
                        legend,
                        xlim = range(sim$sim_communities$parameters[,x]),
                        ylim = range(sim$sim_communities$parameters[,y]),
                        hide_notsel = FALSE, hide_una = FALSE,
                        hide_ref = FALSE, hide_sel = FALSE){
  # Visualize parameters of simulated and selected communities and reference sites.
  # ARGUMENTS:
  # x: name of the variable in x axis. Must match a column in sim and sel_sim
  # y: name of the variable in the y axis. Must match a column in sim and sel_sim
  # sim: list with results of function comSelection.
  # cols: colors of simulated communities, in this order: communities with at least one unavailable species; communities with only available species; communities with at least one unavailable species that satisfy restoration criteria; communities with only available species that satisfy restoration criteria.
  # xlab: x axis label
  # ylab: y axis label
  # xlim: x axis limits
  # ylim: y axis limits
  # hide_una: hide simulated communities with unavailable species?
  # hide_notsel: hide simulated communities not selected?
  # hide_ref: hide reference sites?
  
  if(!require(RColorBrewer)){
    stop("Package RColorBrewer not found")
  }
  
  all <- sim$sim_communities$parameters
  if(!hide_ref){
    ref <- sim$ref_communities$parameters
  }
  if(!hide_sel){
    sel <- sim$sel_communities$parameters 
  }
  
  # Define xlim e ylim:
  if(missing(xlim)){
    xlim <- range(all[,x])
  }
  if(missing(ylim)){
    ylim <- range(all[,y])
  }
  
  par(mfrow=c(1,2), mar=c(5.1, 1, 4.1, 0), oma=c(0,3,0,0))
  plot(0,0, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
       col=rgb(0,0,0,alpha=0))
  mtext(ylab,2,3)
  
  #ALL:
  if(!hide_notsel & !hide_una){
    points(all[,x], all[,y], pch = 19, col=cols[1])
  }
  
  #AVA NOT SEL:
  if(!hide_notsel){
    pos <- all$unavailable==0 #posicao de ava
    points(all[pos,x], all[pos,y], pch=19, col=cols[2])
  } 
  
  #SELECTED:
  if(!hide_sel){
    if(!hide_una){
      points(sel[,x], sel[,y], pch=19, col=cols[3])
      pos <- sel$unavailable == 0 
      points(sel[pos,x], sel[pos,y], pch = 19, col=cols[4])
    } else {
      pos <- sel$unavailable == 0
      points(sel[pos,x], sel[pos,y], pch = 19, col=cols[4])
    }
  }
  
  #REF:
  if(!hide_ref){
    points(ref[,x], ref[,y], pch=19, col=cols[5], cex=1.5)
  }
  
  tsh <- sim$thresholds
  abline(v = tsh[x], lty=2, col=cols[4])
  abline(h = tsh[y], lty=2, col=cols[4])
  
  legend <- c("Unavailable",
              "Available",
              "Unavailable - selected",
              "Available - selected",
              "References")
  pos <- c(!hide_una & !hide_notsel,
           !hide_notsel,
           !hide_una & !hide_sel,
           !hide_sel,
           !hide_ref)
  
  plot.new()
  legend("topleft", col = cols[pos], cex = 1, pch=19,
         legend=legend[pos])
}

viewResults2 <- function(x, y, sim,
                         cols = c('grey', 'black', brewer.pal(3, 'Paired')),
                         xlab, ylab,
                         legend,
                         xlim = range(sim$sim_communities$parameters[,x]),
                         ylim = range(sim$sim_communities$parameters[,y]),
                         hide_notsel = FALSE, hide_una = FALSE,
                         hide_ref = FALSE, hide_sel = FALSE){
  # Visualize parameters of simulated and selected communities and reference sites.
  # ARGUMENTS:
  # x: name of the variable in x axis. Must match a column in sim and sel_sim
  # y: name of the variable in the y axis. Must match a column in sim and sel_sim
  # sim: list with results of function comSelection.
  # cols: colors of simulated communities, in this order: communities with at least one unavailable species; communities with only available species; communities with at least one unavailable species that satisfy restoration criteria; communities with only available species that satisfy restoration criteria.
  # xlab: x axis label
  # ylab: y axis label
  # xlim: x axis limits
  # ylim: y axis limits
  # hide_una: hide simulated communities with unavailable species?
  # hide_notsel: hide simulated communities not selected?
  # hide_ref: hide reference sites?
  
  if(!require(RColorBrewer)){
    stop("Package RColorBrewer not found")
  }
  
  all <- sim$sim_communities$parameters
  if(!hide_ref){
    ref <- sim$ref_communities$parameters
  }
  if(!hide_sel){
    sel <- sim$sel_communities$parameters 
  }
  
  # Define xlim e ylim:
  if(missing(xlim)){
    xlim <- range(all[,x])
  }
  if(missing(ylim)){
    ylim <- range(all[,y])
  }
  
  #par(mfrow=c(1,2), mar=c(5.1, 1, 4.1, 0), oma=c(0,3,0,0))
  plot(0,0, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
       col=rgb(0,0,0,alpha=0))
  #mtext(ylab,2,3)
  
  #ALL:
  if(!hide_notsel & !hide_una){
    points(all[,x], all[,y], pch = 19, col=cols[1])
  }
  
  #AVA NOT SEL:
  if(!hide_notsel){
    pos <- all$unavailable==0 #posicao de ava
    points(all[pos,x], all[pos,y], pch=19, col=cols[2])
  } 
  
  #SELECTED:
  if(!hide_sel){
    if(!hide_una){
      points(sel[,x], sel[,y], pch=19, col=cols[3])
      pos <- sel$unavailable == 0 
      points(sel[pos,x], sel[pos,y], pch = 19, col=cols[4])
    } else {
      pos <- sel$unavailable == 0
      points(sel[pos,x], sel[pos,y], pch = 19, col=cols[4])
    }
  }
  
  #REF:
  if(!hide_ref){
    points(ref[,x], ref[,y], pch=19, col=cols[5], cex=1.5)
  }
  
  tsh <- sim$thresholds
  abline(v = tsh[x], lty=2, col=cols[4])
  abline(h = tsh[y], lty=2, col=cols[4])
  
  legend <- c("Unavailable",
              "Available",
              "Unavailable - selected",
              "Available - selected",
              "References")
  pos <- c(!hide_una & !hide_notsel,
           !hide_notsel,
           !hide_una & !hide_sel,
           !hide_sel,
           !hide_ref)
  
  #plot.new()
  # legend("topleft", col = cols[pos], cex = 1, pch=19,
  #        legend=legend[pos])
}
