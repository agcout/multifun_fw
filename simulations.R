# +++++++++++++++++++++++++++++++++++++++++++++++
# Natural regeneration enhances ecosystem
# multifunctionality but species addition can
# increase it during restoration monitoring
# Coutinho et al.
# +++++++++++++++++++++++++++++++++++++++++++++++

max_add_seq <- seq(0.1,2,0.1) #vector of densities
set.seed(75421) #seed to define seeds
seed_seq <- sample(1:1000, 20) #seeds
run_up_to <- c() #to track completed iterations

#################################################
# > FUNCTIONS AND PACKAGES ######################
#++++++++++++++++++++++++++++++++++++++++++++++++

library(adiv)
library(Select)
library(fundiversity)
library(RColorBrewer)
library(vegan)
library(ComplexUpset)
library(ggplot2)
library(ape)
library(prettyGraphs) 
library(tidyr)

source('framework_functions.R')
library(ecodist)

#################################################
# > DATA ########################################
#++++++++++++++++++++++++++++++++++++++++++++++++

data <- read.csv2('data_secil.csv')
bur_sit <- read.csv2('burned_sites.csv')
all(bur_sit$site == colnames(data)[12:89]) #checking

# Identify undesired species:
und <- data$Undesired

# Identify reference areas:
your_reference_sites <- c("X1VN", "X3.1VN", "X3.2VN", "X4VN", "X5.1VN", 
                          "X5.2VN", "X6VN", "CTM1", "CTM2", "CTM3", 
                          "CTM4", "CTM5", "CTM6", "UECN1", "UECN2",
                          "UECN3", "UEMA1", "UEMA2","UEMA3")
setdiff(your_reference_sites, colnames(data))
REF <- t(data[,your_reference_sites] )
colnames(REF) <- data$all_species

# Identify restored areas:
your_restoration_sites <- c("X1", "X3", "X4", "X6", "X7", "X8", "X9",
                            "X10", "X11", "X13", "X14", "X16", "X17", 
                            "X19", "X20", "X21", "X22", "X27", "X28", 
                            "X29", "X30", "X36", "X37", "X38", "X39",
                            "X43", "X44", "X45", "X46", "X47", "X48", 
                            "X51", "X52", "X53", "X55", "X56", "X57",
                            "X58", "X60", "X61", "X62", "X63", "X65", 
                            "X69", "X71", "X74", "X75", "X76", "X79", 
                            "X80", "PE10", "PE11", "PE13.1", "PE13.2",
                            "PE14", "PE20", "UECR1", "UECR2", "UECR3")

setdiff(your_restoration_sites, colnames(data))
REST <- t(data[,your_restoration_sites] )
colnames(REST) <- data$all_species

# List of species on each set:
spp_all <- data$all_species
ab_sum <- colSums(REF) #total abundance of each species in REF
spp_ref <- spp_all[ab_sum > 0] #species that occur in any reference site
ab_sum <- colSums(REST) #total abundance of each species in REST
spp_rest <- spp_all[ab_sum > 0] #species that occur in any restored site
pos <- !spp_all%in%spp_ref & !spp_all%in%spp_rest
spp_arr <- spp_all[pos] #species of Arrabida Park, but not in restored area
spp_plan <- data$all_species[data$Planted == 1] #planted species
spp_spon <- spp_rest[!spp_rest%in%spp_plan] #spontaneous species
spp_sets <- list(pool = spp_all, ref = spp_ref,
                 rest = spp_rest, arrab = spp_arr,
                 plan = spp_plan, spon = spp_spon
                 )
saveRDS(spp_sets, 'outputs/species_sets/spp_sets.rds')

# select the columns in data with the functional traits
your_traits <- colnames(data[,c('LMA', 'Height', 'Dur_flowering',
                                'Entomophily', 'Zooc', 'Resprouter')])
traits <- data[,c(your_traits)]
row.names(traits) <- data$all_species #rename rows

#################################################
# > SIMULATIONS #################################
#++++++++++++++++++++++++++++++++++++++++++++++++

# >> Ongoing restoration area ###################
rounds <- c(1:length(max_add_seq)) #define which rounds to run
for(SIMULATION_INDEX in rounds){
  INI <- Sys.time()
  
  set.seed(seed_seq[SIMULATION_INDEX])
  
  ###############################################
  # > SIMULATION ################################
  #++++++++++++++++++++++++++++++++++++++++++++++
  
  # >> with all species ###########################
  SR<- c(10, 26) #alfa richness in reference sites
  N <- 10000
  cwm	<- c('LMA', 'Dur_flowering',
           'Entomophily', 'Zooc', 'Resprouter')
  cwv <-  'Height'
  ava <- data$Planted
  max_add <- max_add_seq[SIMULATION_INDEX]
  
  ini <- Sys.time()
  addSim <- addSimulation(trait = traits, ava = ava,
                          und = und,
                          it = N, rich = SR,
                          cwm = cwm, cwv = cwv,
                          ref = REF, rest = REST,
                          max_add = max_add,
                          min_p = 0.01)
  fim <- Sys.time()
  (dif <- fim - ini)
  FILENAME <- paste0('outputs/',SIMULATION_INDEX,'/','addSim.rds')
  saveRDS(addSim, FILENAME)
  addSim <- readRDS(FILENAME)
  
  rownames(addSim$compositions$restored) <- your_restoration_sites
  rownames(addSim$parameters$restored) <- your_restoration_sites
  rownames(addSim$compositions$reference) <- your_reference_sites
  rownames(addSim$parameters$reference) <- your_reference_sites
  
  #richness range in each site:
  sapply(addSim$parameters, FUN=function(x){range(x[,2])})
  
  #Round relative abundances close to zero:
  indexes <- which(!names(addSim$compositions) %in% c('restored', 'reference'))
  for(i in indexes){
    comp <- addSim$compositions[[i]] #composition matrix i
    comp[comp <= 0.01] <- 0 #replace by 0 if <=0.01
    totals <- rowSums(comp) #total in each simulation
    comp <- comp/totals #standardize
    addSim$compositions[[i]] <- comp
  }
  
  # >> with resprouters species #######################
  # Simulation with resprouters (only).
  und_nr <- as.logical(data$Resprouter == 0) #position of not resprouter
  all(data$all_species == colnames(REST)) #checking
  ini <- Sys.time()
  addSim_r <- addSimulation(trait = traits, ava = ava,
                            und = und_nr,
                            it = N, rich = SR,
                            cwm = cwm, cwv = cwv,
                            ref = REF, rest = REST,
                            max_add = max_add,
                            min_p = 0.01)
  fim <- Sys.time()
  (dif <- fim - ini)
  FILENAME <- paste0('outputs/',SIMULATION_INDEX,'/','addSim_r.rds')
  saveRDS(addSim_r, FILENAME)
  addSim_r <- readRDS(FILENAME)
  
  #Round relative abundances close to zero:
  indexes <- which(!names(addSim_r$compositions) %in% c('restored', 'reference'))
  for(i in indexes){
    comp <- addSim_r$compositions[[i]] #composition matrix i
    comp[comp <= 0.01] <- 0 #replace by 0 if <=0.01
    totals <- rowSums(comp) #total in each simulation
    comp <- comp/totals #standardize
    table(rowSums(comp)) #checking
    addSim_r$compositions[[i]] <- comp
  }
  
  # >> Replace simulations ###########################
  # In burned sites, replace original simulations by new simulations with resprouters only.
  bur_sit2 <- bur_sit$site[as.logical(bur_sit$fire)] #burned sites
  all(names(addSim$compositions)[1:59] == rownames(addSim$compositions$restored)) #checking
  is <- which(names(addSim$compositions) %in% bur_sit2) #indices
  #replace compositions:
  addSim$compositions[is] <- addSim_r$compositions[is]
  #replace parameters:
  addSim$parameters[is] <- addSim_r$parameters[is]
  
  # Get maximum values of CWV and standardize. Necessary to standardize simulation without spontaneous.
  #maximos:
  maximos <- sapply(addSim$parameters, FUN = function(x){
    maximos_x <- apply(x[,'Height',drop=F], 2, max)
    return(maximos_x)
  }) #maximum in each simulation
  if(class(maximos)[1] == 'numeric'){
    maximos <- as.matrix(maximos)
  }
  maximo <- apply(maximos, 2, max) #maximo of all simulations
  #standardize:
  for(i in 1:length(addSim$parameters)){
    param_i <- addSim$parameters[[i]][,'Height',drop=F] #parametros site i
    new_val <- mapply('/', as.data.frame(param_i), maximo) #divide pelo maximo
    addSim$parameters[[i]][,'Height'] <- new_val
  }
  
  # >> with planted species #######################
  # Simulation is not necessary. We need only the parameters
  #without spontaneous. addSimulation is used, but only the
  #output of restored sites is used.
  und_s <- as.logical(data$Planted != 1) #position of spontaneous
  all(data$all_species == colnames(REST)) #checagem
  REST_p <- REST #RESTored planted only
  REST_p[,und_s] <- 0 #valor nas colunas undesired = 0
  ini <- Sys.time()
  addSim_p <- addSimulation(trait = traits, ava = ava,
                            und = und,
                            it = 100, rich = SR,
                            cwm = cwm, cwv = cwv,
                            ref = REF, rest = REST_p, #excluir und_s das areas existentes
                            max_add = max_add,
                            min_p = 0.01)
  fim <- Sys.time()
  (dif <- fim - ini)
  FILENAME <- paste0('outputs/',SIMULATION_INDEX,'/','addSim_p.rds')
  saveRDS(addSim_p, FILENAME)
  addSim_p <- readRDS(FILENAME)
  
  p_only <- list()
  p_only$compositions <- addSim_p$compositions$restored
  p_only$parameters <- addSim_p$parameters$restored

  p_only$parameters[,"Height"] <- p_only$parameters[,"Height"]/maximo
  
  # >> Exclude 'Height1' column #####################
  # # Height1 correspond to the CWM of height, which is
  # #used to calculate CWV, but is not used in this study
  # pos1 <- colnames(addSim$parameters$X1) == 'Height1'
  # pos2 <- colnames(addSim$parameters$X1) == 'Height2'
  # 
  # #addSim:
  # sapply(addSim$parameters, ncol)
  # for(i in 1:length(addSim$parameters)){
  #   colnames(addSim$parameters[[i]])[pos2] <- 'Height'
  #   addSim$parameters[[i]] <- addSim$parameters[[i]][,!pos1]
  # }
  # sapply(addSim$parameters, ncol)
  # 
  # #p_only:
  # colnames(p_only$parameters)[pos2] <- 'Height'
  # p_only$parameters <- p_only$parameters[,!pos1]
  # ncol(p_only$parameters)
  # addSim_p$parameters$restored <- p_only$parameters
  
  #################################################
  # > MULTIFUNCTIONALITY CALCULATION ##############
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  # Define thresholds: <<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ref_fun <- addSim$parameters$reference[,c(-1,-2)] #reference function columns
  
  # Threshold (mean value of parameter):
  thrs <- apply(ref_fun, 2, mean)
  saveRDS(thrs, 'outputs/thresholds.rds')
  
  # Multifunctionality calculation: <<<<<<<<<<<<<<<<<<<<<<<<<<
  mult <- multCalculation(x = addSim$parameters, th = thrs)
  addSim$multifunctionality <- mult
  
  mult_p <- multCalculation(x = addSim_p$parameters, th = thrs)
  mult_p$restored[29,] <- 0 #community 29 has no planted species
  addSim_p$multifunctionality <- mult_p
  p_only$multifunctionality <- mult_p$restored
  
  # save data: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  FILENAME <- paste0('outputs/', SIMULATION_INDEX, '/', 'addSim_multifunctionality_restored.txt')
  write.table(addSim$multifunctionality$restored, FILENAME)
  FILENAME <- paste0('outputs/', SIMULATION_INDEX, '/', 'addSim_multifunctionality_reference.txt')
  write.table(addSim$multifunctionality$reference, FILENAME)
  FILENAME <- paste0('outputs/', SIMULATION_INDEX,'/','addSim_p_multifunctionality_restored.txt')
  write.table(addSim_p$multifunctionality$restored, FILENAME)
  
  #################################################
  # > SIMULATION SELECTION ###############
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  # This section is divide in two parts:
  # In unburned areas, the only criteria is higher multifunctionality.
  # In burned areas, the first criteria is higher % of resprouters the second
  #is higher multifunctionality.
  
  # >> with all species: ##########################
  # >>> selection of burned sites: ################
  bur_sit2 <- bur_sit$site[as.logical(bur_sit$fire)] #burned sites
  all(names(addSim$compositions)[1:59] == rownames(addSim$compositions$restored)) #checking
  
  comp_sel_bur <- NULL
  mult_sel_bur <- NULL
  param_sel_bur <- NULL
  index_sel_bur <- NULL #index of selected simulations. For tracking.
  is <- which(names(addSim$compositions) %in% bur_sit2) #indices
  for(i in is){
    #subset tables of ith site:
    comp_i <- addSim$compositions[[i]]
    param_i <- addSim$parameters[[i]]
    mult_i <- addSim$multifunctionality[[i]]
    site <- names(addSim$compositions)[i]
    #add values of existent community:
    comp_i <- rbind(comp_i, addSim$compositions$restored[i,])
    param_i <- rbind(param_i, addSim$parameters$restored[i,])
    mult_i <- rbind(mult_i, addSim$multifunctionality$restored[i,])
    #first filter:
    pos <- mult_i[,"Resprouter"] == 1 #restore resprouter
    if(sum(pos) > 0){
      comp_i <- comp_i[pos,, drop = FALSE]
      mult_i <- mult_i[pos,, drop = FALSE]
      param_i <- param_i[pos,, drop = FALSE]
      index_i <- which(pos)
      
      #second filter:
      pos <- which.max(mult_i[,"multifunctionality"]) #higher MF value
      comp_i <- comp_i[pos,, drop = FALSE]
      mult_i <- mult_i[pos,, drop = FALSE]
      param_i <- param_i[pos,, drop = FALSE]
      index_i <- index_i[pos]
    } else {
      #first filter:
      pos <- which.max(param_i[,"Resprouter"]) #higher % of resprouter
      comp_i <- comp_i[pos,, drop = FALSE]
      mult_i <- mult_i[pos,, drop = FALSE]
      param_i <- param_i[pos,, drop = FALSE]
      index_i <- pos
    }
    #rename row:
    row.names(comp_i) <- site
    row.names(param_i) <- site
    row.names(mult_i) <- site
    names(index_i) <- site
    #add to results:
    comp_sel_bur <- rbind(comp_sel_bur, comp_i)
    mult_sel_bur <- rbind(mult_sel_bur, mult_i)
    param_sel_bur <- rbind(param_sel_bur, param_i)
    index_sel_bur <- c(index_sel_bur, index_i)
  }
  nrow(comp_sel_bur); nrow(mult_sel_bur); nrow(param_sel_bur)
  sum(bur_sit$fire)
  
  resprouter_is <- addSim$parameters$restored[is,"Resprouter"]
  hist(resprouter_is, xlab = 'Resprouter %', main = '')
  
  # >>> selection of remaining sites: #############
  comp_sel <- NULL
  mult_sel <- NULL
  param_sel <- NULL
  index_sel <- NULL
  length(mult)
  final <- length(mult) - 2 #last indeces are restored and reference sites
  is2 <- 1:final
  is2 <- setdiff(is2, is) #all indices, but restored and reference
  for(i in is2){
    comp_i <- addSim$compositions[[i]]
    param_i <- addSim$parameters[[i]]
    mult_i <- addSim$multifunctionality[[i]]
    site <- names(addSim$compositions)[i]
    #add values of existent community:
    comp_i <- rbind(comp_i, addSim$compositions$restored[i,])
    param_i <- rbind(param_i, addSim$parameters$restored[i,])
    mult_i <- rbind(mult_i, addSim$multifunctionality$restored[i,])
    #first filter:
    pos <- mult_i[,"multifunctionality"] == max(mult_i[,"multifunctionality"])
    comp_i <- comp_i[pos,, drop = FALSE]
    mult_i <- mult_i[pos,, drop = FALSE]
    param_i <- param_i[pos,, drop = FALSE]
    index_i <- which(pos)
    #second filter:
    pos <- sample(nrow(param_i))[1] #criteria: random draw
    comp_i <- comp_i[pos,, drop = FALSE]
    mult_i <- mult_i[pos,, drop = FALSE]
    param_i <- param_i[pos,, drop = FALSE]
    index_i <- index_i[pos]
    #rename row:
    row.names(comp_i) <- site
    row.names(param_i) <- site
    row.names(mult_i) <- site
    names(index_i) <- site
    #add to results:
    comp_sel <- rbind(comp_sel, comp_i)
    mult_sel <- rbind(mult_sel, mult_i)
    param_sel <- rbind(param_sel, param_i)
    index_sel <- c(index_sel, index_i)
  }
  #indexsel1 = index_sel
  nrow(comp_sel); nrow(mult_sel); nrow(param_sel)
  
  # >>> bind selected simulations: ################
  #bind:
  comp_sel <- rbind(comp_sel, comp_sel_bur)
  mult_sel <- rbind(mult_sel, mult_sel_bur)
  param_sel <- rbind(param_sel, param_sel_bur)
  index_sel <- c(index_sel, index_sel_bur)
  
  #reorder:
  comp_sel <- comp_sel[names(addSim$compositions)[1:59],]
  param_sel <- param_sel[names(addSim$compositions)[1:59],]
  mult_sel <- mult_sel[names(addSim$compositions)[1:59],]
  index_sel <- index_sel[names(addSim$compositions)[1:59]]
  
  #checking:
  row.names(comp_sel)
  row.names(param_sel)
  row.names(mult_sel)
  names(index_sel)
  
  addSimSel <- list(compositions = comp_sel,
                    parameters = param_sel,
                    multifunctionality = mult_sel)
  
  # >>> calculate species richness: #############
  #how many new species in each simulation? <<<<<
  pos1 <- apply(REST, 1, function(x){
    return(x == 0)
  }) #spp absence (on current composition)
  pos1 <- as.data.frame(pos1)
  pos2 <- apply(addSimSel$compositions, 1, function(x){
    return(x > 0)
  }) #spp ocurrence after adding new individuals
  pos2 <- as.data.frame(pos2)
  pos3 <- mapply('&', pos1, pos2) #existant after adding individuals, but not initially
  all((pos1[,1] & pos2[,1]) == pos3[,1]) #checking
  all((pos1[,10] & pos2[,10]) == pos3[,10]) #checking
  rich_new_spp <- colSums(pos3)
  addSimSel$parameters <- cbind(addSimSel$parameters,
                                richness_new_spp = rich_new_spp)
  
  #species richness per set: <<<<<<<<<<<<<<<<<<<<
  spp_list <- apply(addSimSel$compositions, 1, function(x){
    x <- round(x,5)
    spp_list_x <- names(x)[x>0]
    return(spp_list_x)
  })
  sets_rich <- sapply(spp_list, function(x){
    rest <- intersect(x, spp_rest)
    ref <- intersect(x, spp_ref)
    arrab <- intersect(x, spp_arr)
    plan <- intersect(x, spp_plan)
    spon <- intersect(x, spp_spon)
    out <- c(total = length(x),
             ref = length(ref),
             rest = length(rest),
             arrab = length(arrab),
             plan = length(plan),
             spon = length(spon))
    return(out)
  })
  sets_rich <- t(sets_rich)
  
  # >>> export: #################################
  FILENAME <- paste0('outputs/',SIMULATION_INDEX,'/','addSimSel_multifunctionality.csv')
  write.csv2(addSimSel$multifunctionality, FILENAME)
  
  FILENAME <- paste0('outputs/',SIMULATION_INDEX,'/','addSimSel_parameters.csv')
  write.csv2(addSimSel$parameters, FILENAME, row.names = TRUE)
  
  FILENAME <- paste0('outputs/',SIMULATION_INDEX,'/','addSimSel_composition.csv')
  colnames(addSimSel$compositions)[6] <- 'Aristolochia baetica'
  colnames(addSimSel$compositions)[20] <- 'Bupleurum fruticosum'
  colnames(addSimSel$compositions)[130] <- 'Thymus mastichina'
  write.csv2(addSimSel$compositions, FILENAME, row.names = TRUE)
  
  FILENAME <- paste0('outputs/',SIMULATION_INDEX,'/','addSimSel_sets_richness.csv')
  write.csv2(sets_rich, FILENAME, row.names = TRUE)
  
  run_up_to <- c(run_up_to, SIMULATION_INDEX)
  write.table(run_up_to, 'run_up_to.txt')
  
  rm(list = c('addSim', 'addSim_r', 'addSim_p'))
  FIM <- Sys.time()
  print('SIMULATION:')
  print(SIMULATION_INDEX)
  print(DIFERENCA <- FIM - INI)
  print('+++++++++++++++++++++++++++++++++++++++')
}

# >> New restoration area #######################

# Simulation: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
SR<- c(10, 26) #alfa richness in reference sites
#N <- 10000
N <- 100
cwm	<- c('LMA', 'Dur_flowering',
         'Entomophily', 'Zooc', 'Resprouter')
cwv <-  'Height'
ava <- data$Planted

pos <- und == 0
traits_new <- traits[pos,] #elimina Pinus
ava_new <- ava[pos]
REF_new <- REF[,pos]

ini <- Sys.time()
comSim <- comSimulation(trait = traits_new, ava = ava_new,
                        it = 20000, rich = SR,
                        cwm = cwm, cwv = cwv,
                        ref = REF_new)
fim <- Sys.time()
(dif <- fim - ini)
FILENAME <- paste0('outputs/new_area/comSim.rds')
saveRDS(comSim, FILENAME)
comSim <- readRDS(FILENAME)

# Calculate multifunctionality: <<<<<<<<<<<<<<<<<

# Define thresholds: <<<<<<<<<<<<<<<<<<<<<<<<<<<<
ref_fun <- comSim$ref_communities$parameters[,c(-1,-2)] #reference function columns

thrs <- apply(ref_fun, 2, mean)

mult_n <- multCalculation(x = list(new_areas = comSim$sim_communities$parameters),
                          th = thrs)
comSim$sim_communities$multifunctionality <- mult_n$new_areas

# Simulation selection: <<<<<<<<<<<<<<<<<<<<<<<<<
comSim_m <- comSim$sim_communities$multifunctionality
pos <- comSim_m[,"multifunctionality"] == 6
comSimSel <- comSim
comSimSel$sim_communities$parameters <- comSimSel$sim_communities$parameters[pos,]
comSimSel$sim_communities$composition <- comSimSel$sim_communities$composition[pos,]
comSimSel$sim_communities$multifunctionality <- comSimSel$sim_communities$multifunctionality[pos,]
saveRDS(comSimSel, "outputs/new_area/comSimSel.rds")

# Species richness per set:
spp_list <- apply(comSim$sim_communities$composition[pos,], 1, function(x){
  x <- round(x,5)
  spp_list_x <- names(x)[x>0]
  return(spp_list_x)
})
saveRDS(spp_list, 'outputs/spp_list_newareas.rds')
sets_rich_new <- sapply(spp_list, function(x){
  total <- length(x)
  rest <- length(intersect(x, spp_rest) )
  ref <- length(intersect(x, spp_ref) )
  arrab <- length(intersect(x, spp_arr) )
  plan <- length(intersect(x, spp_plan) )
  not_plan <- total - plan
  out <- c(total = total,
           ref = ref,
           rest = rest,
           arrab = arrab,
           plan = plan,
           not_plan = not_plan
  )
  return(out)
})
sets_rich_new <- t(sets_rich_new)
write.table(sets_rich_new, 'outputs/new_area/sets_rich_new.txt')












# END