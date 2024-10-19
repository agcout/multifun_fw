############################################
# BIND ALL MATRICES ########################
############################################

# Multifunctionality #######################
files <- list.files('outputs')
pos <- files %in% as.character(1:20)
files <- files[pos]
files <- files[order(nchar(files), files)] #https://stackoverflow.com/questions/17531403/how-to-sort-a-character-vector-where-elements-contain-letters-and-numbers
files <- paste0('outputs/', files)
tab_list <- lapply(files, FUN = function(x){
  fn <- paste0(x, '/addSimSel_multifunctionality.csv' )
  tab_x <- read.csv2(fn)
  print(dim(tab_x))
  return(tab_x)
})
tab_all <- do.call(rbind, tab_list)
colnames(tab_all)[1] <- 'Site'
tab_all$Density <- rep(seq(0.1,2,0.1), each = 59)
View(tab_all)

write.table(tab_all, 'outputs/binded_dfs/addSimSel_multifunctionality.txt', row.names = FALSE)

# Parameters ###############################
files <- list.files('outputs')
pos <- files %in% as.character(1:20)
files <- files[pos]
files <- files[order(nchar(files), files)] #https://stackoverflow.com/questions/17531403/how-to-sort-a-character-vector-where-elements-contain-letters-and-numbers
files <- paste0('outputs/', files)
tab_list <- lapply(files, FUN = function(x){
  fn <- paste0(x, '/addSimSel_parameters.csv' )
  tab_x <- read.csv2(fn)
  print(dim(tab_x))
  return(tab_x)
})
tab_all <- do.call(rbind, tab_list)
colnames(tab_all)[1] <- 'Site'
tab_all$Density <- rep(seq(0.1,2,0.1), each = 59)
View(tab_all)

write.table(tab_all, 'outputs/binded_dfs/addSimSel_parameters.txt', row.names = FALSE)

# Composition ###################################
files <- list.files('outputs')
pos <- files %in% as.character(1:20)
files <- files[pos]
files <- files[order(nchar(files), files)] #https://stackoverflow.com/questions/17531403/how-to-sort-a-character-vector-where-elements-contain-letters-and-numbers
files <- paste0('outputs/', files)
tab_list <- lapply(files, FUN = function(x){
  fn <- paste0(x, '/addSimSel_composition.csv' )
  tab_x <- read.csv2(fn)
  print(dim(tab_x))
  return(tab_x)
})
tab_all <- do.call(rbind, tab_list)
colnames(tab_all)[1] <- 'Site'
tab_all <- data.frame(Density = rep(seq(0.1,2,0.1), each = 59),
                      tab_all)
View(tab_all)

write.table(tab_all, 'outputs/binded_dfs/addSimSel_composition.txt', row.names = FALSE)

# Composition new species #######################
addSim <- readRDS('outputs/1/addSim.rds')
REST <- addSim$compositions$restored #composition of restored sites
#REST <- round(REST, 3)
REST[REST>0] <- 1
files <- list.files('outputs')
pos <- files %in% as.character(1:20)
files <- files[pos]
files <- files[order(nchar(files), files)] #https://stackoverflow.com/questions/17531403/how-to-sort-a-character-vector-where-elements-contain-letters-and-numbers
files <- paste0('outputs/', files)
tab_list <- lapply(files, FUN = function(x){
  fn <- paste0(x, '/addSimSel_composition.csv' ) #complete path
  tab_x <- read.csv2(fn, check.names = FALSE) #composition of simulation
  row.names(tab_x) <- tab_x[,1] #add row names
  #tab_x <- round(tab_x[,-1],3) #delete first column (row names)
  tab_x <- tab_x[,-1]
  tab_x[tab_x > 0] <- 1 #transform to occurrence table
  tab_sum <- tab_x + REST #sum to REST
  tab_sum[tab_sum == 2] <- 0 #if it is equal to 2, it is not added
  return(tab_sum)
})
tab_all <- do.call(rbind, tab_list)
tab_all <- data.frame(Density = rep(seq(0.1,2,0.1), each = 59),
                      Site = row.names(tab_list[[1]]),
                      tab_all,
                      check.names = FALSE
                      )
row.names(tab_all) <- 1:nrow(tab_all)

write.table(tab_all, 'outputs/binded_dfs/addSimSel_composition_new_spp.txt', row.names = FALSE)











# end