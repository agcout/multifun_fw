
# FIGURES #######################################
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# > DATA ##########################################
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(dplyr)

dados <- read.table('outputs/binded_dfs/addSimSel_multifunctionality.txt', h=TRUE)
colnames(dados)[8] <- 'Multifunctionality'
dados$Density <- dados$Density * 100

dados2 <- read.table('outputs/binded_dfs/addSimSel_parameters.txt', h=TRUE)

dados3 <- read.table('outputs/binded_dfs/addSimSel_composition.txt', h=TRUE)

dados4 <- read.table('outputs/binded_dfs/addSimSel_composition_new_spp.txt', h=TRUE, check.names = FALSE)

spp_sets <- readRDS('outputs/species_sets/spp_sets.rds')

# > Multifunctionality vs. density and richness #######
# >> Vs. density ############################
# Total LM:
mt_gama <- aggregate(Multifunctionality ~ Density,
                     data = dados, FUN = sum)
max_mult <- 6*59 #maximum total of multifunctionality: 6 function x 59 sites
mt_gama$Multifunctionality <- mt_gama$Multifunctionality/max_mult*100

# Total resprouters:
mt_gama_resp <- aggregate(Resprouter ~ Density,
                          data = dados, FUN = sum)
mt_gama_resp$Resprouter <- mt_gama_resp$Resprouter/max_mult*100

# Border sites, resprouters:
bur_sit <- read.csv2('burned_sites.csv') #burned sites table
bur_sit_names <- bur_sit$site[bur_sit$fire == 1] #burned sites
pos <- dados$Site %in% bur_sit_names
dados_bur <- dados[pos,] #data subset with burned sites only
mt_gama_bur <- aggregate(Resprouter ~ Density,
                         data = dados_bur, FUN = sum)
mt_gama_bur$Resprouter <- mt_gama_bur$Resprouter/8*100

# bind data frames:
mt_gama_all <- cbind(mt_gama,
                     Resprouter = mt_gama_resp$Resprouter,
                     Resprouter_border = mt_gama_bur$Resprouter)

# long format:
mt_gama_lon <- pivot_longer(mt_gama_all,
                            cols = c('Multifunctionality',
                                     'Resprouter',
                                     'Resprouter_border'),
                            names_to = 'MF_type'
                            )
colnames(mt_gama_lon)[3] <- 'MF'
pos <- mt_gama_lon$MF_type == 'Multifunctionality'
mt_gama_lon$MF_type[pos] <- 'LM (all sites)'
pos <- mt_gama_lon$MF_type == 'Resprouter'
mt_gama_lon$MF_type[pos] <- 'Fire resistance (all sites)'
pos <- mt_gama_lon$MF_type == 'Resprouter_border'
mt_gama_lon$MF_type[pos] <- 'Fire resistance (border sites)'

# plot:
cols <- brewer.pal(3,'Dark2')
reference_mf <- (c(6+10+10+4+8+11)/114)*100
current_mf <- (c(33+30+7+40+34)/max_mult)*100
planted_mf <- (c(36+23+1)/max_mult)*100
p1 <- ggplot(mt_gama_lon, aes(x=Density, y=MF, color = MF_type)) +
  geom_line() + geom_point() +
  scale_color_manual(values = cols) +
  ylab("Landscape multifunctionality (LM)") + xlab("Individuals added (%)")+
  theme_bw() +
  geom_abline(slope=0, intercept = reference_mf, lty=2) +
  annotate("text", x=175+6, y=reference_mf+3, label= "Reference") +
  geom_abline(slope=0, intercept = current_mf, lty=2) +
  annotate("text", x=175+10, y=current_mf-3, label= "Current") +
  geom_abline(slope=0, intercept = planted_mf, lty=2) +
  annotate("text", x=175+3, y=planted_mf+3, label= "Planted only") +
  theme(legend.title = element_blank(),
        legend.position=c(.22*1.3,.865),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
        #panel.border = element_blank()
        )
p1
ggsave('outputs/figures/LM_density.png', p1,
       scaling = 0.6,
       width = 4, height = 2.5,
       bg = 'white')

# >> Vs. richness ##########################
#Average:
rich_av <- aggregate(richness_new_spp ~ Density,
                     data = dados2,
                     FUN = mean)
rich_av$MF <- mt_gama$Multifunctionality

#Total:
#needs to be calculated
dados4_s <- split(dados4, dados4$Density) #dados4 split
rich_tot <- sapply(dados4_s, FUN = function(x){
  tot1 <- colSums(x[,-c(1,2)])
  added_alfa <- names(tot1[tot1>0]) #added species in alfa level
  added_gama <- added_alfa[!added_alfa%in%spp_sets$rest] #added species in gama level
  tot2 <- length(added_gama)
  return(tot2)
})
names(rich_tot)

#Figure:
mt_gama$Average <- rich_av$richness_new_spp
mt_gama$Total <- rich_tot
mt_gama_lon <- mt_gama %>% pivot_longer(cols=c('Average', 'Total'),
                                        names_to='Type',
                                        values_to = 'Richness')
mt_gama_lon$Type[mt_gama_lon$Type=='Average'] <- 'Average richness'
mt_gama_lon$Type[mt_gama_lon$Type=='Total'] <- 'Total richness'
cols <- brewer.pal(5,'Dark2')[4:5]
p2 <- ggplot(mt_gama_lon, aes(x=Richness, y=Multifunctionality, color = Type)) +
  #geom_line() +
  geom_path() +
  geom_point() + theme_bw() + ylim(0,100) +
  scale_color_manual(values = cols) + 
  xlab('Species added') +
  theme(legend.title = element_blank(),
        legend.position=c(.8,.1),
        legend.background = element_blank(),
        #legend.box.background = element_rect(colour = "white"),
        legend.box.background = element_blank(),
        axis.title.y = element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
        #panel.border = element_blank()
  )
p2

ggsave('outputs/figures/LM_richness.png', p2,
       scaling = 0.6,
       width = 4, height = 2.5,
       bg = 'white')

# >> Nested figure #####
fator <- 0.9 #correction factor for figures
library(gridExtra)
nested <- grid.arrange(p1,p2, ncol=2)
ggsave('outputs/figures/LM_nested.png', nested,
       width = 10.5*fator, height = 5*fator)

# > Multifunctionality with max_add = 12 ####################################
# Figure showing multifunctionality when max_add = 1.2

library(ComplexUpset)
library(RColorBrewer)
library(ggplot2)
cols <- brewer.pal(6,'Pastel1')[c(1,2,3,4,5,6)]
bkg <- brewer.pal(9,'Greys')[3] #background
resol <- 400
cols_bar <- viridis::viridis(7, direction = -1)

# >> planted and spontaneous: #################
rest_mult <- read.table('outputs/12/addSim_multifunctionality_restored.txt')
rest_mult <- rest_mult[,1:6]
colnames(rest_mult)
colnames(rest_mult) <- c('DR',
                         'RPd',
                         'RPa',
                         'RF',
                         'FR',
                         'VS')
ini <- 9.1; sep <- 0.5; fim <- ini+(6*sep)
xleg <- seq(ini,fim,sep)
yleg <- rep(32,7)
legenda <- data.frame(x = xleg, y = yleg )

setsize <- colSums(rest_mult)
mult_gama <- mean(setsize)/nrow(rest_mult)*100
seg_data <- data.frame(x1=0.5, x2=6.5,
                       y1=mult_gama, y2=mult_gama)

up_current <- upset(rest_mult, 
                    c('FR', 'DR', 'VS', 'RF', 'RPa','RPd'),
                    name = 'Function combination',
                    base_annotations = list(
                      'Intersection size'=
                        intersection_size() 
                      + annotate(geom='text', x=Inf, y=Inf,
                                 label=paste('Total:', nrow(rest_mult)),
                                 vjust=1, hjust=1, size = 3)
                      #LEGEND:
                      + annotate(geom='text', x=xleg[c(1,7)], y=yleg[1]-3,
                                 label=c('0', '6'),size = 3)
                      + annotate(geom='text', x=xleg[4], y=yleg[1]+3,
                                 label='Number of functions:', size = 3)
                      + geom_point(data = legenda, aes(x=x,y=y),
                                   shape = 15, size=4, color=cols_bar[7:1])
                      + ylab('Percentage of sites')
                      + theme(plot.background=element_rect(fill=bkg))
                      + scale_y_continuous(breaks = c(0,11.8,23.6,35.4),
                                           labels=scales::percent_format(
                                             scale=100 / nrow(rest_mult)))
                      + coord_cartesian(ylim = c(0,35))
                    ),
                    height_ratio = 1,
                    set_sizes=(
                      upset_set_size(
                        mapping=aes(y=after_stat(count)/max(after_stat(count))*0.6779661*100)
                      )
                      + geom_text(aes(label=after_stat(count)), hjust=1.1, stat='count')
                      + expand_limits(y=120)
                      + ylab('Percentage of sites')
                      + geom_segment(aes(x=x1, xend=x2, y=y1, yend=y2),
                                     data=seg_data, linetype = "dashed",
                                     size=1.3)
                    ),
                    stripes = cols,
                    sort_intersections_by=c('degree','cardinality'),
                    sort_sets = FALSE,
                    
                    queries=list( #<<<<<<<<<<<<<<<<
                      upset_query(
                        intersect=c('DR', 'VS', 'RF', 'RPa','RPd'),
                        fill=cols_bar[2],
                        color=cols_bar[2],
                        only_components=c('Intersection size')
                      ),

                      upset_query(
                        intersect=c('DR', 'VS', 'RPa','RPd'),
                        fill=cols_bar[3],
                        color=cols_bar[3],
                        only_components=c('Intersection size')
                      ),
                      upset_query(
                        intersect=c('VS', 'RF', 'RPa','RPd'),
                        fill=cols_bar[3],
                        color=cols_bar[3],
                        only_components=c('Intersection size')
                      ),
                      upset_query(
                        intersect=c('DR', 'RF', 'RPa','RPd'),
                        fill=cols_bar[3],
                        color=cols_bar[3],
                        only_components=c('Intersection size')
                      ),

                      upset_query(
                        intersect=c('DR', 'VS', 'RPd'),
                        fill=cols_bar[4],
                        color=cols_bar[4],
                        only_components=c('Intersection size')
                      ),
                      upset_query(
                        intersect=c('VS', 'RPa','RPd'),
                        fill=cols_bar[4],
                        color=cols_bar[4],
                        only_components=c('Intersection size')
                      ),
                      upset_query(
                        intersect=c('RF', 'RPa','RPd'),
                        fill=cols_bar[4],
                        color=cols_bar[4],
                        only_components=c('Intersection size')
                      ),
                      upset_query(
                        intersect=c('DR', 'VS', 'RPa'),
                        fill=cols_bar[4],
                        color=cols_bar[4],
                        only_components=c('Intersection size')
                      ),
                      upset_query(
                        intersect=c('DR', 'RPa','RPd'),
                        fill=cols_bar[4],
                        color=cols_bar[4],
                        only_components=c('Intersection size')
                      ),

                      upset_query(
                        intersect=c('DR', 'VS'),
                        fill=cols_bar[5],
                        color=cols_bar[5],
                        only_components=c('Intersection size')
                      ),
                      upset_query(
                        intersect=c('RPa','RPd'),
                        fill=cols_bar[5],
                        color=cols_bar[5],
                        only_components=c('Intersection size')
                      ),
                      upset_query(
                        intersect=c('VS', 'RPa'),
                        fill=cols_bar[5],
                        color=cols_bar[5],
                        only_components=c('Intersection size')
                      ),
                      upset_query(
                        intersect=c( 'VS', 'RPd'),
                        fill=cols_bar[5],
                        color=cols_bar[5],
                        only_components=c('Intersection size')
                      ),
                      upset_query(
                        intersect=c('RF','RPd'),
                        fill=cols_bar[5],
                        color=cols_bar[5],
                        only_components=c('Intersection size')
                      ),
                      upset_query(
                        intersect=c('DR', 'RPd'),
                        fill=cols_bar[5],
                        color=cols_bar[5],
                        only_components=c('Intersection size')
                      ),
                      upset_query(
                        intersect=c('DR', 'RPa'),
                        fill=cols_bar[5],
                        color=cols_bar[5],
                        only_components=c('Intersection size')
                      ),

                      upset_query(
                        intersect=c('VS'),
                        fill=cols_bar[6],
                        color=cols_bar[6],
                        only_components=c('Intersection size')
                      ),
                      upset_query(
                        intersect=c('RPd'),
                        fill=cols_bar[6],
                        color=cols_bar[6],
                        only_components=c('Intersection size')
                      ),
                      upset_query(
                        intersect=c('RPa'),
                        fill=cols_bar[6],
                        color=cols_bar[6],
                        only_components=c('Intersection size')
                      ),

                      upset_query(
                        intersect=NA,
                        fill=cols_bar[7],
                        color=cols_bar[7],
                        only_components=c('Intersection size')
                      )


                    ) #<<<<<<<<<<<<<<<<
)
up_current

resol <- 400
filename <- paste0('outputs/figures/sets_restored_',
                   resol,'dpi.png')
ggsave(filename, up_current,
       dpi = resol, scale = 2,
       width = 4, height = 2.5)

# >> planted: #################################
plan_mult <- read.table('outputs/12/addSim_p_multifunctionality_restored.txt')
plan_mult <- plan_mult[,1:6]
colnames(plan_mult)
colnames(plan_mult) <- c('DR',
                         'RPd',
                         'RPa',
                         'RF',
                         'FR',
                         'VS')
ini <- 2.1; sep <- 0.11; fim <- ini+(6*sep)
xleg <- seq(ini,fim,sep)
yleg <- rep(32,7)
legenda <- data.frame(x = xleg, y = yleg )

setsize <- colSums(plan_mult)
mult_gama <- mean(setsize)/nrow(plan_mult)*100
seg_data <- data.frame(x1=0.5, x2=6.5,
                       y1=mult_gama, y2=mult_gama)

up_planted <- upset(plan_mult,
                    c('FR', 'DR', 'VS', 'RF', 'RPa','RPd'),
                    name = 'Function combination',
                    base_annotations=list(
                      'Intersection size'=intersection_size() 
                      + annotate(geom='text', x=Inf, y=Inf,
                                 label=paste('Total:', nrow(plan_mult)),
                                 vjust=1, hjust=1, size = 3)
                      + annotate(geom='text', x=xleg[c(1,7)], y=yleg[1]-3,
                                 label=c('0', '6'),size = 3)
                      + annotate(geom='text', x=xleg[4], y=yleg[1]+3,
                                 label='Number of functions:', size = 3)
                      + geom_point(data = legenda, aes(x=x,y=y),
                                   shape = 15, size=4, color=cols_bar[7:1])
                      
                      + ylab('Percentage of sites')
                      + theme(plot.background=element_rect(fill=bkg) )
                      + scale_y_continuous(breaks = c(0,11.8,23.6,35.4),
                                           labels=scales::percent_format(
                                             scale=100 / nrow(plan_mult)))
                      + coord_cartesian(ylim = c(0,35))
                    ),
                    height_ratio = 1,
                    set_sizes=(
                      upset_set_size(
                        mapping=aes(y=after_stat(count)/max(after_stat(count))*0.6101695*100)
                        #mapping=aes(y=..count../max(..count..))
                      )
                      + geom_text(aes(label=after_stat(count)), hjust=1.1, stat='count')
                      + expand_limits(y=120)
                      + ylab('Percentage of sites')
                      + geom_segment(aes(x=x1, xend=x2, y=y1, yend=y2),
                                     data=seg_data, linetype = "dashed",
                                     size=1.3)
                    ),
                    stripes=cols,
                    sort_intersections_by=c('degree','cardinality'),
                    sort_sets = FALSE,
                    
                    queries=list( #<<<<<<<<<<<<<<<<
                      upset_query(
                        intersect=c('DR', 'VS'),
                        fill=cols_bar[5],
                        color=cols_bar[5],
                        only_components=c('Intersection size')
                      ),
                      upset_query(
                        intersect=c('VS'),
                        fill=cols_bar[6],
                        color=cols_bar[6],
                        only_components=c('Intersection size')
                      ),
                      upset_query(
                        intersect=c('RPd'),
                        fill=cols_bar[6],
                        color=cols_bar[6],
                        only_components=c('Intersection size')
                      ),
                      upset_query(
                        intersect=NA,
                        fill=cols_bar[7],
                        color=cols_bar[7],
                        only_components=c('Intersection size')
                      )
                    ) #<<<<<<<<<<<<<<<<
)
up_planted
resol <- 400
filename <- paste0('outputs/figures/sets_planted_',
                   resol,'dpi.png')
ggsave(filename, up_planted,
       dpi = resol, scale = 2,
       width = 4, height = 2.5)

# >> planted, spontaneous, added: #############
pos <- dados$Density == 120
add_mult <- dados[pos,2:7]
colnames(add_mult)
colnames(add_mult) <- c('DR',
                        'RPd',
                        'RPa',
                        'RF',
                        'FR',
                        'VS')
ini <- 3.1; sep <- 0.18; fim <- ini+(6*sep)
xleg <- seq(ini,fim,sep)
yleg <- rep(32,7)
legenda <- data.frame(x = xleg, y = yleg )

setsize <- colSums(add_mult)
mult_gama <- mean(setsize)/nrow(add_mult)*100
seg_data <- data.frame(x1=0.5, x2=6.5,
                       y1=mult_gama, y2=mult_gama)

up_added <- upset(add_mult, 
                  c('FR', 'DR', 'VS', 'RF', 'RPa','RPd'),
                  name = 'Function combination',
                  base_annotations=list(
                    'Intersection size'=intersection_size() 
                    + annotate(geom='text', x=Inf, y=Inf,
                               label=paste('Total:', nrow(add_mult)),
                               vjust=1, hjust=1, size = 3)
                    + ylab('Percentage of sites')
                    + theme(plot.background=element_rect(fill=bkg))
                    + scale_y_continuous(breaks = c(0,11.8,23.6,35.4),
                                         labels=scales::percent_format(
                                           scale=100 / nrow(add_mult)))
                    + coord_cartesian(ylim = c(0,35))
                    #LEGENDA:
                    + annotate(geom='text', x=xleg[c(1,7)], y=yleg[1]-3,
                               label=c('0', '6'),size = 3)
                    + annotate(geom='text', x=xleg[4], y=yleg[1]+3,
                               label='Number of functions:', size = 3)
                    + geom_point(data = legenda, aes(x=x,y=y),
                                 shape = 15, size=4, color=cols_bar[7:1])
                  ),
                  height_ratio = 1,
                  set_sizes=(
                    upset_set_size(
                      mapping=aes(y=after_stat(count)/max(after_stat(count))*0.9491525*100)
                    )
                    + geom_text(aes(label=after_stat(count)), hjust=1.1, stat='count')
                    + expand_limits(y=120)
                    + ylab('Percentage of sites')
                    + geom_segment(aes(x=x1, xend=x2, y=y1, yend=y2),
                                   data=seg_data, linetype = "dashed",
                                   size=1.3)
                    + geom_segment(aes(x=x1, xend=x2, y=y1, yend=y2),
                                   data=seg_data, linetype = "dashed",
                                   size=1.3)
                  ),
                  stripes=cols,
                  sort_intersections_by=c('degree','cardinality'),
                  sort_sets = FALSE,
                  
                  queries=list( #<<<<<<<<<<<<<<<<
                    upset_query(
                      intersect=c('FR', 'DR', 'VS', 'RF', 'RPa','RPd'),
                      fill=cols_bar[1],
                      color=cols_bar[1],
                      only_components=c('Intersection size')
                    ), #6 functions
                    
                    upset_query(
                      intersect=c('DR', 'VS', 'RF', 'RPa','RPd'),
                      fill=cols_bar[2],
                      color=cols_bar[2],
                      only_components=c('Intersection size')
                    ), #5 functions
                    upset_query(
                      intersect=c('FR', 'DR', 'VS', 'RPa','RPd'),
                      fill=cols_bar[2],
                      color=cols_bar[2],
                      only_components=c('Intersection size')
                    ), #5 functions
                    
                    upset_query(
                      intersect=c('DR', 'VS', 'RPa','RPd'),
                      fill=cols_bar[3],
                      color=cols_bar[3],
                      only_components=c('Intersection size')
                    ), #4 functions
                    
                    upset_query(
                      intersect=c('VS', 'RPa', 'RPd'),
                      fill=cols_bar[4],
                      color=cols_bar[4],
                      only_components=c('Intersection size')
                    ), #3 functions
                    upset_query(
                      intersect=c('DR', 'VS', 'RF'),
                      fill=cols_bar[4],
                      color=cols_bar[4],
                      only_components=c('Intersection size')
                    ), #3 functions
                    

                    upset_query(
                      intersect=c('VS', 'RPd'),
                      fill=cols_bar[5],
                      color=cols_bar[5],
                      only_components=c('Intersection size')
                    ) #2 functions

                  ) #<<<<<<<<<<<<<<<
)
up_added
resol <- 400
filename <- paste0('outputs/figures/sets_added_',
                   resol,'dpi.png')
ggsave(filename, up_added,
       dpi = resol, scale = 2,
       width = 4, height = 2.5)

# >> reference: ###############################
ref_mult <- read.table('outputs/5/addSim_multifunctionality_reference.txt')
ref_mult <- ref_mult[,1:6]
colnames(ref_mult)
colnames(ref_mult) <- c('DR',
                        'RPd',
                        'RPa',
                        'RF',
                        'FR',
                        'VS')
#para legenda:
ini <- 7.1; sep <- 0.3; fim <- ini+(6*sep)
xleg <- seq(ini,fim,sep)
yleg <- rep(3,7)
legenda <- data.frame(x = xleg, y = yleg )

setsize <- colSums(ref_mult)
mult_gama <- mean(setsize)/nrow(ref_mult)*100
seg_data <- data.frame(x1=0.5, x2=6.5,
                       y1=mult_gama, y2=mult_gama)

up_refer <- upset(ref_mult,
                  c('FR', 'DR', 'VS', 'RF', 'RPa','RPd'),
                  name = 'Function combination',
                  base_annotations=list(
                    'Intersection size'=intersection_size() 
                    + annotate(geom='text', x=Inf, y=Inf,
                               label=paste('Total:', nrow(ref_mult)),
                               vjust=1, hjust=1, size = 3)
                    + ylab('Percentage of sites')
                    + theme(plot.background=element_rect(fill=bkg))
                    + scale_y_continuous(breaks = c(0,0.95,1.9,2.85),
                                         labels=scales::percent_format(
                                           scale=100 / nrow(ref_mult)))
                    #LEGENDA:
                    + annotate(geom='text', x=xleg[c(1,7)], y=yleg[1]-0.25,
                               label=c('0', '6'),size = 3)
                    + annotate(geom='text', x=xleg[4], y=yleg[1]+0.25,
                               label='Number of functions:', size = 3)
                    + geom_point(data = legenda, aes(x=x,y=y),
                                 shape = 15, size=4, color=cols_bar[7:1])
                  ),
                  height_ratio = 1,
                  set_sizes=(
                    upset_set_size(
                      mapping=aes(y=after_stat(count)/max(after_stat(count))*0.5526316*100*1.05)
                    )
                    + geom_text(aes(label=after_stat(count)), hjust=1.1, stat='count')
                    + expand_limits(y=70)
                    + ylab('Percentage of sites')
                    + geom_segment(aes(x=x1, xend=x2, y=y1, yend=y2),
                                   data=seg_data, linetype = "dashed",
                                   size=1.3)
                  ),
                  stripes=cols,
                  sort_intersections_by=c('degree','cardinality'),
                  sort_sets = FALSE,
                  queries=list( #<<<<<<<<<<<<<
                    upset_query(
                      intersect=c('FR', 'DR', 'RF', 'RPa'),
                      fill=cols_bar[3],
                      color=cols_bar[3],
                      only_components=c('Intersection size')
                    ),
                    upset_query(
                      intersect=c('FR', 'DR', 'VS', 'RPa'),
                      fill=cols_bar[3],
                      color=cols_bar[3],
                      only_components=c('Intersection size')
                    ),
                    upset_query(
                      intersect=c('FR', 'DR', 'VS', 'RF'),
                      fill=cols_bar[3],
                      color=cols_bar[3],
                      only_components=c('Intersection size')
                    ),
                    
                    upset_query(
                      intersect=c('RF', 'RPa','RPd'),
                      fill=cols_bar[4],
                      color=cols_bar[4],
                      only_components=c('Intersection size')
                    ),
                    upset_query(
                      intersect=c('VS', 'RPa','RPd'),
                      fill=cols_bar[4],
                      color=cols_bar[4],
                      only_components=c('Intersection size')
                    ),
                    upset_query(
                      intersect=c('FR', 'RPa','RPd'),
                      fill=cols_bar[4],
                      color=cols_bar[4],
                      only_components=c('Intersection size')
                    ),
                    upset_query(
                      intersect=c('FR', 'RF', 'RPa'),
                      fill=cols_bar[4],
                      color=cols_bar[4],
                      only_components=c('Intersection size')
                    ),
                    upset_query(
                      intersect=c('FR', 'DR', 'RF'),
                      fill=cols_bar[4],
                      color=cols_bar[4],
                      only_components=c('Intersection size')
                    ),
                    
                    upset_query(
                      intersect=c('FR', 'RF'),
                      fill=cols_bar[5],
                      color=cols_bar[5],
                      only_components=c('Intersection size')
                    ),
                    upset_query(
                      intersect=c('RPa','RPd'),
                      fill=cols_bar[5],
                      color=cols_bar[5],
                      only_components=c('Intersection size')
                    ),
                    upset_query(
                      intersect=c('FR','RPd'),
                      fill=cols_bar[5],
                      color=cols_bar[5],
                      only_components=c('Intersection size')
                    ),
                    
                    upset_query(
                      intersect=c('DR'),
                      fill=cols_bar[6],
                      color=cols_bar[6],
                      only_components=c('Intersection size')
                    ),
                    upset_query(
                      intersect=c('VS'),
                      fill=cols_bar[6],
                      color=cols_bar[6],
                      only_components=c('Intersection size')
                    ),
                    
                    upset_query(
                      intersect=NA,
                      fill=cols_bar[7],
                      color=cols_bar[7],
                      only_components=c('Intersection size')
                    )
                    
                    
                  ) #<<<<<<<<<<<<<
)
up_refer
resol <- 400
filename <- paste0('outputs/figures/sets_reference_',
                   resol,'dpi.png')
ggsave(filename, up_refer,
       dpi = resol, scale = 2,
       width = 4, height = 2.5)

# > Histograms ref trait distribution: ##########
addSim <- readRDS('outputs/1/addSim.rds')
thrs <- readRDS('thresholds.rds')
thrs[6] <- thrs[6] * 547556.9
ref_par <- addSim$parameters$reference
ref_par2 <- ref_par[,c(3:7,9)]
ref_par2[,3:5] <- ref_par2[,3:5]*100
thrs[3:5] <- thrs[3:5]*100
colnames(ref_par2)[c(2,4,6)] <- c('Flowering duration',
                                  'Zoochory',
                                  'Height')
labels <- parse(text=c('mg/mm^2', 'Months', 'Percentage',
                       'Percentage', 'Percentage', 'm^2'))

resol <- 800
png("outputs/figures/histograms.png", width = resol*5, height = resol*5,
    res = resol)
par(mfrow=c(3,2), mar=c(5,4,2,1))
for(i in 1:ncol(ref_par2)){
  print(i)
  main <- colnames(ref_par2)[i]
  hist(ref_par2[,i], xlab='', main = main)
  mtext(labels[i], 1, 2.5, cex=0.7)
  abline(v=thrs[i], col='red', lwd=2)
}
dev.off()












# end
