overall_start_time = Sys.time()

library(AnVIL)
library(tidyverse)
library(gridExtra)
library(ggrepel)
library(Seurat)
library(janitor)
library(hdf5r)
library(magick)
library(MASS); select = dplyr::select; summarize = dplyr::summarize

# statistical parameters
ci_size = .95
alpha = 1 - ci_size # 0.05
qnormval = qnorm(1-alpha/2) # ~1.96

# function
percent = function(proportion,digits=2) {
  return ( gsub(' ','',paste(formatC(proportion*100, digits=digits, format='fg'),"%",sep="") ) )
}

meta = read_tsv('data_in/meta.tsv', col_types=cols()) %>%
  mutate(region_long = case_when(reg=='cb' ~ 'cerebellum',
                                 reg %in% c('fc','ss') ~ 'cortex',
                                 reg=='th' ~ 'thalamus')) %>%
  mutate(aggr_dataset=paste0(species,'_',region_long,'_full'))

aggrs = c("mouse_cerebellum_full", 
          "mouse_cortex_full", 
          "mouse_thalamus_full", 
          "cyno_cerebellum_full",
          "cyno_cortex_full") # unique(meta$aggr_dataset)

####
# ANALYTICAL DATASET AND FIGURE 1
####


h5_mouse_cortex_full     = Read10X_h5(paste0('data_in/mouse_cortex_full/count/filtered_feature_bc_matrix.h5'),use.names=T)
h5_mouse_cerebellum_full = Read10X_h5(paste0('data_in/mouse_cerebellum_full/count/filtered_feature_bc_matrix.h5'),use.names=T)
h5_mouse_thalamus_full   = Read10X_h5(paste0('data_in/mouse_thalamus_full/count/filtered_feature_bc_matrix.h5'),use.names=T)
h5_cyno_cortex_full      = Read10X_h5(paste0('data_in/cyno_cortex_full/count/filtered_feature_bc_matrix.h5'),use.names=T)
h5_cyno_cerebellum_full  = Read10X_h5(paste0('data_in/cyno_cerebellum_full/count/filtered_feature_bc_matrix.h5'),use.names=T)

zeroes = tibble(aggr=aggrs, 
                rowcount=integer(length(aggrs)),
                colcount=integer(length(aggrs)),
                p_zero=numeric(length(aggrs)))

panels = {}
panel = 1
disp_panel = 1

for (aggr in aggrs) {
  
  this_start_time = Sys.time()
  
  cat(file=stderr(), paste0('Working on ',aggr,' ...'))
  
  ####
  # Analytical dataset
  ####
  
  species = gsub('_.*','',aggr)
  aggr_disp = gsub('cyno','cynomolgus',gsub('_',' ',gsub('_full','',aggr)))
  prnp_rowname = case_when(species == 'cyno' ~ 'PRNP', species == 'mouse' ~ 'Prnp')
  malat1_rowname = case_when(species == 'cyno' ~ '', species == 'mouse' ~ 'Malat1')
  rnaseh1_rowname = case_when(species == 'cyno' ~ 'RNASEH1', species == 'mouse' ~ 'Rnaseh1')
  
  h5_path = paste0('data_in/',aggr,'/count/filtered_feature_bc_matrix.h5')
  cluster_path = paste0('data_in/',aggr,'/count/analysis/clustering/gene_expression_graphclust/clusters.csv')
  umap_path = paste0('data_in/',aggr,'/count/analysis/umap/gene_expression_2_components/projection.csv')
  aggregation_path = paste0('data_in/',aggr,'/aggregation.csv')
  start_time = Sys.time()
  cat(file=stderr(), '\nReading h5 file: ', h5_path)
  h5 = get(paste0('h5_',aggr)) # Read10X_h5(h5_path,use.names=T)
  clust = read_csv(cluster_path, col_types=cols()) %>% clean_names()
  aggregation = read_csv(aggregation_path, col_types=cols()) %>% clean_names() %>% mutate(h5_num=row_number())
  umap = read_csv(umap_path, col_types=cols()) %>% clean_names()
  
  amt = read_csv(paste0('data_in/cla_',aggr,'.csv'), col_types=cols()) %>% clean_names()
  mar = read_tsv(paste0('data_in/mar_',aggr,'.tsv'), col_types=cols())
  met = read_tsv(paste0('data_in/met_',aggr,'.tsv'), col_types=cols())
  amt$assignments[is.na(amt$assignments)] = 'exclude'
  amt$assignments[grepl('^exclude',amt$assignments)] = 'exclude'
  amt$amtfct = factor(amt$assignments, levels=c('exclude',rev(met$assignments)),  ordered=T)
  
  p_zero = sum(h5==0)/prod(dim(h5))
  zeroes$p_zero[zeroes$aggr==aggr] = NA #p_zero
  zeroes$rowcount[zeroes$aggr==aggr] = dim(h5)[1]
  zeroes$colcount[zeroes$aggr==aggr] = dim(h5)[2]
  
  # export analytical dataset
  tibble(barcode = colnames(h5),
         total_umis = colSums(h5),
         unique_genes = colSums(h5 > 0),
         malat1_umis = as.numeric(NA),
         prnp_umis = h5[prnp_rowname,],
         rnaseh1_umis = h5[rnaseh1_rowname,]) %>%
    inner_join(select(clust, barcode, cluster), by='barcode') %>%
    mutate(h5_num = as.numeric(substr(barcode, 18, 19))) %>%
    inner_join(select(aggregation, h5_num, sample_id), by='h5_num') %>%
    inner_join(meta, by=c('sample_id'='sample')) %>%
    rename(sample=sample_id) %>%
    inner_join(umap, by='barcode') -> analytical
  
  analytical %>% group_by(cluster) %>%
    summarize(.groups='keep', n_cells=n(),
              mean_umis = mean(total_umis), mean_genes = mean(unique_genes)) %>%
    ungroup() %>% View()
  
  # malat1_umis = case_when(malat1_rowname=='' ~ as.numeric(NA),
  #                          malat1_rowname!='' ~ h5[malat1_rowname,]),
  # the above throws an "invalid character indexing" error even though the h5[malat1_rowname,] should never be reached if blank. 
  # according to docs this is because https://dplyr.tidyverse.org/reference/case_when.html
  # "case_when() evaluates **all** RHS expressions, and then constructs its
  # result by extracting the selected (via the LHS expressions) parts."
  # instead, do this separately:
  if (malat1_rowname != '') analytical$malat1_umis = h5[malat1_rowname,]
  
  write_tsv(analytical, paste0('data_out/ana_',aggr,'.tsv'), na='')
  
  
  ####
  # Create UMAP and DotPlots for Figure 1
  ####
  
  seur = CreateSeuratObject(h5)
  Idents(seur) = factor(amt$amtfct)
  if ('exclude' %in% amt$assignments) {
    seur = subset(seur, idents='exclude', invert=T)
  }
  
  umap %>%
    inner_join(amt, by='barcode') %>%
    inner_join(met, by='assignments') -> umapr
  if (aggr=='mouse_cortex_full') {
    single_label = 'excitatory'
  } else {
    single_label = c('granule', 'excitatory','inhibitory')  
  }
  umapr %>%
    filter(abs(umap_1) > 3 | abs(umap_2) > 3) %>%
    filter(!(assignments %in% single_label)) %>%
    mutate(quadrant = paste0(umap_1 > 0, '-', umap_2 > 0)) %>%
    group_by(assignments, quadrant) %>%
    summarize(.groups='keep', 
              umap1_mean=mean(umap_1), 
              umap2_mean=mean(umap_2),
              n_cells = n()) %>%
    ungroup() %>%
    filter(n_cells > 10) %>%
    filter(!(assignments %in% single_label)) -> umap_centroids1
  umap_centroids1 %>%
    group_by(assignments) %>%
    summarize(.groups='keep', 
              n_quads = n(),
              dist=sqrt(max(diff(umap1_mean))^2 + max(diff(umap2_mean))^2)) -> quadcount
  umap_centroids1 %>%
    filter(assignments %in% quadcount$assignments[quadcount$n_quads > 1 & quadcount$dist > 5]) -> umap_centroids1
  if (aggr=='mouse_cortex_full') {
    umap_centroids1 %>%
      filter(!(assignments %in% 'astrocyte' & umap1_mean < 0 & umap2_mean > 0)) -> umap_centroids1
  }
  umapr %>%
    filter(!(assignments %in% quadcount$assignments[quadcount$n_quads > 1 & quadcount$dist > 5])) %>%
    group_by(assignments) %>%
    summarize(.groups='keep', umap1_mean=mean(umap_1), umap2_mean=mean(umap_2), n_cells=n()) %>%
    ungroup() %>% mutate(quadrant='') -> umap_centroids2
  rbind(umap_centroids1, umap_centroids2) -> umap_centroids
  
  this_umap_plot = ggplot(umapr, aes(umap_1, umap_2, colour=color)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    geom_point(size=0.5, shape=20) +
    scale_color_identity() +
    labs(x = 'UMAP 1', y='UMAP 2') +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank()) +
    geom_text_repel(data=umap_centroids, point.size = NA, size=3, min.segment.length = 0,
                    aes(umap1_mean, umap2_mean, label=assignments, color='black')) +
    ggtitle(aggr_disp) +
    labs(tag = LETTERS[disp_panel])
  panels[[panel]] = this_umap_plot
  panel = panel + 1
  disp_panel = disp_panel + 1
  
  
  ymeta = tibble(celltypelabel = levels(amt$amtfct)[!grepl('exclude',levels(amt$amtfct))]) %>%
    filter(celltypelabel %in% amt$assignments) %>%
    mutate(ymid = row_number()) %>%
    mutate(ybottom = ymid-0.5, ytop=ymid+0.5) %>%
    inner_join(met, by=c('celltypelabel'='assignments')) %>%
    mutate(xmin=0, xmax=1) %>%
    rename(typecolor=color)
  
  bottom_margin = case_when(species=='mouse' ~ 26, species=='cyno' ~ 30)

  dot_plot_y_blocks = ggplot(ymeta, aes(xmin=xmin, xmax=xmax, ymin=ybottom, ymax=ytop, colour=typecolor, label=celltypelabel)) + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank()) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "#FFFFFF00")) +
    labs(x='', y='') +
    scale_x_continuous(name="") + 
    scale_y_continuous(name="") +
    theme(plot.margin = margin(r = -11, l =-12, t=-6, b=bottom_margin)) +
    geom_rect(aes(colour=typecolor, fill=typecolor), alpha=1) +
    geom_text(aes(x=xmax,y=ymid,label=celltypelabel,colour='#000000',hjust=1),size=3) +
    scale_colour_identity(aesthetics=c('colour','fill')) 
  panels[[panel]] = dot_plot_y_blocks
  panel = panel + 1
  
  this_dot_plot = DotPlot(seur, features=mar$gene, dot.scale=4)  + 
    xlab('') + ylab('') + 
    theme(axis.text.y = element_blank(),
          axis.ticks = element_blank()) + # theme(axis.text.y = element_text(size=8)) +
    theme(plot.margin = margin(l = -12, t=1, b=1, r=0.1)) +
    theme(axis.text.x = element_text(angle = 45, hjust=1, size=7, face='italic')) +
    theme(legend.position="none")
  panels[[panel]] = this_dot_plot
  panel = panel + 1
  
  elapsed = Sys.time() - this_start_time
  
  cat(file=stderr(), paste0('completed in ',round(as.numeric(elapsed),1),' ',units(elapsed),'.\n'))
}

write_tsv(zeroes, 'data_out/zeroes.tsv')

resx=300
png(paste0('data_out/figure_1.png'),width=6.5*resx,height=12.0*resx,res=resx)
grid.arrange(grobs=panels, ncol=3, widths=c(1,.4,1.1))
dev.off()

overall_elapsed = (Sys.time() - overall_start_time)
cat(file=stderr(), paste0('done.\nAll tasks completed in ',round(as.numeric(overall_elapsed),1),' ',units(elapsed),'.'))


####
# EXCITATORY / INHIBITORY CLUSTERS
####

aggr = 'mouse_cortex_full'

ana = read_tsv(paste0('data_out/ana_',aggr,'.tsv'), col_types=cols())
amt = read_csv(paste0('data_in/cla_',aggr,'.csv'), col_types=cols()) %>% clean_names()
mar = read_tsv(paste0('data_in/mar_',aggr,'.tsv'), col_types=cols())
met = read_tsv(paste0('data_in/met_',aggr,'.tsv'), col_types=cols())
amt$assignments[is.na(amt$assignments)] = 'exclude'
amt$assignments[grepl('^exclude',amt$assignments)] = 'exclude'
amt$amtfct = factor(amt$assignments, levels=c('exclude',rev(met$assignments)),  ordered=T)
amt$subtype = as.character(amt$amtfct)

umap = read_csv(paste0('data_in/',aggr,'/count/analysis/umap/gene_expression_2_components/projection.csv'), col_types=cols()) %>% clean_names()
h5 = get(paste0('h5_',aggr))

ex6 = read_csv('data_in/ex6_recluster.csv', 
               col_types=cols()) %>% 
  clean_names() %>%
  filter(!(barcode %in% amt$barcode[amt$amtfct=='exclude']))
exc_mar = read_csv('data_in/mouse_excitatory_layer_markers.csv', 
                   col_types=cols()) %>% clean_names()

as.matrix(h5[exc_mar$name,ex6$barcode]) %>%
  t() %>%
  as_tibble() %>%
  add_column(barcode = ex6$barcode) %>%
  inner_join(ex6, by='barcode') -> exc_exp

exc_exp %>%
  mutate(cluster = as.numeric(gsub('Cluster ','',ex6_recluster))) %>%
  inner_join(select(ana, barcode, total_umis), by='barcode') %>%
  inner_join(select(amt, barcode, subtype), by='barcode') %>%
  filter(subtype=='excitatory') %>%
  pivot_longer(Slc17a7:Foxp2) %>%
  group_by(cluster, name) %>%
  summarize(.groups='keep',
            n_cells = length(unique(barcode)),
            rpm = 1e6*sum(value)/sum(total_umis)) %>%
  ungroup() %>%
  filter(n_cells > 500) -> exc_clust

write_tsv(exc_clust, 'data_out/mouse_exc_clust.tsv')

in6 = read_csv('data_in/in6_recluster.csv', 
               col_types=cols()) %>% 
  clean_names() %>%
  filter(!(barcode %in% amt$barcode[amt$amtfct=='exclude']))
inh_mar = read_csv('data_in/mouse_inhibitory_markers.csv', 
                   col_types=cols()) %>% clean_names()

as.matrix(h5[inh_mar$name,in6$barcode]) %>%
  t() %>%
  as_tibble() %>%
  add_column(barcode = in6$barcode) %>%
  inner_join(in6, by='barcode') -> inh_exp

inh_exp %>%
  mutate(cluster = as.numeric(gsub('Cluster ','',in6_recluster))) %>%
  inner_join(select(ana, barcode, total_umis), by='barcode') %>%
  inner_join(select(amt, barcode, assignments), by='barcode') %>%
  filter(assignments=='inhibitory') %>%
  pivot_longer(Vip:Calb2) %>%
  group_by(cluster, name) %>%
  summarize(.groups='keep',
            n_cells = length(unique(barcode)),
            rpm = 1e6*sum(value)/sum(total_umis)) %>%
  ungroup() %>%
  filter(n_cells > 100) -> inh_clust

write_tsv(inh_clust, 'data_out/mouse_inh_clust.tsv')






