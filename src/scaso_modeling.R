options(stringsAsFactors=F)
if(interactive()) setwd('~/d/sci/src/scaso')
suppressMessages(library(tidyverse))
suppressMessages(library(janitor))
suppressMessages(library(MASS)); select = dplyr::select; summarize = dplyr::summarize

# statistical parameters
ci_size = .95
alpha = 1 - ci_size # 0.05
qnormval = qnorm(1-alpha/2) # ~1.96

# metadata
tx_meta = read_tsv('data/meta/tx_meta.tsv', col_types=cols())
samples = read_tsv('data/meta/meta.tsv', col_types=cols()) %>%
  mutate(region_long = case_when(reg=='cb' ~ 'cerebellum',
                                 reg %in% c('fc','ss') ~ 'cortex',
                                 reg=='th' ~ 'thalamus')) %>%
  mutate(aggr_dataset=paste0(species,'_',region_long,'_full'))
aggrs = unique(samples$aggr_dataset)

for (aggr in aggrs) {
  ### by assignment section
  
  if (exists('aco')) rm('aco')
  if (exists('afi')) rm('afi')
  
  amt_path = paste0('data/analytical/cla_',aggr,'.csv')
  met_path = paste0('data/analytical/met_',aggr,'.tsv')
  ana_path = paste0('data/analytical/ana_',aggr,'.tsv')
  
  read_csv(amt_path, col_types=cols()) %>%
    mutate(assignments = replace_na(assignments,'exclude')) %>%
    filter(!grepl('^exclude',assignments)) %>%
    distinct(assignments) -> actual_subtypes
  
  met = read_tsv(met_path, col_types=cols()) %>%
    filter(assignments %in% actual_subtypes$assignments) %>%
    mutate(x=row_number()) %>%
    mutate(y = max(x) + 1 - x)
  
  amt = read_csv(amt_path, col_types=cols()) %>%
    clean_names() %>%
    mutate(assignments = replace_na(assignments,'exclude')) %>%
    filter(!grepl('^exclude',assignments)) %>%
    left_join(met, by='assignments') %>%
    rename(subtype = assignments) 
  
  ana = read_tsv(ana_path,col_types=cols())

  ana %>%
    inner_join(amt, by='barcode') %>%
    group_by(celltype, subtype, x, y) %>%
    summarize(.groups='keep',
              n_cells=n(),
              sum_total_umis = sum(total_umis),
              ctl_cell_ratio = mean(tx %in% c('pbs0','acsf')),
              ctl_read_ratio = sum(total_umis[tx %in% c('pbs0','acsf')])/sum(total_umis)) %>%
    ungroup() %>%
    rename(n_umis = sum_total_umis) -> subtypes
    
  subtypes %>%
    group_by(celltype) %>%
    summarize(.groups='keep',
              midx=mean(x), minx=min(x), maxx=max(x),
              midy=mean(y), miny=min(y), maxy=max(y)) -> celltypes
  
  ana %>%
    inner_join(amt, by='barcode') %>%
    group_by(h5_num, sample, tx, weeks, celltype, subtype, x, y) %>%
    summarize(.groups='keep',
              n_cells = n(),
              total_umis = sum(total_umis),
              prnp_umis = sum(prnp_umis),
              malat1_umis = sum(malat1_umis),
              rnaseh1_umis = sum(rnaseh1_umis)) %>%
    ungroup() %>%
    mutate(subtype_fct = as.factor(subtype)) %>%
    mutate(txc = case_when(tx == 'acsf' ~ '_ctl',
                           tx == 'pbs0' ~ '_ctl',
                           TRUE ~ tx)) %>%
    mutate(prnp_rpm = 1e6 * prnp_umis / total_umis) %>%
    mutate(total_umis_pseudocount = ifelse(total_umis==0, 1, total_umis)) -> by_amt
  
  samples %>%
    filter(aggr_dataset == aggr) %>%
    filter(!(tx %in% c('acsf','pbs0'))) %>%
    group_by(reg, tx, weeks) %>%
    summarize(.groups='keep',
              n_active=n()) %>%
    ungroup() %>%
    mutate(model = paste0(weeks,'-',tx)) -> models
  
  for (m in 1:nrow(models)) {
    
    subs = by_amt %>%
      filter(tx %in% c(models$tx[m], 'acsf', 'pbs0')) %>%
      filter(weeks == models$weeks[m])
    
    if (models$tx[m] %in% c('aso1','aso6','aso7')) {
      fit = glm.nb(prnp_umis ~ subtype_fct + subtype_fct:txc + offset(log(total_umis_pseudocount)), data=subs)
    } else if (models$tx[m] %in% c('asom')) {
      fit = glm.nb(malat1_umis ~ subtype_fct + subtype_fct:txc + offset(log(total_umis_pseudocount)), data=subs)
    }
    
    tx_colors = tx_meta %>% distinct(color, txc)
    
    # note that the join solely on subtype here creates a brittleness where subtype
    # has to be unique, so you can't have celltype/subtype combinations such as
    # "neuron/other" and "glia/other", you'll get a many-to-many join here
    
    subs %>%
      group_by(subtype, txc) %>%
      summarize(.groups='keep', 
                txc_cells = sum(n_cells)) %>%
      ungroup() %>%
      group_by(subtype) %>%
      mutate(total_cells = sum(txc_cells)) %>%
      ungroup() %>%
      select(txc, subtype, n_cells=txc_cells, total_cells) -> cell_counts
    
    as_tibble(summary(fit)$coefficients, rownames='coefname') %>%
      select(coefname, estimate=Estimate, se=`Std. Error`, z=`z value`, p=`Pr(>|z|)`) %>%
      mutate(subtype = as.factor(gsub('\\(Intercept\\)',levels(subs$subtype_fct)[1],gsub('subtype_fct','',gsub(':.*','',coefname))))) %>%
      mutate(txc = ifelse(grepl('txc',coefname),gsub('.*txc','',coefname),'_ctl')) %>%
      mutate(log_normed = ifelse(txc!='_ctl', estimate, 0)) %>%
      mutate(normed = exp(log_normed),
             l95 = exp(log_normed - qnormval*se),
             u95 = exp(log_normed + qnormval*se)) %>%
      inner_join(tx_colors, by='txc') %>%
      inner_join(subtypes[,c('celltype','subtype','x','y')], by=c('subtype')) %>%
      inner_join(cell_counts, by=c('txc','subtype')) -> ba_coefs
    
    subs$fit_resid = fit$residuals
    
    subs %>%
      select(-x, -y) %>%
      inner_join(select(ba_coefs, -n_cells, -total_cells), by=c('subtype_fct'='subtype', 'celltype', 'txc')) %>%
      mutate(point_estimate = exp(log_normed + fit_resid)) %>%
      select(-coefname, -estimate, -se, -z, -p, -log_normed, -normed, -l95, -u95) -> ba_fits
    
    ba_coefs$model = models$model[m]
    ba_fits$model = models$model[m]
    
    if(exists('aco')) {
      aco = rbind(aco, ba_coefs)
    } else {
      aco = ba_coefs
    }
    
    if (exists('afi')) {
      afi = rbind(afi, ba_fits)
    } else {
      afi = ba_fits
    }
    
  }
  
  write_tsv(aco, paste0('data/analytical/aco_',aggr,'.tsv'))
  write_tsv(afi, paste0('data/analytical/afi_',aggr,'.tsv'))
  
  rm(aco)
  rm(afi)
  
}

