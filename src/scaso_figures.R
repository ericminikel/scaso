##############
# STARTUP
##############

overall_start_time = Sys.time()
cat(file=stderr(), 'Loading dependencies...')

options(stringsAsFactors=F)
if(interactive()) setwd('~/d/sci/src/scaso')
suppressMessages(library(tidyverse))
suppressMessages(library(janitor))
suppressMessages(library(weights))
suppressMessages(library(openxlsx))
suppressMessages(library(Hmisc))
suppressMessages(library(MASS)); select = dplyr::select; summarize = dplyr::summarize


##############
# FUNCTIONS & CONSTANTS
##############

cat(file=stderr(), 'done.\nSetting constants and functions...')

ctl_color = '#A7A7A7'
mal_color = '#FF00FF'

percent = function(x, digits=0, signed=F) gsub(' ','',paste0(ifelse(x <= 0, '', ifelse(signed, '+', '')),formatC(100*x,format='f',digits=digits),'%'))

upper = function(x, ci=0.95) { 
  alpha = 1 - ci
  sds = qnorm(1-alpha/2)
  mean(x) + sds*sd(x)/sqrt(sum(!is.na(x)))
}
lower = function(x, ci=0.95) { 
  alpha = 1 - ci
  sds = qnorm(1-alpha/2)
  mean(x) - sds*sd(x)/sqrt(sum(!is.na(x)))
}

alpha = function(rgb_hexcolor, proportion) {
  hex_proportion = sprintf("%02x",round(proportion*255))
  rgba = paste(rgb_hexcolor,hex_proportion,sep='')
  return (rgba)
}

ci_alpha = 0.35 # degree of transparency for shading confidence intervals in plot

clipdist = function(x, minx, maxx) {
  return (pmin(maxx,pmax(minx,x)))
}

rbind_files = function(path, grepstring) {
  all_files = list.files(path, full.names=T)
  these_files = all_files[grepl(grepstring,all_files)]
  if (exists('rbound_table')) rm('rbound_table')
  for (this_file in these_files) {
    this_tbl = read_delim(this_file, col_types=cols()) %>% clean_names()
    this_tbl$file = gsub('.*\\/','',gsub('\\.[tc]sv','',this_file))
    if (exists('rbound_table')) {
      rbound_table = rbind(rbound_table, this_tbl)
    } else {
      rbound_table = this_tbl
    }
  }
  return (rbound_table)
}

meta = read_tsv('data/meta/meta.tsv', col_types=cols())
tx_meta = read_tsv('data/meta/tx_meta.tsv', col_types = cols())
region_meta = read_tsv('data/meta/region_meta.tsv', col_types = cols())

##############
# OUTPUT STREAMS
##############

cat(file=stderr(), 'done.\nCreating output streams...'); flush.console()

text_stats_path = 'display_items/stats_for_text.txt'
write(paste('Last updated: ',Sys.Date(),'\n',sep=''),text_stats_path,append=F) # start anew - but all subsequent writings will be append=T

supplement_path = 'display_items/supplement.xlsx'
supplement = createWorkbook()
# options("openxlsx.numFmt" = "0.00") # this looks better for residuals but terrible for p values and weeks post-dose
supplement_directory = tibble(name=character(0), title=character(0))
write_supp_table = function(tbl, title='') {
  # write Excel sheet for supplement
  table_number = length(names(supplement)) + 1
  table_name = paste0('s',formatC(table_number,'d',digits=0,width=2,flag='0'))
  addWorksheet(supplement,table_name)
  bold_style = createStyle(textDecoration = "Bold")
  writeData(supplement,table_name,tbl,headerStyle=bold_style,withFilter=T)
  freezePane(supplement,table_name,firstRow=T)
  saveWorkbook(supplement,supplement_path,overwrite = TRUE)
  
  # also write tab-sep version for GitHub repo
  write_tsv(tbl,paste0('display_items/table-',table_name,'.tsv'), na='')
  
  # and save the title in the directory tibble for later
  assign('supplement_directory',
         supplement_directory %>% add_row(name=table_name, title=title),
         envir = .GlobalEnv)
}


####
# SUPPLEMENTARY TABLES
####

# QC metrics

metrics = rbind_files('data/metrics/', 'csv')
metrics %>%
  inner_join(meta, by=c('file'='sample')) %>%
  inner_join(tx_meta, by='tx') %>%
  inner_join(region_meta, by='reg') %>%
  select(animal_short, species, treatment=disp, weeks, region, estimated_number_of_cells:median_umi_counts_per_cell) -> metrics_out

write_supp_table(metrics_out, 'Quality control metrics for single-nucleus sequencing data.')

write(paste("QC metrics | mean number of reads: ",formatC(mean(metrics$number_of_reads),format='d',big.mark=','),'\n',sep=''),text_stats_path,append=T)
write(paste("QC metrics | mean number of cells: ",formatC(mean(metrics$estimated_number_of_cells),format='d',big.mark=','),'\n',sep=''),text_stats_path,append=T)
write(paste("QC metrics | mean reads per cell: ",formatC(mean(metrics$mean_reads_per_cell),format='d',big.mark=','),'\n',sep=''),text_stats_path,append=T)
write(paste("QC metrics | mean of median UMI count per cell: ",formatC(mean(metrics$median_umi_counts_per_cell),format='d',big.mark=','),'\n',sep=''),text_stats_path,append=T)
write(paste("QC metrics | mean of median genes per cell: ",formatC(mean(metrics$median_genes_per_cell),format='d',big.mark=','),'\n',sep=''),text_stats_path,append=T)

# qPCR vs. SC measures of bulk target engagement

qpcr_all = read_tsv('data/other/qpcr_all.tsv', col_types=cols())

all_ana = rbind_files('data/analytical/','ana_.*full') %>%
  mutate(aggr = gsub('ana_','',file)) %>%
  select(-file) %>%
  relocate(aggr)

all_ana %>%
  group_by(aggr, tx, weeks, animal_short, reg) %>%
  summarize(.groups='keep',
            prnp_rpm = 1e6*mean(prnp_umis/total_umis),
            malat1_rpm = 1e6*mean(malat1_umis/total_umis)) %>%
  ungroup() %>%
  group_by(aggr, weeks, reg) %>%
  mutate(prnp_rel = prnp_rpm / mean(prnp_rpm[tx %in% c('acsf','pbs0')]),
         malat1_rel = malat1_rpm / mean(malat1_rpm[tx %in% c('acsf','pbs0')])) %>%
  ungroup() -> sc_indiv_all

qpcr_all %>%
  pivot_wider(names_from=analyte, values_from=mean_rel) %>%
  rename(malat1_qpcr = Malat1, prnp_qpcr = Prnp) -> qpcr_pivoted

sc_indiv_all %>%
  inner_join(tx_meta, by='tx') %>%
  inner_join(meta, by=c('animal_short','tx','weeks','reg')) %>%
  inner_join(region_meta, by=c('reg')) %>%
  left_join(qpcr_pivoted, by=c('animal_short'='animal','reg'='reg')) %>%
  select(aggr, species, treatment=disp, tx, weeks, reg, region, animal_short, 
         sc_prnp_rpm=prnp_rpm, sc_prnp_rel=prnp_rel, prnp_qpcr,
         sc_malat1_rpm=malat1_rpm, sc_malat1_rel=malat1_rel, malat1_qpcr, color) -> sc_vs_qpcr_indiv
  
sc_vs_qpcr_indiv %>%
  select(-tx, -reg, -color) -> sc_vs_qpcr_indiv_out

write_supp_table(sc_vs_qpcr_indiv_out, 'Bulk tissue target engagement quantified by single cell sequencing and by qPCR, all samples all conditions.')

sc_vs_qpcr_indiv %>%
  group_by(aggr, weeks, tx, treatment, reg, region, color) %>%
  summarize(.groups='keep', 
            prnp_qpcr_mean=mean(prnp_qpcr), 
            prnp_qpcr_l95=lower(prnp_qpcr), 
            prnp_qpcr_u95=upper(prnp_qpcr), 
            prnp_sc_mean=mean(sc_prnp_rel), 
            prnp_sc_l95=lower(sc_prnp_rel), 
            prnp_sc_u95=upper(sc_prnp_rel), 
            malat1_qpcr_mean= mean(malat1_qpcr), 
            malat1_qpcr_l95= lower(malat1_qpcr), 
            malat1_qpcr_u95= upper(malat1_qpcr), 
            malat1_sc_mean=mean(sc_malat1_rel), 
            malat1_sc_l95=lower(sc_malat1_rel), 
            malat1_sc_u95=upper(sc_malat1_rel)) %>%
  ungroup() -> sc_vs_qpcr_smry

sc_vs_qpcr_smry %>%
  select(-tx, -reg, -color) -> sc_vs_qpcr_smry_out

write_supp_table(sc_vs_qpcr_smry_out, 'Bulk tissue target engagement quantified by single cell sequencing and by qPCR, summarized for all conditions.')

# binomial model results

afis = rbind_files('data/analytical/', 'afi_') %>%
  mutate(aggr = gsub('afi_','',file)) %>%
  mutate(region = gsub('(mouse|cyno)_','',gsub('_full','',aggr))) %>%
  mutate(species = gsub('_.*','',aggr)) %>%
  select(-file) %>%
  inner_join(select(meta, sample, animal_short), by=c('sample')) %>%
  select(-sample) %>%
  relocate(aggr, animal_short, species, region, model, celltype, subtype, txc)

afis %>%
  inner_join(select(tx_meta, tx, disp), by='tx') %>%
  select(-txc) %>%
  rename(treatment=disp) %>%
  mutate(model = gsub('aso7','ason',model)) %>%
  select(aggr, animal_short, species, region, treatment, weeks, model, celltype, subtype, h5_num, weeks, n_cells, total_umis, prnp_umis, malat1_umis, rnaseh1_umis, prnp_rpm, total_umis_pseudocount, fit_resid, color, x, y, point_estimate) -> afis_out 

write_supp_table(afis_out, 'Negative binomial model individual point estimates for all samples, conditions.')

acos = rbind_files('data/analytical/', 'aco_') %>%
  mutate(aggr = gsub('aco_','',file)) %>%
  mutate(region = gsub('(mouse|cyno)_','',gsub('_full','',aggr))) %>%
  mutate(species = gsub('_.*','',aggr)) %>%
  select(-file) %>%
  relocate(aggr, species, region, model, celltype, subtype, txc)

tx_meta %>% mutate(treatment = gsub('PBS|aCSF','control',disp)) %>% distinct(txc, treatment) -> txc_map
acos %>%
  inner_join(txc_map, by='txc') %>%
  select(-txc) %>%
  mutate(model = gsub('aso7','ason',model)) %>%
  mutate(weeks = as.integer(gsub('-.*','',model))) %>%
  select(aggr, species, region, treatment, weeks, model, celltype, subtype, coefname, estimate, se, z, p, log_normed, normed, l95, u95, color, x, y, n_cells, total_cells) -> acos_out

write_supp_table(acos_out, 'Negative binomial model coefficients for all conditions.')

acos %>%
  filter(txc != '_ctl') %>%
  inner_join(tx_meta, by=c('txc'='tx')) %>%
  rename(residual=normed, treatment=disp) %>%
  mutate(weeks = as.integer(gsub('-.*','',model))) %>%
  select(species, region, treatment, weeks, celltype, subtype, residual, l95, u95, total_cells) -> acos_minimal

write_supp_table(acos_minimal, 'Negative binomial model coefficients simplified summary.')

##############
# TABLE 2
##############

cat(file=stderr(), 'done.\nCreating Table 2...')

### summary table for Table 2
meta %>%
  inner_join(tx_meta, by='tx') %>%
  inner_join(region_meta, by='reg') %>%
  rename(treatment=disp) %>%
  group_by(species, treatment, region, weeks) %>%
  summarize(.groups='keep', n=n()) %>%
  ungroup() %>%
  arrange(desc(species), region, weeks, treatment) -> table_2

write_tsv(table_2, 'display_items/table-2.tsv')

zeroes = read_tsv('data/analytical/zeroes.tsv', col_types=cols())

write_supp_table(zeroes, 'Proportion of UMI count matrix cells that are zeroes.')

##############
# FIGURE 2
##############

cat(file=stderr(), 'done.\nCreating Figure 2...')

resx=300
png(paste0('display_items/figure-2.png'),width=6.5*resx,height=7.5*resx,res=resx)

layout_matrix = matrix(c(1,1,1,2,2,2,2,
                         seq(3,16,2),
                         seq(4,16,2),
                         seq(17,30,2),
                         seq(18,30,2),
                         c(seq(31,38,2),rep(39,3)),
                         c(seq(32,38,2),rep(39,3)),
                         40,40,40,41,41,41,41), nrow=8, byrow=T)
                         #45,45,45,46,46,46,46), nrow=8, byrow=T)
layout(layout_matrix,heights=c(2.5,rep(1,6),2.5),widths=c(.2,rep(1,6)))
#layout.show(41)
panel = 1

all_mouse_ana = rbind_files('data/analytical/','ana_mouse.*full') 
write(paste("Proportion of all UMIs that are Malat1 in all mouse samples: ",percent(sum(all_mouse_ana$malat1_umis) / sum(all_mouse_ana$total_umis),digits=1),'\n',sep=''),text_stats_path,append=T)

ana = read_tsv('data/analytical/ana_mouse_cerebellum_full.tsv', col_types=cols())
met = read_tsv('data/analytical/met_mouse_cerebellum_full.tsv', col_types=cols()) %>%
  mutate(order = row_number())
amt = read_csv('data/analytical/cla_mouse_cerebellum_full.csv', col_types=cols()) %>% 
  clean_names()

amt  %>%
  inner_join(ana, by='barcode') %>%
  inner_join(met, by=c('assignments')) %>%
  rename(subtype = assignments) %>%
  rename(subtype_color = color) %>%
  mutate(subtype_color = replace_na(subtype_color, '#A7A7A7')) -> model_data

aco = read_tsv(paste0('data/analytical/aco_mouse_cerebellum_full.tsv'), col_types=cols()) %>%
  filter(model=='12-asom')
afi = read_tsv(paste0('data/analytical/afi_mouse_cerebellum_full.tsv'), col_types=cols()) %>%
  filter(model=='12-asom')

par(mar=c(3,3,2.5,1))
xlims = c(0, 1.25)
ylims = c(0, 1.25)
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F)
axis(side=1, at=0:5/4, labels=percent(0:5/4), line=-1.1, lwd=0, cex.axis=0.7)
axis(side=1, at=0:5/4, labels=NA, tck=-0.02)
mtext(side=1, line=1, text=substitute(paste('bulk qPCR ', italic('Malat1'))), cex=0.6)
axis(side=2, at=0:5/4, labels=percent(0:5/4), las=2, line=-0.75, lwd=0, cex.axis=0.7)
axis(side=2, at=0:5/4, labels=NA, tck=-0.02)
mtext(side=2, line=1.6, text=substitute(paste('single cell ', italic('Malat1'))), cex=0.6)
abline(h=1, lty=3)
abline(v=1, lty=3)
abline(a=0, b=1, lwd=0.5, col='#A7A7A7')

indiv_subs = sc_vs_qpcr_indiv %>% filter(weeks==12 & reg=='cb' & tx %in% c('pbs0','asom'))
points(x=indiv_subs$malat1_qpcr, y=indiv_subs$sc_malat1_rel, pch=19, col=alpha(indiv_subs$color, ci_alpha))
smry_subs = sc_vs_qpcr_smry %>% filter(weeks==12 & reg=='cb' & tx %in% c('pbs0','asom'))
segments(x0=smry_subs$malat1_qpcr_l95, x1=smry_subs$malat1_qpcr_u95, y0=smry_subs$malat1_sc_mean, lwd=1.5, col=smry_subs$color)
segments(x0=smry_subs$malat1_qpcr_mean, y0=smry_subs$malat1_sc_l95, y1=smry_subs$malat1_sc_u95, lwd=1.5, col=smry_subs$color)

write(paste('Malat1 UMI/cell in PBS animals: mean ',formatC(mean(model_data$malat1_umis[model_data$tx=='pbs0']), format='f', digits=1),
            ' median: ',formatC(median(model_data$malat1_umis[model_data$tx=='pbs0']), format='f', digits=1),'\n',sep=''),text_stats_path,append=T)

write(paste('Prnp UMI/cell in PBS animals: mean ',formatC(mean(model_data$prnp_umis[model_data$tx=='pbs0']), format='f', digits=2),
            ' median: ',formatC(median(model_data$prnp_umis[model_data$tx=='pbs0']), format='f', digits=2),'\n',sep=''),text_stats_path,append=T)

write(paste('Total UMI/cell in PBS animals: mean ',formatC(mean(model_data$total_umis[model_data$tx=='pbs0']), format='f', digits=2),
            ' median: ',formatC(median(model_data$total_umis[model_data$tx=='pbs0']), format='f', digits=2),'\n',sep=''),text_stats_path,append=T)

write(paste('Prnp UMI/cell in ASO 6 animals: mean ',formatC(mean(model_data$prnp_umis[model_data$tx=='aso6']), format='f', digits=2),
            ' median: ',formatC(median(model_data$prnp_umis[model_data$tx=='aso6']), format='f', digits=2),'\n',sep=''),text_stats_path,append=T)

write(paste('Total UMI/cell in ASO 6 animals: mean ',formatC(mean(model_data$total_umis[model_data$tx=='aso6']), format='f', digits=2),
            ' median: ',formatC(median(model_data$total_umis[model_data$tx=='aso6']), format='f', digits=2),'\n',sep=''),text_stats_path,append=T)


par(xpd=T)
tx_meta %>%
  filter(tx %in% smry_subs$tx) -> leg
legend(x=.75, y=1.6 ,bty='n',legend=leg$disp,col=leg$color,pch=19,cex=0.75)
par(xpd=F)
mtext(LETTERS[panel], side=3, cex=1.5, adj = 0.0, line = 0.4)
panel = panel + 1
mtext('C', side=1, cex=1.5, adj = -0.2, line = 1.7)

this_subtype = 'all'
model_data %>%
  filter(tx %in% c('asom','pbs0')) %>%
  filter(weeks==12) %>%
  filter(this_subtype=='all' | subtype==this_subtype) -> this_subtype_data

binsize = 10

if (nrow(this_subtype_data) < 500) binsize = binsize*10
asom_hist = hist(this_subtype_data$malat1_umis[this_subtype_data$tx=='asom'], breaks=seq(0,1e5,binsize), plot=F)
pbs0_hist = hist(this_subtype_data$malat1_umis[this_subtype_data$tx=='pbs0'], breaks=seq(0,1e5,binsize), plot=F)

xmax = as.numeric(ceiling(quantile(this_subtype_data$malat1_umis[this_subtype_data$tx=='pbs0'], .9)/binsize)*binsize)
ymax = max(asom_hist$counts)

xlims = c(0,xmax*1.2)
xbigs = seq(0,xmax,binsize*5)
xats = seq(0,xmax,binsize)
ylims = c(0,ymax)
plot(NA, NA, xlim=xlims, ylim = ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
rect(xleft=pbs0_hist$mids-binsize/2, xright=pbs0_hist$mids+binsize/2, ybottom=rep(0,length(pbs0_hist)), ytop=pbs0_hist$counts, col=alpha(ctl_color, ci_alpha), border=NA)
rect(xleft=asom_hist$mids-binsize/2, xright=asom_hist$mids+binsize/2, ybottom=rep(0,length(asom_hist)), ytop=asom_hist$counts, col=alpha(mal_color, ci_alpha), border=NA)

axis(side=1, at=0:10*100, labels=NA, tck=-0.02)
axis(side=1, at=0:10*100, line=-1, lwd=0, cex.axis=0.6)
mtext(side=1, line=1.0, text=substitute(paste(italic('Malat1'),' UMIs per cell')), cex=0.6)
axis(side=2, labels=NA, tck=-0.025)
axis(side=2, line=-0.5, lwd=0, las=2, cex.axis=0.6)
mtext(side=2, line=1.75, text='cells', cex=0.6)

par(xpd=T)
arrows(x0=70,x1=50,y0=ymax,y1=980,angle=45,code=2,length=0.025,lwd=3)
par(xpd=F)
legend('topright', bty='n',legend=leg$disp,col=leg$color,pch=15,cex=0.75)
mtext(LETTERS[panel], side=3, cex=1.5, adj = 0.0, line = 0.4)
panel = panel + 1

cb12_smry = tibble(subtype=character(0),
                   malat1_resfit=numeric(0),
                   prnp_resfit=numeric(0),
                   n_cells=numeric(0))

lm_summary_tbl = tibble(subtype=character(0),
                        tx=character(0),
                        parameter=character(0),
                        estimate=numeric(0),
                        se=numeric(0),
                        t=numeric(0),
                        p=numeric(0))

last_supertype = '' # supertype means neuron, glia, other
for (this_subtype in met$assignments[met$assignments %in% model_data$subtype]) {
  # this_subtype = 'Golgi'
  
  this_supertype = met$celltype[met$assignments==this_subtype]
  if (this_supertype != last_supertype) {
    par(mar=c(0,0,0,0))
    plot(NA, NA, xlim=0:1, ylim=0:1, xaxs='i', yaxs='i', axes=F, ann=F)
    par(xpd=T)
    mtext(side=4, padj=0, line=-1.5, text=this_supertype, at=0, cex=0.8)
    par(xpd=F)
    plot(NA, NA, xlim=0:1, ylim=0:1, xaxs='i', yaxs='i', axes=F, ann=F)
    last_supertype = this_supertype
  }
  
  model_data %>%
    filter(tx %in% c('asom','pbs0')) %>%
    filter(weeks==12) %>%
    filter(this_subtype=='all' | subtype==this_subtype) -> this_subtype_data
  
    binsize = max(10^floor(log10(max(this_subtype_data$malat1_umis))-2),10)
    par(mar=c(0,0,1.5,0))
  
  if (nrow(this_subtype_data) < 225) binsize = binsize*10
  asom_hist = hist(this_subtype_data$malat1_umis[this_subtype_data$tx=='asom'], breaks=seq(0,1e5,binsize), plot=F)
  pbs0_hist = hist(this_subtype_data$malat1_umis[this_subtype_data$tx=='pbs0'], breaks=seq(0,1e5,binsize), plot=F)
  
  xmax = as.numeric(ceiling(quantile(this_subtype_data$malat1_umis[this_subtype_data$tx=='pbs0'], .9)/binsize)*binsize)
  ymax = max(asom_hist$counts)
  
  xlims = c(0,xmax*1.2)
  xbigs = seq(0,xmax,binsize*5)
  xats = seq(0,xmax,binsize)
  ylims = c(0,ymax)
  #yincrement =  max(10^floor(log10(ymax)-1),1)
  #ybigs = seq(0, ymax,yincrement)
  plot(NA, NA, xlim=xlims, ylim = ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
  axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
  rect(xleft=pbs0_hist$mids-binsize/2, xright=pbs0_hist$mids+binsize/2, ybottom=rep(0,length(pbs0_hist)), ytop=pbs0_hist$counts, col=alpha(ctl_color, ci_alpha), border=NA)
  rect(xleft=asom_hist$mids-binsize/2, xright=asom_hist$mids+binsize/2, ybottom=rep(0,length(asom_hist)), ytop=asom_hist$counts, col=alpha(mal_color, ci_alpha), border=NA)
  residual = mean(this_subtype_data$malat1_umis[this_subtype_data$tx=='asom'])/mean(this_subtype_data$malat1_umis[this_subtype_data$tx=='pbs0'])
  residual_adjusted = mean(this_subtype_data$malat1_umis[this_subtype_data$tx=='asom']/this_subtype_data$total_umis[this_subtype_data$tx=='asom'])/
    mean(this_subtype_data$malat1_umis[this_subtype_data$tx=='pbs0']/this_subtype_data$total_umis[this_subtype_data$tx=='pbs0'])
  residual_fit = aco$normed[aco$subtype==this_subtype & aco$txc=='asom']
  
  axis(side=1, labels=NA, tck=0.025)
  axis(side=1, line=-2, padj=0, lwd=0, cex.axis=0.3)
  axis(side=2, labels=NA, tck=0.025)
  axis(side=2, line=-1.25, hadj=0, lwd=0, las=2, cex.axis=0.3)
  mtext(side=3, line=0, text=paste0(' ',this_subtype, ' ', percent(round(residual_fit*100)/100)), cex=0.6)
  axis(side=4, at=ylims, lwd.ticks=0, labels=NA)
  
  model_data %>%
    filter(tx %in% c('asom','pbs0','aso6')) %>%
    filter(weeks==12) %>%
    filter(this_subtype=='all' | subtype==this_subtype) -> this_subtype_both_asos
  n_cells_total = nrow(this_subtype_both_asos)
  cb12_smry = rbind(cb12_smry, cbind.data.frame(subtype=this_subtype, malat1_resfit=residual_fit, prnp_resfit=as.numeric(NA), n_cells=n_cells_total))
  
  model_data %>%
    filter(tx %in% c('asom','pbs0')) %>%
    filter(weeks==12) %>%
    filter(this_subtype=='all' | subtype==this_subtype) %>%
    inner_join(tx_meta, by='tx') %>%
    rename(tx_color=color) -> this_subtype_data
  
  par(mar=c(1.5,0,0,0))
  ymax = as.numeric(ceiling(quantile(this_subtype_data$malat1_umis[this_subtype_data$tx=='pbs0'], .9)/binsize)*binsize)
  xmax = as.numeric(ceiling(quantile(this_subtype_data$total_umis[this_subtype_data$tx=='pbs0'], .9)/binsize)*binsize)
  xlims = c(0,xmax)
  ylims = c(0,ymax)
  plot(NA, NA, xlim = xlims, ylim = ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  points(x=this_subtype_data$total_umis, this_subtype_data$malat1_umis, col=alpha(this_subtype_data$tx_color,ci_alpha), pch=20, cex=0.1)
  for (this_tx in unique(this_subtype_data$tx)) {
    subs = subset(this_subtype_data, tx==this_tx)
    m = lm(malat1_umis ~ total_umis, data=subs)
    abline(m, col=subs$tx_color[1])
    
    summary(m)$coefficients %>% 
      as_tibble(rownames='varname') %>%
      rename(estimate=Estimate, se=`Std. Error`, t=`t value`, p=`Pr(>|t|)`) %>%
      mutate(subtype=this_subtype, 
             tx=this_tx) %>%
      mutate(parameter = case_when(varname=='(Intercept)' ~ 'intercept',
                                   varname=='total_umis' ~ 'slope')) %>%
      select(subtype, tx, parameter, estimate, se, t, p) -> this_lm_summary
    lm_summary_tbl = rbind(lm_summary_tbl,
                           this_lm_summary)
  }
  axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
  axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
  axis(side=1, labels=NA, tck=0.025)
  axis(side=1, line=-2, padj=0, lwd=0, cex.axis=0.3)
  axis(side=2, labels=NA, tck=0.025)
  axis(side=2, line=-1.25, hadj=0, lwd=0, las=2, cex.axis=0.3)
  axis(side=4, at=ylims, lwd.ticks=0, labels=NA)
  
  if (this_subtype=='microglia') {
    par(mar=c(0,1,1.5,0))
    plot(NA, NA, xlim=c(0,1), ylim=c(0,1), axes=F, ann=F, xaxs='i', yaxs='i')
    mtext(side=3, line=0, text='key', font=3, cex=0.6)
    mtext(side=3, line=-0.05, text='___', font=3, cex=0.6)
    axis(side=1, labels=NA, tck=0.025)
    axis(side=2, labels=NA, tck=0.025)
    mtext(side=1, line=-1, text=expression(paste(italic('Malat1'),' UMIs')), cex=0.45)
    mtext(side=2, line=-0.75, text='cells', cex=0.45)
    legend('topright', bty='n',legend=leg$disp,col=leg$color,pch=15,cex=0.75)
    #box('plot',lwd=1.25)
    
    par(mar=c(1.5,1,0,0))
    plot(NA, NA, xlim=c(0,1), ylim=c(0,1), axes=F, ann=F, xaxs='i', yaxs='i')
    segments(x0=0,y0=0,x1=.5,y1=.75,col=leg$color[leg$tx=='pbs0'])
    segments(x0=0,y0=0,x1=.5,y1=.38,col=leg$color[leg$tx=='asom'])
    text(x=c(.5,.5),y=c(.75,.38),labels=c('PBS fit','ASO fit'),pos=4,col=leg$color[match(leg$tx,c('pbs0','asom'))],cex=0.6)
    axis(side=1, labels=NA, tck=0.025)
    axis(side=2, labels=NA, tck=0.025)
    mtext(side=1, line=-1, text='total UMIs', cex=0.45)
    mtext(side=2, line=-0.75, text=expression(paste(italic('Malat1'),' UMIs')), cex=0.45)
    
  } 
  
}

write_supp_table(lm_summary_tbl, 'Parameters of linear models malat1_umis ~ total_umis by cell type.')

# aco/afi plots for Malat1
aco %>%
  distinct(x, y, celltype, subtype) -> xyleg
xyleg %>%
  group_by(celltype) %>%
  summarize(.groups='keep',
            miny=min(y),
            maxy=max(y),
            midy = (min(y)+max(y))/2,
            minx=min(x),
            maxx=max(x),
            midx = (min(x)+max(x))/2) -> tranches
par(mar=c(4,3,1.5,0.5))
ylims = c(0, 1.4)
xlims = range(afi$x) + c(-0.5, 0.5)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=2, at=0:8/4, labels=NA, tck=-0.02)
axis(side=2, at=0:8/4, labels=percent(0:8/4), lwd=0, line=-0.75, las=2,cex.axis=0.6)
mtext(side=2, line=1.6, text=expression(paste('residual ',italic('Malat1'))), cex=0.6)
abline(h=1, lwd=0.25, lty=3)
abline(h=0.5, lwd=0.125, lty=3)
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
par(xpd=T)
text(x=xyleg$x+0.05, y=rep(-0.07,nrow(xyleg)), labels=xyleg$subtype, srt=75, adj=1, cex=0.6)
par(xpd=F)
tranche_line = 4.0
overhang_left = 0.5
overhang_right = 0.1
for (i in 1:nrow(tranches)) {
  axis(side=1, line=tranche_line, at=c(tranches$minx[i], tranches$maxx[i]) + c(-1,1)*c(overhang_left,overhang_right), tck=0.03, labels=NA)
  mtext(side=1, line=tranche_line, at=tranches$midx[i], text=tranches$celltype[i], cex=0.6)
}
points(x=afi$x, y=afi$point_estimate, pch=19, col=alpha(afi$color, ci_alpha))
barwidth=0.4
segments(x0=aco$x-barwidth, x1=aco$x+barwidth, y0=aco$normed, col=aco$color, lwd=1)
arrows(x0=aco$x, y0=aco$l95, y1=aco$u95, angle=90, length=0.03, code=3, col=aco$color, lwd=0.5)
panel = 4
mtext(LETTERS[panel], side=3, cex=1.5, adj = 0.0, line = 0.3)
panel = panel + 1

aco_prnp = read_tsv(paste0('data/analytical/aco_mouse_cerebellum_full.tsv'), col_types=cols()) %>%
  filter(model=='12-aso6') %>% filter(txc=='aso6')
cb12_smry$prnp_resfit = aco_prnp$normed[match(cb12_smry$subtype, aco_prnp$subtype)]

met = read_tsv('data/analytical/met_mouse_cerebellum_full.tsv',col_types=cols())
cb12_smry$color = met$color[match(cb12_smry$subtype, met$assignments)]

par(mar=c(3,3,2.5,1))
xlims = c(0,1.05)
ylims = c(0,1.05)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=0:5/4, labels=NA, tck=-0.02)
axis(side=1, at=0:5/4, labels=percent(0:5/4), lwd=0, line=-0.75, cex.axis=0.8)
abline(h=1, lty=3)
abline(v=1, lty=3)
axis(side=2, at=0:5/4, labels=NA, tck=-0.02)
axis(side=2, at=0:5/4, labels=percent(0:5/4), las=2, lwd=0, line=-0.75, cex.axis=0.8)
points(x=cb12_smry$malat1_resfit, y=cb12_smry$prnp_resfit, pch=19, col=cb12_smry$color, cex=log10(cb12_smry$n_cells)-2)
mtext(side=1, line=1.5, text=substitute(paste(italic('Malat1'),' residual')), cex=0.8)
mtext(side=2, line=1.6, text=substitute(paste(italic('Prnp'),' residual')), cex=0.8)
m = lm(prnp_resfit ~ malat1_resfit, weights=n_cells, data=cb12_smry)
abline(m, col='#A7A7A7')
mtext(LETTERS[panel], side=3, cex=1.5, adj = 0.0, line = 0.4)
panel = panel + 1

malat1_prnp_wtd_cor = wtd.cor(cb12_smry$malat1_resfit, cb12_smry$prnp_resfit, weight = cb12_smry$n_cells)

write(paste("Malat1 vs. Prnp residual in 12 week cerebellum weighted Pearson's correlation: rho = ",formatC(malat1_prnp_wtd_cor['Y','correlation'],format='f',digits=2),' P = ',formatC(malat1_prnp_wtd_cor['Y','p.value'],format='e',digits=1),'\n',sep=''),text_stats_path,append=T)

ctl_color = '#A7A7A7'
aso6_color = '#41B6C4'
fit1_color = '#373737'
hypoth1_color = '#D43D1A'
hypoth2_color = '#9F79EE'

this_subtype = 'astrocyte'

model_data %>%
  filter(tx %in% c('aso6','pbs0')) %>%
  filter(weeks==12) %>%
  filter(this_subtype=='all' | subtype==this_subtype) -> this_subtype_data

par(mar=c(3,3,2.5,1))

prnp_breaks = c(0:10,100)

aso6_hist = hist(this_subtype_data$prnp_umis[this_subtype_data$tx=='aso6'], breaks=prnp_breaks, include.lowest=T, right=F, plot=F)
pbs0_hist = hist(this_subtype_data$prnp_umis[this_subtype_data$tx=='pbs0'], breaks=prnp_breaks, include.lowest=T, right=F, plot=F)

xmax = 5
ymax = 0.6
xats = 0:10
yats = 1:3/4
xlims = c(0,xmax)
ylims = c(0,ymax)
plot(NA, NA, xlim=xlims, ylim = ylims, axes=F, ann=F, xaxs='i', yaxs='i')
offset_aso6 = -0.15
offset_pbs0 = -0.3
hwidth = 0.15
rect(xleft=pbs0_hist$mids+offset_pbs0-hwidth/2, xright=pbs0_hist$mids+offset_pbs0+hwidth/2, ybottom=rep(0,length(pbs0_hist$density)), ytop=pbs0_hist$density, col=ctl_color , border=NA)
rect(xleft=aso6_hist$mids+offset_aso6-hwidth/2, xright=aso6_hist$mids+offset_aso6+hwidth/2, ybottom=rep(0,length(aso6_hist$density)), ytop=aso6_hist$density, col=aso6_color, border=NA)
residual = mean(this_subtype_data$prnp_umis[this_subtype_data$tx=='aso6'])/mean(this_subtype_data$prnp_umis[this_subtype_data$tx=='pbs0'])
residual_adjusted = mean(this_subtype_data$prnp_umis[this_subtype_data$tx=='aso6']/this_subtype_data$total_umis[this_subtype_data$tx=='aso6'])/
  mean(this_subtype_data$prnp_umis[this_subtype_data$tx=='pbs0']/this_subtype_data$total_umis[this_subtype_data$tx=='pbs0'])

pbs_nb = glm.nb(prnp_umis ~ 1, data=subset(this_subtype_data, tx=='pbs0'))

set.seed(1)
binom_sampling_n = 100000
nb_pbs = rnbinom(n=binom_sampling_n, size=pbs_nb$theta, mu=fitted(pbs_nb)[1])
# mean(nb_pbs) / mean(this_subtype_data$prnp_umis[this_subtype_data$tx=='pbs0']) # check that this is ~1

fit1 = nb_pbs
fit1_hist = hist(fit1, breaks=prnp_breaks, include.lowest=T, right=F, plot=F)

hypoth1 = c(sample(this_subtype_data$prnp_umis[this_subtype_data$tx=='pbs0'], replace=T, size = round(binom_sampling_n*residual)),
            rep(0,round(binom_sampling_n*(1-residual))))
hypoth1_hist = hist(hypoth1, breaks=prnp_breaks, include.lowest=T, right=F, plot=F)
# this should be similar to residual without adjustment for total UMIs (which is what is plotted in the histo)
hypoth1_residual = mean(hypoth1)/mean(this_subtype_data$prnp_umis[this_subtype_data$tx=='pbs0'])
hypoth2 = rnbinom(n=binom_sampling_n, size=pbs_nb$theta, mu=fitted(pbs_nb)[1]*residual)
hypoth2_hist = hist(hypoth2, breaks=prnp_breaks, include.lowest=T, right=F, plot=F)
hypoth2_residual = mean(hypoth2)/mean(this_subtype_data$prnp_umis[this_subtype_data$tx=='pbs0'])

write(paste('Residual Prnp under hypothesis 1 and 2 respectively: ',percent(hypoth1_residual,digits=1),' and ',percent(hypoth1_residual,digits=1),'\n',sep=''),text_stats_path,append=T)

offset_fit1 = 0.0
offset_hypoth1 = 0.15
offset_hypoth2 = 0.3

rect(xleft=fit1_hist$mids+offset_fit1-hwidth/2, xright=fit1_hist$mids+offset_fit1+hwidth/2, ybottom=rep(0,length(fit1_hist$density)), ytop=fit1_hist$density, col=fit1_color, density=30)
rect(xleft=hypoth1_hist$mids+offset_hypoth1-hwidth/2, xright=hypoth1_hist$mids+offset_hypoth1+hwidth/2, ybottom=rep(0,length(hypoth1_hist$density)), ytop=hypoth1_hist$density, col=hypoth1_color, density=30)
rect(xleft=hypoth2_hist$mids+offset_hypoth2-hwidth/2, xright=hypoth2_hist$mids+offset_hypoth2+hwidth/2, ybottom=rep(0,length(hypoth2_hist$density)), ytop=hypoth2_hist$density, col=hypoth2_color, density=30)

axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
axis(side=1, at=xats, labels=NA, tck=-0.025)
mtext(side=1, at=xats+0.5, text=xats, line=0.25, cex=0.8)
axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
axis(side=2, at=yats, labels=NA, tck=-0.025)
axis(side=2, at=yats, labels=percent(yats), line=-0.5, lwd=0, las=2)
#mtext(side=2, at=yats, text=percent(yats), line=-1, cex=0.4, las=2)
mtext(side=1, line=1.5, text=substitute(paste(italic('Prnp'),' UMIs per cell')), cex=0.75)
mtext(side=2, line=2.5, text='% of cells', cex=0.75)

par(xpd=T)
legend(x=2.2, y=.62,
       c('PBS actual','ASO 6 actual','NB fit for PBS','Hypothesis 1: 100% knockdown in subset of cells','Hypothesis 2: equal knockdown in all cells'),
       fill=c(ctl_color, aso6_color, fit1_color, hypoth1_color, hypoth2_color),
       border=c(ctl_color, aso6_color, fit1_color, hypoth1_color, hypoth2_color),
       density=c(NA,NA,30,30,30), cex=0.75)
par(xpd=F)
mtext(LETTERS[panel], side=3, cex=1.5, adj = 0.0, line = 0.4)

hist_models = tibble(objname = c('pbs0_hist', 'aso6_hist', 
                                 'fit1_hist', 'hypoth1_hist', 'hypoth2_hist'),
                     description = c('PBS actual','ASO 6 actual','NB fit for PBS','Hypothesis 1: 100% knockdown in subset of cells','Hypothesis 2: equal knockdown in all cells'))
if (exists('histogram_summary')) rm('histogram_summary')
for (i in 1:nrow(hist_models)) {
  histogram = get(hist_models$objname[i])
  this_tbl = cbind(model=rep(hist_models$description[i],11),
                   umis = histogram$breaks[1:11],
                   count = histogram$counts[1:11],
                   density = histogram$density[1:11])
  if (exists('histogram_summary')) {
    histogram_summary = rbind(histogram_summary, this_tbl)
  } else {
    histogram_summary = this_tbl
  }
}
histogram_summary %>% as_tibble() -> histogram_summary

write_supp_table(histogram_summary, "Histograms of Prnp UMI counts in astrocytes under different models.")

unnecessary_message_end_of_figure_2 = dev.off()


##############
# FIGURE 3
##############

cat(file=stderr(), 'done.\nCreating Figure 3...')

resx=300
png(paste0('display_items/figure-3.png'),width=6.5*resx,height=7*resx,res=resx)

layout_matrix = matrix(c(1,2,2,3,3,3,4,4,
                         5,5,5,5,8,8,8,8,
                         5,5,5,5,9,9,9,9,
                         6,6,6,6,10,10,10,10,
                         6,6,6,6,11,11,11,11,
                         7,7,7,7,12,13,14,15,
                         7,7,7,7,12,13,14,15), nrow=7, byrow=T)
layout(layout_matrix, 
       widths=c(.5,1.25,1.25,1,.5,1,1.25,1.25),
       heights=c(1.3,1,1,1,1,1.2,1.2))
#layout.show(15)

panel = 1

par(mar=c(2,0,2,0))
xlims = c(0, 1.25)
ylims = c(0, 1.25)
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F)
axis(side=4, at=0:5/4, labels=percent(0:5/4), las=2, line=-1.25, hadj=1, lwd=0, cex.axis=0.7)
axis(side=4, at=0:5/4, labels=NA, tck=0.02)
mtext(side=4, line=-2.75, text=substitute(paste('single cell ', italic('Prnp'))), cex=0.6)


for (aggr in c('mouse_cortex_full','mouse_thalamus_full','mouse_cerebellum_full')) {

  this_reg = case_when(aggr=='mouse_cortex_full' ~ 'ss',
                  aggr=='mouse_thalamus_full' ~ 'th',
                  aggr=='mouse_cerebellum_full' ~ 'cb')
  this_reg_disp = gsub('mouse_','',gsub('_full','',aggr))

  par(mar=c(2,0,2,2))
  xlims = c(0, 1.25)
  ylims = c(0, 1.25)
  plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F)
  axis(side=1, at=0:5/4, labels=percent(0:5/4), line=-1.1, lwd=0, cex.axis=0.7)
  axis(side=1, at=0:5/4, labels=NA, tck=-0.02)
  mtext(side=1, line=1, text=substitute(paste('bulk qPCR ', italic('Prnp'))), cex=0.6)
  # axis(side=2, at=0:5/4, labels=percent(0:5/4), las=2, line=-0.75, lwd=0, cex.axis=0.7)
  axis(side=2, at=0:5/4, labels=NA, tck=-0.02)
  # mtext(side=2, line=1.4, text=substitute(paste('single cell ', italic('Prnp'))), cex=0.6)
  abline(h=1, lty=3)
  abline(v=1, lty=3)
  abline(a=0, b=1, lwd=0.5, col='#A7A7A7')

  indiv_subs = sc_vs_qpcr_indiv %>% filter(weeks==3 & reg==this_reg & tx %in% c('pbs0','aso6'))
  points(x=indiv_subs$prnp_qpcr, y=indiv_subs$sc_prnp_rel, pch=19, col=alpha(indiv_subs$color, ci_alpha))
  smry_subs = sc_vs_qpcr_smry %>% filter(weeks==3 & reg==this_reg & tx %in% c('pbs0','aso6'))
  segments(x0=smry_subs$prnp_qpcr_l95,  x1=smry_subs$prnp_qpcr_u95, y0=smry_subs$prnp_sc_mean, lwd=1.5, col=smry_subs$color)
  segments(x0=smry_subs$prnp_qpcr_mean, y0=smry_subs$prnp_sc_l95,   y1=smry_subs$prnp_sc_u95,  lwd=1.5, col=smry_subs$color)
  
  mtext(side=3, line=0, text=this_reg_disp, cex=0.6)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.1, line = 0.1)
  panel = panel + 1
}

aggr = 'mouse_cortex_full'
read_tsv(paste0('data/analytical/aco_',aggr,'.tsv'), col_types=cols()) %>%
  distinct(y, celltype, subtype) -> yleg
yleg %>%
  group_by(celltype) %>%
  summarize(.groups='keep',
            miny=min(y),
            maxy=max(y),
            midy = (min(y)+max(y))/2) -> tranches

for (aggr in c('mouse_cortex_full','mouse_thalamus_full','mouse_cerebellum_full')) {
  
  this_tx = 'aso6'
  this_model = '3-aso6'
  
  afi = read_tsv(paste0('data/analytical/afi_',aggr,'.tsv'), col_types=cols()) %>%
    filter(model==this_model)
  aco = read_tsv(paste0('data/analytical/aco_',aggr,'.tsv'), col_types=cols()) %>%
    filter(model==this_model) 
  
  aco %>%
    distinct(y, celltype, subtype) -> yleg
  yleg %>%
    group_by(celltype) %>%
    summarize(.groups='keep',
              miny=min(y),
              maxy=max(y),
              midy = (min(y)+max(y))/2) -> tranches
  
  if (aggr=='mouse_cortex_full') {
    par(mar=c(0.5,11,2,1))
  } else if (aggr=='mouse_thalamus_full') {
    par(mar=c(0.5,11,0,1))
  } else if (aggr=='mouse_cerebellum_full') {
    par(mar=c(3,11,0,1))
  }
  
  xlims = c(0, 1.5)
  ylims = range(afi$y) + c(-0.5, 0.5)
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=0:8/4, labels=NA, tck=-0.02)
  abline(v=1, lwd=0.25, lty=3)
  abline(v=0.5, lwd=0.125, lty=3)
  axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
  mtext(side=2, at=yleg$y, text=yleg$subtype, line=0.25, las=2, cex=0.75)
  tranche_line = 8
  overhang = 0.3
  for (i in 1:nrow(tranches)) {
    axis(side=2, line=tranche_line, at=c(tranches$miny[i], tranches$maxy[i]) + c(-1,1)*overhang, tck=0.05, labels=NA)
    mtext(side=2, line=tranche_line + 0.25, at=tranches$midy[i], text=tranches$celltype[i], cex=0.8)
  }
  points(y=afi$y, x=afi$point_estimate, pch=19, col=alpha(afi$color, ci_alpha))
  barwidth=0.4
  segments(y0=aco$y-barwidth, y1=aco$y+barwidth, x0=aco$normed, col=aco$color, lwd=2)
  arrows(y0=aco$y, x0=aco$l95, x1=aco$u95, angle=90, length=0.03, code=3, col=aco$color, lwd=2)
  
  reg = gsub('mouse_','',gsub('_full','',aggr))
  mtext(side=2, line=tranche_line+1.25, text=reg, font=1, cex=1)
  
  if (aggr=='mouse_cortex_full') {
    mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.1, line = 0.2)
    panel = panel + 1
  } else if (aggr=='mouse_cerebellum_full') {
    axis(side=1, at=0:8/4, labels=percent(0:8/4), lwd=0, line=-0.75)
    mtext(side=1, line=1.6, text=substitute(paste('residual ',italic('Prnp'))), cex=0.7)
  }
}

exc_clust = read_tsv('data/analytical/mouse_exc_clust.tsv', col_types=cols())
exc_meta = tibble(cluster=c(10,15,4,7,3,8,13,6,5,2,9,14,12,11,1),
                  x = 1:15)

exc_marker = tibble(gene=c('Cux2','Rorb','Sulf2','Foxp2'),
                    layer = c('L2/3','L4','L5','L6'),
                    color=c("#FF7777","#FFAA00","#00BB66","#6600FF"))
exc_clust %>%
  rename(gene=name) %>%
  inner_join(exc_meta, by='cluster') %>%
  group_by(gene) %>%
  mutate(pmax = rpm/max(rpm)) %>%
  ungroup() %>%
  arrange(x) %>%
  inner_join(exc_marker, by='gene') -> exc

exc_barcodes = read_csv('data/analytical/ex6_recluster.csv', 
                        col_types=cols()) %>% 
  clean_names() %>%
  mutate(cluster = as.numeric(gsub('Cluster ','',ex6_recluster))) %>%
  select(barcode, cluster)
ana = read_tsv('data/analytical/ana_mouse_cortex_full.tsv', col_types=cols())
amt = read_csv('data/analytical/cla_mouse_cortex_full.csv', col_types=cols()) %>% clean_names()

ana %>%
  select(-cluster) %>%
  filter(weeks==3) %>%
  inner_join(amt, by='barcode') %>%
  filter(assignments=='excitatory') %>%
  inner_join(exc_barcodes, by='barcode') %>%
  group_by(animal_short, tx, cluster) %>%
  summarize(.groups='keep',
            prnp_rpm = 1e6*mean(prnp_umis/total_umis, na.rm=T)) %>%
  ungroup() %>%
  group_by(cluster) %>%
  mutate(saline_mean = mean(prnp_rpm[tx=='pbs0'])) %>%
  mutate(rel = prnp_rpm/saline_mean) %>%
  ungroup() %>%
  inner_join(exc_meta, by='cluster') %>%
  inner_join(select(tx_meta, tx, color), by='tx') -> exc_residual_indivs

exc_residual_indivs %>%
  group_by(tx, cluster, color) %>%
  summarize(.groups='keep',
            residual = mean(rel),
            l95 = lower(rel),
            u95 = upper(rel)) %>%
  ungroup() %>%
  inner_join(exc_meta, by='cluster') %>%
  arrange(x) -> exc_residual_smry

exc_tx_meta = tx_meta %>% filter(tx %in% exc_residual_indivs$tx)

par(mar=c(0,5,2,6))
xlims = c(0.5, 15.5)
ylims = c(0, 1.5)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
axis(side=1, exc_meta$x, labels=NA, tck=-0.02)
axis(side=2, at=0:6/4, labels=NA, tck=-0.02)
axis(side=2, at=0:6/4, labels=percent(0:6/4), las=2, lwd=0, line=-0.75)
mtext(side=2, line=2.5, text='Prnp', cex=0.7, font=3)
abline(h=1, lty=3)
points(x=exc_residual_indivs$x, y=exc_residual_indivs$rel, pch=19, col=alpha(exc_residual_indivs$color,ci_alpha))
for (this_tx in c('aso6','pbs0')) {
  subs = subset(exc_residual_smry, tx==this_tx)
  par(xpd=T)
  points( x=  subs$x, y=  subs$residual, type='l', lwd=2, col=subs$color)
  polygon(x=c(subs$x, rev(subs$x)), y=c(subs$l95, rev(subs$u95)), col=alpha(subs$color, ci_alpha), border=NA)
  par(xpd=F)
}
par(xpd=T)
legend(x=max(xlims), y=max(ylims), exc_tx_meta$disp, col=exc_tx_meta$color, pch=19, lwd=2, bty='n')
par(xpd=F)
mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.15, line = 0.4)
panel = panel + 1

par(mar=c(2,5,1,6))
xlims = c(0.5, 15.5)
ylims = c(0, 1.05)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, lwd.ticks=0, labels=NA)
axis(side=1, exc_meta$x, labels=NA, tck=-0.02)
# axis(side=3, exc_meta$x, labels=exc_meta$cluster, line=-1, lwd=0, cex.axis=0.5)
axis(side=2, at=0:5/4, labels=NA, tck=-0.02)
axis(side=2, at=0:5/4, labels=percent(0:5/4), las=2, lwd=0, line=-0.75)
mtext(side=2, line=2.5, text='marker', cex=0.7)
for (mkr in exc_marker$gene) {
  subs = subset(exc, gene==mkr)
  x = seq(1,15,.1)
  points(x=x, y=predict(loess(pmax ~ x, data=subs, span=.5), x), type='l', lwd=2, col=subs$color)
  points(subs$x, subs$pmax, pch=19, col=alpha(subs$color,ci_alpha))
}
par(xpd=T)
mtext(side=1, line=0.5, at=c(3,6,11,14.5), text=exc_marker$layer, col=exc_marker$color, cex=0.9)
legend(x=max(xlims), y=max(ylims), exc_marker$gene, col=exc_marker$color, lwd=2, bty='n', text.font=3)
par(xpd=F)


exc %>%
  select(cluster, gene, rpm, n_cells) %>%
  pivot_wider(values_from=rpm, names_from=gene) %>%
  inner_join(filter(exc_residual_smry, tx=='aso6'), by='cluster') %>%
  select(x, cluster, n_cells,
         cux2_rpm = Cux2, rorb_rpm = Rorb, sulf2_rpm = Sulf2, foxp2_rpm = Foxp2,
         prnp_residual = residual, l95, u95) -> exc_supp
write_supp_table(exc_supp, 'Cortical excitatory neuron clusters: marker expression and residual Prnp.')

inh_clust = read_tsv('data/analytical/mouse_inh_clust.tsv', col_types=cols())
inh_meta = tibble(cluster=c(12,9,7,10,1,3,4,5,6,11,8,13,2),
                  x = 1:13)

inh_marker = tibble(gene=c('Vip','Sst','Pvalb','Npy'),
                    color=c("#a29d37","#b01d24","#223d7f","#0d92c5"))
inh_clust %>%
  rename(gene=name) %>%
  inner_join(inh_meta, by='cluster') %>%
  group_by(gene) %>%
  mutate(pmax = rpm/max(rpm)) %>%
  ungroup() %>%
  arrange(x) %>%
  inner_join(inh_marker, by='gene') -> inh

inh_barcodes = read_csv('data/analytical/in6_recluster.csv', 
                        col_types=cols()) %>% 
  clean_names() %>%
  mutate(cluster = as.numeric(gsub('Cluster ','',in6_recluster))) %>%
  select(barcode, cluster)

ana = read_tsv('data/analytical/ana_mouse_cortex_full.tsv', col_types=cols())
amt = read_csv('data/analytical/cla_mouse_cortex_full.csv', col_types=cols()) %>% clean_names()


ana %>%
  select(-cluster) %>%
  filter(weeks==3) %>%
  inner_join(amt, by='barcode') %>%
  filter(assignments=='inhibitory') %>%
  inner_join(inh_barcodes, by='barcode') %>%
  group_by(animal_short, tx, cluster) %>%
  summarize(.groups='keep',
            prnp_rpm = 1e6*mean(prnp_umis/total_umis, na.rm=T)) %>%
  ungroup() %>%
  group_by(cluster) %>%
  mutate(saline_mean = mean(prnp_rpm[tx=='pbs0'])) %>%
  mutate(rel = prnp_rpm/saline_mean) %>%
  ungroup() %>%
  inner_join(inh_meta, by='cluster') %>%
  inner_join(select(tx_meta, tx, color), by='tx') -> inh_residual_indivs

inh_residual_indivs %>%
  group_by(tx, cluster, color) %>%
  summarize(.groups='keep',
            residual = mean(rel),
            l95 = lower(rel),
            u95 = upper(rel)) %>%
  ungroup() %>%
  inner_join(inh_meta, by='cluster') %>%
  arrange(x) -> inh_residual_smry

inh_tx_meta = tx_meta %>% filter(tx %in% inh_residual_indivs$tx)


par(mar=c(0,5,2,6))
xlims = c(0.5, 13.5)
ylims = c(0, 1.5)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, lwd.ticks=0, labels=NA)
axis(side=1, inh_meta$x, labels=NA, tck=-0.02)
axis(side=2, at=0:6/4, labels=NA, tck=-0.02)
axis(side=2, at=0:6/4, labels=percent(0:6/4), las=2, lwd=0, line=-0.75)
mtext(side=2, line=2.5, text='Prnp', cex=0.7, font=3)
abline(h=1, lty=3)
points(x=inh_residual_indivs$x, y=inh_residual_indivs$rel, pch=19, col=alpha(inh_residual_indivs$color,ci_alpha))
for (this_tx in c('aso6','pbs0')) {
  subs = subset(inh_residual_smry, tx==this_tx)
  par(xpd=T)
  points( x=  subs$x, y=  subs$residual, type='l', lwd=2, col=subs$color)
  polygon(x=c(subs$x, rev(subs$x)), y=c(subs$l95, rev(subs$u95)), col=alpha(subs$color, ci_alpha), border=NA)
  par(xpd=F)
}
par(xpd=T)
legend(x=max(xlims), y=max(ylims), inh_tx_meta$disp, col=inh_tx_meta$color, pch=19, lwd=2, bty='n')
par(xpd=F)
mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.15, line = 0.5)
panel = panel + 1


par(mar=c(2,5,1,6))
xlims = c(0.5, 13.5)
ylims = c(0, 1.05)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, lwd.ticks=0, labels=NA)
axis(side=1, inh_meta$x, labels=NA, tck=-0.02)
axis(side=2, at=0:5/4, labels=NA, tck=-0.02)
axis(side=2, at=0:5/4, labels=percent(0:5/4), las=2, lwd=0, line=-0.75)
mtext(side=2, line=2.5, text='marker', cex=0.7)
for (mkr in inh_marker$gene) {
  subs = subset(inh, gene==mkr)
  x = seq(1,15,.1)
  points(x=x, y=predict(loess(pmax ~ x, data=subs, span=.5), x), type='l', lwd=2, col=subs$color)
  points(subs$x, subs$pmax, pch=19, col=alpha(subs$color,ci_alpha))
}
par(xpd=T)
mtext(side=1, line=0.5, at=c(2.5,5.5,10,13), text=inh_marker$gene, col=inh_marker$color, cex=0.9, font=3)
legend(x=max(xlims), y=max(ylims), inh_marker$gene, col=inh_marker$color, lwd=2, bty='n',text.font=3)
par(xpd=F)


inh %>%
  select(cluster, gene, rpm, n_cells) %>%
  pivot_wider(values_from=rpm, names_from=gene) %>%
  inner_join(filter(inh_residual_smry, tx=='aso6'), by='cluster') %>%
  select(x, cluster, n_cells,
         vip_rpm = Vip, sst_rpm = Sst, pvalb_rpm = Pvalb, npy_rpm = Npy,
         prnp_residual = residual, l95, u95) %>%
  arrange(x) -> inh_supp
write_supp_table(inh_supp, 'Cortical inhibitory neuron clusters: marker expression and residual Prnp.')


if (exists('rrs')) { rm(rrs) }
for (aggr in c('mouse_cortex_full','mouse_thalamus_full','mouse_cerebellum_full')) {
  
  this_tx = 'aso6'
  this_model = '3-aso6'
  
  ana = read_tsv(paste0('data/analytical/ana_',aggr,'.tsv'), col_types=cols())
  amt = read_csv(paste0('data/analytical/cla_',aggr,'.csv'), col_types=cols()) %>% clean_names()
  met = read_tsv(paste0('data/analytical/met_',aggr,'.tsv'), col_types=cols()) %>% clean_names()
  
  aco = read_tsv(paste0('data/analytical/aco_',aggr,'.tsv'), col_types=cols()) %>%
    filter(model %in% this_model) %>%
    mutate(weeks=as.numeric(gsub('-.*','',model)))
  
  ana %>%
    inner_join(amt, by='barcode') %>%
    rename(subtype=assignments) %>%
    filter(!grepl('exclude',subtype)) %>%
    group_by(subtype) %>%
    summarize(.groups='keep',
              rnaseh1_rpm = 1e6*mean(rnaseh1_umis/total_umis),
              umipc = mean(total_umis),
              prnp_basal = 1e6*mean(prnp_umis[tx=='pbs0']/total_umis[tx=='pbs0'])) %>%
    ungroup() -> stats
  
  stats %>%
    inner_join(aco, by='subtype') %>%
    filter(txc=='aso6') %>%
    select(subtype, n_cells, normed, rnaseh1_rpm, umipc, prnp_basal) %>%
    inner_join(met, by=c('subtype'='assignments')) -> rr_mashup
  
  rr_mashup$region = gsub('mouse_','',gsub('_full','',aggr))
  
  if (exists('rrs')) {
    rrs = rbind(rrs, rr_mashup)
  } else {
    rrs = rr_mashup
  }
  
}

par(mar=c(3,0.25,2,0.5))
xlims = range(rrs$rnaseh1_rpm) + c(-2, 2)
ylims = c(0, 1.00)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
mtext(side=4,line=0,adj=1,at=0:4/4, text=percent(0:4/4),cex=0.8,las=2)
mtext(side=4,line=-3,text=substitute(paste('',italic('Prnp'))), cex=0.8)

xlims = range(rrs$rnaseh1_rpm) + c(-2, 2)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=0:6*5, labels=NA, tck=-0.02)
axis(side=1, at=0:3*10, labels=NA, tck=-0.05)
axis(side=1, at=0:3*10, lwd=0, line=-0.5)
abline(h=1, lty=3)
axis(side=2, at=0:5/4, labels=NA, tck=-0.05)
# axis(side=2, at=0:5/4, labels=percent(0:5/4), las=2, lwd=0, line=-0.75)
points(x=rrs$rnaseh1_rpm, y=rrs$normed, cex=log10(rrs$n_cells)-2, col=rrs$color, pch=19)
mtext(side=1, line=1.6, text=substitute(paste(italic('Rnaseh1'), ' UPM')), cex=0.7)
mtext(LETTERS[panel], side=3, cex=1.5, adj = 0.0, line = 0.5)
panel = panel + 1


xlims = c(10^3, 10^5)
xats = rep(1:9,3) * 10^rep(3:5,each=9)
xbigs = 10^(3:5)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i', log='x')
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
axis(side=1, at=xats,  labels=NA, tck=-0.02)
axis(side=1, at=xbigs, labels=NA, tck=-0.05)
axis(side=1, at=xbigs, labels=c('1K','10K','100K'), lwd=0, line=-0.5)
abline(h=1, lty=3)
axis(side=2, at=0:5/4, labels=NA, tck=-0.05)
# axis(side=2, at=0:5/4, labels=percent(0:5/4), las=2, lwd=0, line=-0.75)
points(x=rrs$umipc, y=rrs$normed, cex=log10(rrs$n_cells)-2, col=rrs$color, pch=19)
mtext(side=1,  line=1.6, text=substitute(paste('UMIs / cell')), cex=0.7)
mtext(LETTERS[panel], side=3, cex=1.5, adj = 0.0, line = 0.5)
panel = panel + 1


xlims = c(30, 2000)
xats = rep(1:9,3) * 10^rep(1:3,each=9)
xbigs = 10^(1:3)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i', log='x')
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
axis(side=1, at=xats,  labels=NA, tck=-0.02)
axis(side=1, at=xbigs, labels=NA, tck=-0.05)
axis(side=1, at=xbigs, labels=c('10','100','1K'), lwd=0, line=-0.5)
abline(h=1, lty=3)
axis(side=2, at=0:5/4, labels=NA, tck=-0.05)
# axis(side=2, at=0:5/4, labels=percent(0:5/4), las=2, lwd=0, line=-0.75)
points(x=rrs$prnp_basal, y=rrs$normed, cex=log10(rrs$n_cells)-2, col=rrs$color, pch=19)
mtext(side=1,  line=1.6, text=substitute(paste('basal ',italic('Prnp'), ' UPM')), cex=0.7)
mtext(LETTERS[panel], side=3, cex=1.5, adj = 0.0, line = 0.5)
panel = panel + 1

rnaseh1_normed_wc = wtd.cor(rrs$rnaseh1_rpm, y=rrs$normed, weight = rrs$n_cells)
umipc_normed_wc = wtd.cor(rrs$umipc, y=rrs$normed, weight = rrs$n_cells)
prnp_normed_wc = wtd.cor(rrs$prnp_basal, y=rrs$normed, weight = rrs$n_cells)

write(paste('Rnaseh1 vs. residual weighted correlation: rho = ',formatC(rnaseh1_normed_wc['Y','correlation'], format='f', digits=3),
            ' P = ',formatC(rnaseh1_normed_wc['Y','p.value'], format='f', digits=3),'\n',sep=''),text_stats_path,append=T)

write(paste('UMI/cell vs. residual weighted correlation: rho = ',formatC(umipc_normed_wc['Y','correlation'], format='f', digits=3),
            ' P = ',formatC(umipc_normed_wc['Y','p.value'], format='f', digits=3),'\n',sep=''),text_stats_path,append=T)

write(paste('Basal Prnp vs. residual weighted correlation: rho = ',formatC(prnp_normed_wc['Y','correlation'], format='f', digits=3),
            ' P = ',formatC(prnp_normed_wc['Y','p.value'], format='f', digits=3),'\n',sep=''),text_stats_path,append=T)


unnecessary_message_end_of_figure_3 = dev.off()











##############
# FIGURE 4
##############

cat(file=stderr(), 'done.\nCreating Figure 4...')

resx=300
png(paste0('display_items/figure-4.png'),width=6.5*resx,height=2.0*resx,res=resx)

layout_matrix = matrix(c(1,1,1,1,1,1,2:11,12,13,13,13,14,14,14,15,15,15),nrow=1,byrow=F)
layout_widths = rep(1,dim(layout_matrix)[2])
layout_widths[layout_matrix[1,]==12] = 1.4
layout(layout_matrix, widths=layout_widths)

panel = 1

ana = read_tsv(    'data/analytical/ana_mouse_cortex_full.tsv', col_types=cols())
met = read_tsv('data/analytical/met_mouse_cortex_full.tsv', col_types=cols()) %>%
  mutate(order = row_number())
amt = read_csv('data/analytical/cla_mouse_cortex_full.csv', col_types=cols()) %>% 
  clean_names()

par(mar=c(2,2.5,3,1))
xlims = c(0, 1.25)
ylims = c(0, 1.25)
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F)
axis(side=1, at=0:5/4, labels=percent(0:5/4), line=-1.1, lwd=0, cex.axis=0.7)
axis(side=1, at=0:5/4, labels=NA, tck=-0.02)
mtext(side=1, line=1, text=substitute(paste('bulk qPCR ', italic('Prnp'))), cex=0.6)
axis(side=2, at=0:5/4, labels=percent(0:5/4), las=2, line=-0.75, lwd=0, cex.axis=0.7)
axis(side=2, at=0:5/4, labels=NA, tck=-0.02)
mtext(side=2, line=1.4, text=substitute(paste('single cell ', italic('Prnp'))), cex=0.6)
abline(h=1, lty=3)
abline(v=1, lty=3)
abline(a=0, b=1, lwd=0.5)

indiv_subs = sc_vs_qpcr_indiv %>% filter(weeks %in% c(2,12) & reg=='ss' & tx %in% c('pbs0','aso1','aso6'))
points(x=indiv_subs$prnp_qpcr, y=indiv_subs$sc_prnp_rel, pch=19, col=alpha(indiv_subs$color, ci_alpha))
smry_subs = sc_vs_qpcr_smry %>% filter(weeks %in% c(2,12) & reg=='ss' & tx %in% c('pbs0','aso1','aso6'))
segments(x0=smry_subs$prnp_qpcr_l95,  x1=smry_subs$prnp_qpcr_u95, y0=smry_subs$prnp_sc_mean, lwd=1.5, col=smry_subs$color)
segments(x0=smry_subs$prnp_qpcr_mean, y0=smry_subs$prnp_sc_l95,   y1=smry_subs$prnp_sc_u95,  lwd=1.5, col=smry_subs$color)
smry_subs %>%
  select(tx, weeks, prnp_sc_mean, prnp_qpcr_mean, color) %>%
  group_by(tx) %>%
  pivot_wider(names_from=weeks, values_from=c(prnp_sc_mean, prnp_qpcr_mean)) %>%
  filter(tx!='pbs0') -> bulk_sc_pivot
arrow_offset = 0.02
arrows(x0=bulk_sc_pivot$prnp_qpcr_mean_2 + arrow_offset, 
       y0=bulk_sc_pivot$prnp_sc_mean_2 + arrow_offset,
       x1=bulk_sc_pivot$prnp_qpcr_mean_12 - arrow_offset, 
       y1=bulk_sc_pivot$prnp_sc_mean_12 - arrow_offset,
       col=bulk_sc_pivot$color,
       code=2, angle=30, length=0.05, lwd=1.5)
label_offset = 0.03
smry_subs %>%
  filter(tx != 'pbs0') -> bulk_v_sc_active
text(x=bulk_v_sc_active$prnp_qpcr_mean, 
     y=bulk_v_sc_active$prnp_sc_mean + label_offset,
     labels=paste0(bulk_v_sc_active$weeks, ' weeks'),
     col = bulk_v_sc_active$color,
     pos=2,
     srt=-15, cex=0.45)
tx_meta %>%
  filter(tx %in% smry_subs$tx) -> leg
par(xpd=T)
legend(x=.75, y=1.6,
       leg$disp,
       col=leg$color,
       pch=19,
       lwd=2, 
       bty='n', cex=0.6)
par(xpd=F)
mtext(LETTERS[panel], side=3, cex=1.5, adj = 0.0, line = 0.4)
panel = panel + 1

aggr = 'mouse_cortex_full'

read_tsv(paste0('data/analytical/aco_',aggr,'.tsv'), col_types=cols()) %>%
  distinct(y, celltype, subtype) %>%
  arrange(desc(y)) -> yleg
yleg %>%
  group_by(celltype) %>%
  summarize(.groups='keep',
            miny=min(y),
            maxy=max(y),
            midy = (min(y)+max(y))/2) -> tranches

different_asos = c('aso6','aso1')
different_weeks = c(2, 12)


washout_models = paste0(rep(different_weeks,each=2),'-',rep(different_asos,2))

afi = read_tsv(paste0('data/analytical/afi_',aggr,'.tsv'), col_types=cols()) %>%
  filter(model %in% washout_models)
aco = read_tsv(paste0('data/analytical/aco_',aggr,'.tsv'), col_types=cols()) %>%
  filter(model %in% washout_models) %>%
  mutate(weeks=as.numeric(gsub('-.*','',model)))


# y axis legend
par(mar=c(6,0,3,0))
xlims = c(0, 1)
ylims = c(0, 1.5)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=4, line=0, tck=0.02, at=0:10/10,labels=NA)
axis(side=4, line=0, tck=0.05, at=0:3/2,labels=NA)
axis(side=4, line=-1.25, hadj=1, at=0:3/2,lwd=0, labels=percent(0:3/2), las=2, cex.axis=0.6)
mtext(side=4, line=-2.5, text=substitute(paste('residual ',italic('Prnp'))), cex=0.6)
mtext(LETTERS[panel], side=3, cex=1.5, adj = 0.0, line = 0.4)
panel = panel + 1

for (y in 1:nrow(yleg)) {
  this_subtype = yleg$subtype[y]
  acosub = aco %>% filter(subtype %in% this_subtype) %>% arrange(txc, weeks)
  afisub = afi %>% filter(subtype %in% this_subtype) %>% arrange(txc, weeks)
  
  par(mar=c(6,0,3,0))
  xlims = c(-1,14.5)
  ylims = c(0, 1.5)
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  abline(h=1, lty=3)
  axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
  axis(side=4, at=ylims, lwd.ticks=0, labels=NA)
  abline(h=c(0.5,1), lwd=0.25, lty=3)
  axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
  points(x=afisub$weeks, y=afisub$point_estimate, pch=20, col=alpha(afisub$color, ci_alpha))
  barwidth=1

  for (this_txc in unique(acosub$txc)) {
    acosubsub = acosub %>% filter(txc==this_txc)
    points(x=acosubsub$weeks, y=acosubsub$normed, col=acosubsub$color, lwd=2, type='l')
    polygon(x=c(acosubsub$weeks,rev(acosubsub$weeks)), y=c(acosubsub$l95,rev(acosubsub$u95)), col=alpha(acosubsub$color,ci_alpha), border=acosubsub$color, lwd=0.25, lty=3)
  }
  
  axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
  axis(side=1, at=different_weeks, labels=NA, tck=0.05)
  axis(side=1, at=different_weeks, line=-2, lwd=0, cex.axis=0.5)
  
  mtext(side=1, line=.25, las=2, text=this_subtype, cex=0.6)
}


aco %>%
  filter(txc!='_ctl') %>%
  select(txc, weeks, subtype, normed) %>%
  pivot_wider(values_from = normed, names_from=weeks) %>%
  mutate(recovery = `12`-`2`) %>%
  arrange(txc, desc(recovery)) %>%
  rename(residual_02wk = `2`, residual_12wk = `12`) -> recov_data

write_supp_table(recov_data, 'Washout of ASO 1 and ASO 6 between 2 and 12 weeks in mouse cortex by cell type.')

make_wa = function(aggr, reg, tx, weeks) {
  this_reg = reg
  this_tx = tx
  this_weeks = weeks
  
  amt = read_csv(paste0('data/analytical/cla_',aggr,'.csv'), col_types=cols()) %>% clean_names()
  ana = read_tsv(paste0('data/analytical/ana_',aggr,'.tsv'), col_types=cols()) %>%
    rename(wks=weeks, txx=tx) %>% # avoid name collision with function's arguments
    filter(wks %in% this_weeks) %>%
    filter(txx %in% c('pbs0',this_tx))
  amt %>%
    rename(subtype=assignments) %>%
    inner_join(ana, by='barcode') %>%
    group_by(subtype) %>%
    summarize(.groups='keep',
              total_cells = n(),
              rnaseh1_rpm = 1e6*mean(rnaseh1_umis/total_umis),
              prnp_rpm_basal = 1e6*mean(prnp_umis[txx=='pbs0']/total_umis[txx=='pbs0']),
              umipc = sum(total_umis)/n())  %>%
    ungroup() -> subtype_stats
  
  afi = read_tsv(paste0('data/analytical/afi_',aggr,'.tsv'), col_types=cols()) %>%
    filter(model %in% paste0(this_weeks,'-',this_tx))
  
  aco = read_tsv(paste0('data/analytical/aco_',aggr,'.tsv'), col_types=cols()) %>%
    filter(model %in% paste0(this_weeks,'-',this_tx))
  
  read_tsv(paste0('data/analytical/ana_',aggr,'.tsv'), col_types=cols()) %>% 
    group_by(sample, tx, weeks, reg) %>%
    summarize(.groups='keep',
              prnp_rpm = 1e6*mean(prnp_umis/total_umis)) %>%
    ungroup() -> step1
  
  step1 %>%
    filter(tx == 'pbs0') %>%
    group_by(weeks, reg) %>%
    summarize(.groups='keep',
              ctl_mean = mean(prnp_rpm)) %>%
    ungroup() -> ctl_means
  
  step1 %>%
    inner_join(ctl_means, by=c('weeks','reg')) %>%
    mutate(residual = prnp_rpm / ctl_mean) %>%
    mutate(subtype = '_bulk') %>%
    select(sample, weeks, tx, subtype, residual) -> bulk_data
  
  afi %>%
    filter(tx==this_tx) %>%
    filter(weeks %in% this_weeks) %>%
    inner_join(select(aco, subtype, txc, model, normed), by=c('subtype','txc','model')) %>%
    rename(residual = normed) %>%
    select(sample, color, weeks, tx, txc, subtype, n_cells, residual, prnp_umis, total_umis) -> sc_data
  
  sc_data %>%
    filter(weeks == min(weeks)) %>%
    arrange(sample, tx, weeks, subtype) %>%
    inner_join(select(bulk_data, sample, bulk_resid = residual), by='sample') %>%
    mutate(rres = residual - bulk_resid) %>%
    group_by(tx, subtype, color) %>%
    summarize(.groups='keep', 
              baseline_rres = mean(rres),
              baseline_resid = mean(residual)) %>%
    ungroup() -> rres
  
  sc_data %>% 
    mutate(timepoint = min_rank(weeks)) %>%
    mutate(timepoint = case_when(timepoint > 1 ~ 2, timepoint==1 ~ 1)) %>%
    group_by(tx, subtype, timepoint) %>%
    summarize(.groups='keep', mean_resid = mean(residual)) %>%
    ungroup() %>%
    pivot_wider(names_from = timepoint, values_from=mean_resid, names_prefix = 't') %>%
    mutate(delta_resid = t2 - t1) -> delta_res
  
  subtype_stats %>%
    filter(total_cells >= 100) %>%
    inner_join(rres, by=c('subtype')) %>%
    inner_join(delta_res, by=c('tx','subtype')) %>%
    select(tx, color, subtype, total_cells, baseline_resid, baseline_rres, delta_resid, umipc, prnp_rpm_basal, rnaseh1_rpm) -> wa
  
  return(wa)
}

wa = rbind(make_wa('mouse_cortex_full', 'ss', 'aso6', c(2,12)),
           make_wa('mouse_cortex_full', 'ss', 'aso1', c(2,12)))




# y axis legend
par(mar=c(2,0,3,0))
xlims = c(0, 1)
ylims = c(0, 1)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=4, line=0, tck=0.02, at=0:10/10,labels=NA)
axis(side=4, line=0, tck=0.05, at=0:2/2,labels=NA)
axis(side=4, line=-1.3, hadj=1, at=0:2/2,lwd=0, labels=percent(0:2/2), las=2, cex.axis=0.6)
mtext(side=4, line=-2.2, text='recovery from 2 to 12 weeks', cex=0.6)


par(mar=c(2,0,3,0.5))
xlims = c(100,1000)
ylims = c(0, 1)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i', log='x')
axis(side=1,at=1:10*100,labels=NA, tck=-0.02)
axis(side=1,at=c(1,10)*100,labels=NA, tck=-0.05)
axis(side=1,at=c(1,10)*100,labels=c('100','1K'),lwd=0,line=-1, cex.axis=0.6)
axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
axis(side=2, line=0, tck=-0.02, at=0:10/10,labels=NA)
axis(side=2, line=0, tck=-0.05, at=0:2/2,labels=NA)
mtext(side=1, line=0.75, cex=0.5, text=substitute(paste('basal ',italic('Prnp'),' UPM')))
points(wa$prnp_rpm_basal, wa$delta_resid, cex=(log10(wa$total_cells)-2.5), pch=19, col=wa$color)
wa %>% distinct(tx, color) %>% mutate(disp = gsub('ASO','ASO ',toupper(tx))) -> leg
mtext(LETTERS[panel], side=3, cex=1.5, adj = 0, line = 0.4)
panel = panel + 1


par(mar=c(2,0.5,3,0.5))
xlims = c(10,30)
ylims = c(0, 1)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1,at=10:30,labels=NA, tck=-0.02)
axis(side=1,at=1:3*10,labels=NA, tck=-0.05)
axis(side=1,at=1:3*10,lwd=0,line=-1, cex.axis=0.6)
mtext(side=1, line=0.75, cex=0.5, text=substitute(paste(italic('Rnaseh1'),' UPM')))
axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
axis(side=2, line=0, tck=-0.02, at=0:10/10,labels=NA)
axis(side=2, line=0, tck=-0.05, at=0:2/2,labels=NA)
points(wa$rnaseh1_rpm, wa$delta_resid, cex=(log10(wa$total_cells)-2.5), pch=19, col=wa$color)
wa %>% distinct(tx, color) %>% mutate(disp = gsub('ASO','ASO ',toupper(tx))) -> leg
mtext(LETTERS[panel], side=3, cex=1.5, adj = 0, line = 0.4)
panel = panel + 1

par(mar=c(2,1,3,0.5))
xlims = c(0, 1)
ylims = c(0, 1)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=0:10/10, labels=NA, tck=-0.02)
axis(side=1, line=-1, at=0:10/10, labels=percent(0:10/10), lwd=0, cex.axis=0.6)
mtext(side=1, line=0.75, cex=0.5, text=substitute(paste('2wk residual ',italic('Prnp'))))
axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
axis(side=2, line=0, tck=-0.02, at=0:10/10,labels=NA)
axis(side=2, line=0, tck=-0.05, at=0:2/2,labels=NA)
points(wa$baseline_resid, wa$delta_resid, cex=.5*(log10(wa$total_cells)-2), pch=19, col=wa$color)
abline(a=1,b=-1, lwd=0.5, lty=3)
abline(lm(delta_resid ~ baseline_resid, weights=total_cells, data=subset(wa, tx=='aso1')), col=wa$color[wa$tx=='aso1'][1])
abline(lm(delta_resid ~ baseline_resid, weights=total_cells, data=subset(wa, tx=='aso6')), col=wa$color[wa$tx=='aso6'][1])
wa %>% distinct(tx, color) %>% mutate(disp = gsub('ASO','ASO ',toupper(tx))) -> leg
mtext(LETTERS[panel], side=3, cex=1.5, adj = 0, line = 0.4)
panel = panel + 1

rnaseh1_wa_wc = wtd.cor(wa$rnaseh1_rpm, wa$delta_resid, weight = wa$total_cells)
prnp_wa_wc    = wtd.cor(wa$prnp_rpm_basal, wa$delta_resid, weight = wa$total_cells)
rres_wa_wc_aso1    = wtd.cor(wa$baseline_resid[wa$tx=='aso1'], wa$delta_resid[wa$tx=='aso1'], weight = wa$total_cells[wa$tx=='aso1'])
rres_wa_wc_aso6    = wtd.cor(wa$baseline_resid[wa$tx=='aso6'], wa$delta_resid[wa$tx=='aso6'], weight = wa$total_cells[wa$tx=='aso6'])

write(paste('Rnaseh1 vs. recovery weighted correlation: rho = ',formatC(rnaseh1_wa_wc['Y','correlation'], format='f', digits=3),
            ' P = ',formatC(rnaseh1_wa_wc['Y','p.value'], format='f', digits=3),'\n',sep=''),text_stats_path,append=T)

write(paste('Basal Prnp vs. recovery weighted correlation: rho = ',formatC(prnp_wa_wc['Y','correlation'], format='f', digits=3),
            ' P = ',formatC(prnp_wa_wc['Y','p.value'], format='f', digits=3),'\n',sep=''),text_stats_path,append=T)

write(paste('2-wk knockdown vs. recovery weighted correlation, ASO 1: rho = ',formatC(rres_wa_wc_aso1['Y','correlation'], format='f', digits=3),
            ' P = ',formatC(rres_wa_wc_aso1['Y','p.value'], format='f', digits=3),'\n',sep=''),text_stats_path,append=T)

write(paste('2-wk knockdown vs. recovery weighted correlation, ASO 6: rho = ',formatC(rres_wa_wc_aso6['Y','correlation'], format='f', digits=3),
            ' P = ',formatC(rres_wa_wc_aso6['Y','p.value'], format='f', digits=3),'\n',sep=''),text_stats_path,append=T)


unnecessary_message_end_of_figure_4 = dev.off()




















##############
# FIGURE 5
##############

cat(file=stderr(), 'done.\nCreating Figure 5...')

resx=300
png(paste0('display_items/figure-5.png'),width=6.5*resx,height=5*resx,res=resx)

layout_matrix = matrix(c(1,2,2,3,5,
                         1,2,2,3,5,
                         4,4,4,4,5,
                         4,4,4,4,5,
                         7,7,8,8,6,
                         7,7,8,8,6,
                         7,7,8,8,6,
                         7,7,8,8,6,
                         7,7,8,8,6), nrow=9, byrow=T)
layout(layout_matrix, widths=c(.1, .15, .05, .2, .5))
panel = 1



par(mar=c(2,0,2,0))
xlims = c(0, 1.25)
ylims = c(0, 1.25)
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F)
axis(side=4, at=0:5/4, labels=percent(0:5/4), las=2, line=-1.25, hadj=1, lwd=0, cex.axis=0.7)
axis(side=4, at=0:5/4, labels=NA, tck=0.02)
mtext(side=4, line=-2.75, text=substitute(paste('single cell ', italic('PRNP'))), cex=0.6)

for (aggr in c('cyno_cortex_full','cyno_cerebellum_full')) {
  
  this_reg = case_when(aggr=='cyno_cortex_full' ~ 'fc',
                       aggr=='cyno_cerebellum_full' ~ 'cb')
  this_reg_disp = gsub('cyno_','',gsub('_full','',aggr))
  
  par(mar=c(2,0,2,2))
  xlims = c(0, 1.25)
  ylims = c(0, 1.25)
  plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F)
  axis(side=1, at=0:5/4, labels=percent(0:5/4), line=-1.1, lwd=0, cex.axis=0.7)
  axis(side=1, at=0:5/4, labels=NA, tck=-0.02)
  mtext(side=1, line=1, text=substitute(paste('bulk qPCR ', italic('PRNP'))), cex=0.6)
  axis(side=2, at=0:5/4, labels=NA, tck=-0.02)
  abline(h=1, lty=3)
  abline(v=1, lty=3)
  abline(a=0, b=1, lwd=0.5, col='#A7A7A7')
  
  indiv_subs = sc_vs_qpcr_indiv %>% filter(weeks==13 & reg==this_reg & tx %in% c('acsf','aso7'))
  points(x=indiv_subs$prnp_qpcr, y=indiv_subs$sc_prnp_rel, pch=19, col=alpha(indiv_subs$color, ci_alpha))
  smry_subs = sc_vs_qpcr_smry %>% filter(weeks==13 & reg==this_reg & tx %in% c('acsf','aso7'))
  segments(x0=smry_subs$prnp_qpcr_l95,  x1=smry_subs$prnp_qpcr_u95, y0=smry_subs$prnp_sc_mean, lwd=1.5, col=smry_subs$color)
  segments(x0=smry_subs$prnp_qpcr_mean, y0=smry_subs$prnp_sc_l95,   y1=smry_subs$prnp_sc_u95,  lwd=1.5, col=smry_subs$color)
  
  mtext(side=3, line=0, text=this_reg_disp, cex=0.6)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.1, line = 0.1)
  panel = panel + 1
}


read_tsv('data/other/elisa_080_summary.tsv', col_types=cols()) %>% mutate(plate=80) -> elisa80
read_tsv('data/other/elisa_097_summary.tsv', col_types=cols()) %>% mutate(plate=97) -> elisa97
rbind(elisa80, elisa97) -> elisa

read_tsv('data/meta/cynos.tsv', col_types=cols()) -> cy_meta

elisa %>%
  filter(!grepl('QC',sample)) %>%
  mutate(animal=substr(sample,1,5),
         tissue=substr(sample,7,14)) %>%
  inner_join(cy_meta,by='animal') -> cynos

tissue_meta = tibble(tissue=c('cere ctx','frt ctx','csf'),
                     dispshort=c('Cb','FC','CSF'),
                     displong=c('cerebellum','cortex','CSF'),
                     tissue_x=c(2,1,3))

cynos %>%
  filter(disp=='aCSF') %>%
  group_by(tissue) %>%
  summarize(.groups='keep', 
            acsf_mean=mean(ngml_av)) -> acsf_means

tx_x_offset = 0.15
cynos %>%
  inner_join(tissue_meta,by='tissue') %>%
  inner_join(tx_meta, by=c('disp'='disp')) %>%
  inner_join(acsf_means, by='tissue') %>%
  mutate(x = tissue_x + ifelse(tx=='aso7',1,-1)*tx_x_offset,
         rel = ngml_av/acsf_mean) -> cyno_elisa

cyno_elisa %>%
  group_by(x, tissue, displong, tx, color) %>%
  summarize(.groups='keep',
            rel_mean = mean(rel),
            rel_l95 = mean(rel) - 1.96*sd(rel)/sqrt(n()),
            rel_u95 = mean(rel) + 1.96*sd(rel)/sqrt(n())) %>%
  ungroup() -> cyno_elisa_smry


par(mar=c(1.5,3,2,1))
xlims = c(0.25, 3.75)
ylims = c(0,1.5)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, lwd.ticks=0, labels=NA)
axis(side=2, at=0:8/2, labels=NA, tck=-0.02)
axis(side=2, at=0:8/2, labels=percent(0:8/2), lwd=0, line=-0.75, las=2, cex.axis=0.8)
mtext(side=2, line=2, cex=0.75, text='PrP protein')
abline(h=1, lty=3, lwd=0.5)
mtext(side=1, at=tissue_meta$tissue_x, text=tissue_meta$displong, line=0.1, cex=0.6)
points(cyno_elisa$x, cyno_elisa$rel, pch=19, col=alpha(cyno_elisa$color, ci_alpha))
barwidth=0.14
segments(x0=cyno_elisa_smry$x-barwidth, x1=cyno_elisa_smry$x+barwidth, y0=cyno_elisa_smry$rel_mean, lwd=1.5, col=cyno_elisa_smry$color)
arrows(x0=cyno_elisa_smry$x, y0=cyno_elisa_smry$rel_l95, y1=cyno_elisa_smry$rel_u95, code=3, angle=90, length=0.05, lwd=1.5, col=cyno_elisa_smry$color)
par(xpd=T)
tx_meta %>%
  filter(tx %in% c('aso7','acsf')) -> leg
legend(x=3,y=1.8,leg$disp,pch=16,col=leg$color,text.col=leg$color,bty='n',cex=0.7)
par(xpd=F)
mtext(LETTERS[panel], side=3, cex=1.5, adj = 0, line = 0.4)
panel = panel + 1

cyno_elisa %>%
  select(-animal, -tissue) %>%
  rename(animal=animal_id_short, 
         treatment=disp, 
         tissue=displong,
         prp_ngml = ngml_av,
         prp_relative = rel) %>%
  select(animal, sex, treatment, tissue, prp_ngml, prp_relative) -> cynos_out

write_supp_table(cynos_out, 'PrP ELISA results for cynomolgus macaques, indivudals.')


cyno_elisa_smry %>%
  inner_join(tx_meta, by='tx') %>%
  select(tissue=displong, treatment=disp, residual=rel_mean, l95=rel_l95, u95=rel_u95) -> cynos_smry_out

write_supp_table(cynos_smry_out, 'PrP ELISA results for cynomolgus macaques, summarized.')

aggr = 'cyno_cortex_full'
read_tsv(paste0('data/analytical/aco_',aggr,'.tsv'), col_types=cols()) %>%
  distinct(y, celltype, subtype) -> yleg
yleg %>%
  group_by(celltype) %>%
  summarize(.groups='keep',
            miny=min(y),
            maxy=max(y),
            midy = (min(y)+max(y))/2) -> tranches

for (aggr in c('cyno_cortex_full','cyno_cerebellum_full')) {
  
  this_weeks = 13
  this_tx = 'aso7'
  this_model = '13-aso7'
  
  afi = read_tsv(paste0('data/analytical/afi_',aggr,'.tsv'), col_types=cols()) %>%
    filter(model==this_model)
  aco = read_tsv(paste0('data/analytical/aco_',aggr,'.tsv'), col_types=cols()) %>%
    filter(model==this_model) 
  
  aco %>%
    distinct(y, celltype, subtype) -> yleg
  yleg %>%
    group_by(celltype) %>%
    summarize(.groups='keep',
              miny=min(y),
              maxy=max(y),
              midy = (min(y)+max(y))/2) -> tranches
  
  if (aggr=='cyno_cortex_full') {
    par(mar=c(0.5,9,2,1))
  } else if (aggr=='cyno_cerebellum_full') {
    par(mar=c(3,9,0,1))
  }
  
  xlims = c(0, 1.5)
  ylims = range(afi$y) + c(-0.5, 0.5)
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=0:8/4, labels=NA, tck=-0.02)
  abline(v=1, lwd=0.25, lty=3)
  abline(v=0.5, lwd=0.125, lty=3)
  axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
  mtext(side=2, at=yleg$y, text=yleg$subtype, line=0.25, las=2, cex=0.7)
  tranche_line = 7.5
  overhang = 0.4
  for (i in 1:nrow(tranches)) {
    axis(side=2, line=tranche_line, at=c(tranches$miny[i], tranches$maxy[i]) + c(-1,1)*overhang, tck=0.02, labels=NA)
    mtext(side=2, line=tranche_line + 0.25, at=tranches$midy[i], text=tranches$celltype[i], cex=0.8)
  }
  points(y=afi$y, x=afi$point_estimate, pch=19, col=alpha(afi$color, ci_alpha))
  barwidth=0.4
  segments(y0=aco$y-barwidth, y1=aco$y+barwidth, x0=aco$normed, col=aco$color, lwd=2)
  arrows(y0=aco$y, x0=aco$l95, x1=aco$u95, angle=90, length=0.03, code=3, col=aco$color, lwd=2)
  
  reg = gsub('cyno_','',gsub('_full','',aggr))
  mtext(side=2, line=tranche_line+1.25, text=reg, font=1, cex=1)
  
  if (aggr=='cyno_cortex_full') {
    mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.1, line = 0.2)
    panel = panel + 1
  } else if (aggr=='cyno_cerebellum_full') {
    axis(side=1, at=0:8/4, labels=percent(0:8/4), lwd=0, line=-0.75, cex=0.9)
    mtext(side=1, line=1.6, text=substitute(paste('residual ',italic('PRNP'))), cex=0.7)
  }
}


par(mar=c(2.5,3,2,1))

this_model = '13-aso7'
aggr = 'cyno_cortex_full'
met = read_tsv(paste0('data/analytical/met_',aggr,'.tsv'), col_types = cols())
aco_cy = read_tsv(paste0('data/analytical/aco_',aggr,'.tsv'), col_types=cols()) %>%
  filter(model==this_model) %>%
  filter(txc=='aso7')

this_model = '2-aso6'
aggr = 'mouse_cortex_full'

aco_mo = read_tsv(paste0('data/analytical/aco_',aggr,'.tsv'), col_types=cols()) %>%
  filter(model==this_model) %>%
  filter(txc=='aso6') %>%
  mutate(subtype = case_when(subtype %in% c('pericyte','fibroblast') ~ 'pericyte/fibroblast',
                             TRUE ~ subtype)) %>%
  group_by(subtype) %>%
  summarize(.groups='keep', 
            normed = sum(normed*total_cells)/sum(total_cells),
            total_cells = sum(total_cells)) %>%
  ungroup()

aco_cy %>%
  inner_join(aco_mo, by='subtype', suffix=c('_cy','_mo')) %>%
  mutate(total_cells = total_cells_cy + total_cells_mo) %>%
  inner_join(met, by=c('subtype'='assignments'), suffix=c('_aco','_met')) %>%
  select(subtype, normed_cy, normed_mo, total_cells, color=color_met) -> cortex_cy_mo

cortex_offsets = tibble(subtype=c('astrocyte','inhibitory','OPC','oligodendrocyte','microglia','endothelial stalk','pericyte/fibroblast','excitatory'),
                        x_offset = c(0,0,0,0,0,0,0,0),
                        y_offset = c(0,0,0,0,0,0,0,-0.02),
                        pos = c(4,1,4,3,1,4,4,2))

cortex_cy_mo %>% 
  left_join(cortex_offsets, by='subtype') %>% 
  mutate(x_offset = replace_na(x_offset, 0)) %>%
  mutate(y_offset = replace_na(y_offset, 0)) %>%
  mutate(pos = replace_na(pos, 2)) -> cortex_cy_mo

plot(NA, NA, xlim=c(0,1), ylim=c(0,1), axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=0:10/10, labels=NA, tck=-0.02)
axis(side=1, at=0:2/2, labels=NA,tck=-0.05)
axis(side=1, at=0:2/2, labels=percent(0:2/2), line=-0.75, lwd=0)
mtext(side=1, line=1.2, text='cynomolgus residual', cex=0.6)
axis(side=2, at=0:10/10, labels=NA, tck=-0.02)
axis(side=2, at=0:2/2, labels=NA,tck=-0.05)
axis(side=2, at=0:2/2, labels=percent(0:2/2), las=2, line=-0.5, lwd=0)
mtext(side=2, line=2.1, text='mouse residual', cex=0.6)
abline(a=0, b=1, lwd=0.5, col='#707070')
points(cortex_cy_mo$normed_mo, cortex_cy_mo$normed_cy, cex=log10(cortex_cy_mo$total_cells)-2, pch=19, col=cortex_cy_mo$color)
#text(cortex_cy_mo$normed_mo, cortex_cy_mo$normed_cy, labels=paste0(cortex_cy_mo$subtype, '  '), adj=1, srt=35, cex=0.6)
text(cortex_cy_mo$normed_mo + cortex_cy_mo$x_offset, cortex_cy_mo$normed_cy + cortex_cy_mo$y_offset, labels=cortex_cy_mo$subtype, pos=cortex_cy_mo$pos, cex=0.6)
mtext(side=3, line=0.0, cex=0.6, text='cortex')
mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.1, line = 0.2)
panel = panel + 1

cortex_cy_mo_wc = wtd.cor(cortex_cy_mo$normed_mo, cortex_cy_mo$normed_cy, weight = cortex_cy_mo$total_cells)

write(paste("Cortex cyno vs. mouse weighted Pearson's correlation: rho = ",formatC(cortex_cy_mo_wc['Y','correlation'], format='f', digits=3),
            ' P =  ',formatC(cortex_cy_mo_wc['Y','p.value'], format='e', digits=2),'\n',sep=''),text_stats_path,append=T)


this_model = '13-aso7'
aggr = 'cyno_cerebellum_full'
met = read_tsv(paste0('data/analytical/met_',aggr,'.tsv'), col_types = cols())
aco_cy = read_tsv(paste0('data/analytical/aco_',aggr,'.tsv'), col_types=cols()) %>%
  filter(model==this_model) %>%
  filter(txc=='aso7') 

this_model = '3-aso6'
aggr = 'mouse_cerebellum_full'
aco_mo = read_tsv(paste0('data/analytical/aco_',aggr,'.tsv'), col_types=cols()) %>%
  filter(model==this_model) %>%
  filter(txc=='aso6') %>%
  mutate(subtype = case_when(subtype %in% c('pericyte','fibroblast') ~ 'pericyte/fibroblast',
                             TRUE ~ subtype)) %>%
  group_by(subtype) %>%
  summarize(.groups='keep', 
            normed = sum(normed*total_cells)/sum(total_cells),
            total_cells = sum(total_cells)) %>%
  ungroup()

aco_cy %>%
  inner_join(aco_mo, by='subtype', suffix=c('_cy','_mo')) %>%
  mutate(total_cells = total_cells_cy + total_cells_mo) %>%
  inner_join(met, by=c('subtype'='assignments'), suffix=c('_aco','_met')) %>%
  select(subtype, normed_cy, normed_mo, total_cells, color=color_met) -> cerebellum_cy_mo


cerebellum_offsets = tibble(subtype=c('UBC','OPC','oligodendrocyte','astrocyte','endothelial stalk','microglia','Bergmann glia','pericyte/fibroblast','Golgi','granule'),
                            y_offset = c(0.00, 0.00, 0.02, 0.01,0.00,0.02,0.01,0.00,0.00,0.00),
                            x_offset = c(-0.03,0.02,-0.02,0,0,0,0,0,0.02,0),
                            pos=c(4,2,4,1,1,1,3,3,2,4))

cerebellum_cy_mo %>% 
  left_join(cerebellum_offsets, by='subtype') %>% 
  mutate(y_offset = replace_na(y_offset, 0)) %>%
  mutate(x_offset = replace_na(x_offset, 0)) %>%
  mutate(pos=replace_na(pos, 2)) -> cerebellum_cy_mo

plot(NA, NA, xlim=c(0,1), ylim=c(0,1), axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=0:10/10, labels=NA, tck=-0.02)
axis(side=1, at=0:2/2, labels=NA,tck=-0.05)
axis(side=1, at=0:2/2, labels=percent(0:2/2), line=-0.75, lwd=0)
mtext(side=1, line=1.2, text='cynomolgus residual', cex=0.6)
axis(side=2, at=0:10/10, labels=NA, tck=-0.02)
axis(side=2, at=0:2/2, labels=NA,tck=-0.05)
axis(side=2, at=0:2/2, labels=percent(0:2/2), las=2, line=-0.5, lwd=0)
mtext(side=2, line=2.1, text='mouse residual', cex=0.6)
abline(a=0, b=1, lwd=0.5, col='#707070')
points(cerebellum_cy_mo$normed_mo, cerebellum_cy_mo$normed_cy, cex=log10(cerebellum_cy_mo$total_cells)-2, pch=19, col=cerebellum_cy_mo$color)
text(x=cerebellum_cy_mo$normed_mo + cerebellum_cy_mo$x_offset, 
     y=cerebellum_cy_mo$normed_cy + cerebellum_cy_mo$y_offset, 
     labels=cerebellum_cy_mo$subtype, 
     pos=cerebellum_cy_mo$pos, 
     cex=0.6, col='#000000')
mtext(side=3, line=0.0, cex=0.6, text='cerebellum')
mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.1, line = 0.2)
panel = panel + 1

cerebellum_cy_mo_wc = wtd.cor(cerebellum_cy_mo$normed_mo, cerebellum_cy_mo$normed_cy, weight = cerebellum_cy_mo$total_cells)

write(paste("Cerebellum cyno vs. mouse weighted Pearson's correlation: rho = ",formatC(cerebellum_cy_mo_wc['Y','correlation'], format='f', digits=3),
            ' P =  ',formatC(cerebellum_cy_mo_wc['Y','p.value'], format='e', digits=2),'\n',sep=''),text_stats_path,append=T)



unnecessary_message_end_of_figure_5 = dev.off()





##############
# FIGURE 6
##############

cat(file=stderr(), 'done.\nCreating Figure 6...')

resx=300
png(paste0('display_items/figure-6.png'),width=6.5*resx,height=7*resx,res=resx)

layout_matrix = matrix(c(1,1,
                         2,2,
                         3,4,
                         5,5), nrow=4, byrow=T)
layout(layout_matrix, heights=c(0.8,0.7,1,1.25))
panel = 1



mets = rbind_files('data/analytical/','met_')
subtype_meta = mets %>%
  rename(subtype=assignments) %>%
  distinct(subtype, color) %>%
  group_by(subtype) %>%
  slice(1)


acos %>%
  filter(txc != '_ctl') %>%
  mutate(region_sorter = case_when(region == 'cortex' ~ 1,
                                   region == 'thalamus' ~ 2,
                                   region == 'cerebellum' ~ 3)) %>%
  mutate(weeks = as.integer(gsub('-.*','',model))) %>%
  arrange(region_sorter, desc(species), txc, weeks) %>%
  distinct(region_sorter, region, species, txc, weeks) %>%
  mutate(x = row_number()) -> adr_meta

acos %>%
  select(-x, -y) %>%
  filter(txc != '_ctl') %>%
  mutate(weeks = as.integer(gsub('-.*','',model))) %>%
  rename(tx_color = color) %>%
  group_by(region, weeks, txc, model) %>%
  ungroup() %>%
  inner_join(sc_vs_qpcr_smry, by=c('aggr','region','weeks','txc'='tx')) %>%
  select(-color) %>%
  mutate(target_bulk_residual = case_when(txc == 'asom' ~ malat1_sc_mean,
                                          txc %in% c('aso1','aso6','aso7') ~ prnp_sc_mean)) %>%
  mutate(difference = normed - target_bulk_residual,
         difference_l95 = l95 - target_bulk_residual,
         difference_u95 = u95 - target_bulk_residual) %>%
  inner_join(adr_meta, by=c('species','region','txc','weeks')) %>%
  select(region, model, txc, weeks, x, subtype, normed, difference, difference_l95, difference_u95, total_cells) %>%
  inner_join(subtype_meta, by='subtype') %>%
  arrange(x, desc(total_cells)) -> adr

adr %>%
  group_by(region) %>%
  summarize(.groups='keep', minx=min(x), maxx=max(x), midx=mean(x)) -> x_tranches

adr %>%
  distinct(x, region, txc, weeks) %>%
  inner_join(tx_meta, by=c('txc'='txc')) %>%
  mutate(xlab = paste0(disp,'\n',weeks,' wk')) -> xlabs



par(mar=c(0,4,2,1))
xlims = range(adr$x) + c(-0.5, 0.5)
ylims = c(-0.5, 0.6)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
axis(side=2, at=-10:10/10, labels=NA, tck=-0.02)
axis(side=2, at=-10:10/10, labels=percent(-10:10/10, signed=T), line=-0.75, lwd=0, las=2, cex.axis=0.6)
mtext(side=2, line=2.5, text='difference from overall residual', cex=0.6)
mtext(side=2, line=1.6, at=c(-.33, .33), text=c('\u2190deeper KD','weaker KD\u2192'), cex=0.55)
abline(h=0, lty=3)
points(x=adr$x, y=adr$difference, col=adr$color, pch=19, cex=pmax(log10(adr$total_cells)-2,0.05))
mtext(LETTERS[panel], side=3, cex=1.5, adj = 0.0, line = 0.2)
panel = panel + 1

adr %>%
  group_by(x, region, txc, weeks) %>%
  summarize(.groups='keep',
            wsd = sqrt(wtd.var(x=normed, weights=total_cells))) %>%
  ungroup() -> adr_sd

adr_sd %>%
  group_by(region) %>%
  summarize(.groups='keep', 
            mwsd = mean(wsd)) -> adr_sd_smry

write(paste("Mean weighted standard deviations by tissue: ",paste(adr_sd_smry$region, percent(adr_sd_smry$mwsd), sep=': ', collapse=', '),'\n',sep=''),text_stats_path,append=T)

par(mar=c(5,4,1,1))
xlims = range(adr$x) + c(-0.5, 0.5)
ylims = c(0,0.17)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
axis(side=2, at=0:100/100, labels=NA, tck=-0.02)
axis(side=2, at=0:20/20, labels=percent(0:20/20), line=-0.75, lwd=0, las=2, cex.axis=0.6)
mtext(side=2, line=1.2, text='wtd SD',cex=0.6)
abline(h=0, lty=3)
barwidth = 0.4
rect(xleft=adr_sd$x-barwidth, xright=adr_sd$x+barwidth, ybottom=rep(0, nrow(adr_sd)), ytop=adr_sd$wsd, col='#7C7C7C', border=NA)
par(xpd=T)
x_offset = 0.1
y_offset = -0.02
text(x=xlabs$x + x_offset, y=rep(min(ylims),nrow(xlabs)) + y_offset, srt=60, adj=1, labels=xlabs$xlab, cex=0.8)
tranche_line = 3.6
overhang_left = 0.4
overhang_right = 0.4
for (i in 1:nrow(x_tranches)) {
  axis(side=1, line=tranche_line, at=c(x_tranches$minx[i], x_tranches$maxx[i]) + c(-1,1)*c(overhang_left, overhang_right), tck=0.05, labels=NA)
  mtext(side=1, line=tranche_line, at=x_tranches$midx[i]-overhang_left/2+overhang_right/2, text=x_tranches$region[i], cex=0.6)
}
par(xpd=F)

adr %>%
  group_by(region, subtype) %>%
  summarize(.groups='keep',
            mean_difference = mean(difference),
            mean_absolute_diference = mean(abs(difference))) -> adr_by_type

write_supp_table(adr_by_type, 'Mean differences from bulk tissue residual target by cell type.')

acos %>%
  filter(txc != '_ctl') %>%
  mutate(subtype = case_when(subtype %in% c('pericyte','fibroblast') ~ 'pericyte/fibroblast',
                             TRUE ~ subtype)) %>%
  filter(region %in% c('cerebellum','cortex')) %>%
  group_by(species, region, txc, model, celltype, subtype) %>%
  summarize(.groups='keep',
            normed = sum(normed*n_cells)/sum(n_cells),
            n_cells = sum(n_cells)) %>% # note order is important here, because otherwise this n_cells assignment would affect previous line
  ungroup() %>%
  mutate(cgram_id = paste(species,region,model,sep='-')) %>%
  select(cgram_id, species, region, txc, model, celltype, subtype, normed, n_cells) -> cgram_data

cgram_data %>%
  distinct(cgram_id, region, species, model) %>%
  mutate(weeks = as.numeric(gsub('-.*','',model)),
         tx = gsub('.+-','',model)) %>%
  arrange(region, desc(species), desc(tx), weeks) %>%
  group_by(region) %>%
  mutate(y = max(row_number()) - row_number() + 1,
         x = row_number()) %>%
         #x = max(row_number()) - row_number() + 1) %>%
  ungroup()  -> cgram_meta

cgram_meta %>%
  select(region, cgram_id, x) %>%
  inner_join(select(cgram_meta, region, cgram_id, y), by=c('region'), suffix=c('_left','_right')) %>%
  group_by(region) %>%
  filter(y + x <= max(y+x)/2+1) %>%
  ungroup() %>%
  mutate(row = row_number()) %>%
  add_column(rho = as.numeric(NA), p=as.numeric(NA)) -> cgram

for (i in 1:nrow(cgram)) {
  if (cgram$cgram_id_left[i] == cgram$cgram_id_right[i]) {
    cgram$rho[i] = 1
    cgram$p[i] = 0
    next
  }
  left = cgram_data %>% filter(cgram_id == cgram$cgram_id_left[i])
  right = cgram_data %>% filter(cgram_id == cgram$cgram_id_right[i])
  left %>%
    inner_join(right, by='subtype', suffix=c('_left','_right')) %>%
    mutate(total_cells = n_cells_left + n_cells_right) %>%
    select(subtype, total_cells, normed_left, normed_right) -> this_data
  cor_obj = wtd.cor(x=this_data$normed_left, y=this_data$normed_right, weight = this_data$total_cells)
  cgram$rho[i] = cor_obj['Y','correlation']
  cgram$p[i] = cor_obj['Y','p.value']
}

# from https://github.com/ericminikel/prp_mrm/blob/master/src/mrm_figures.R
color_ramp = c("#230007","#67000d", "#a50f15", "#cb181d", "#ef3b2c", "#fb6a4a", "#fc9272", "#fcbba1", "#fee0d2", "#fff5f0",'#FFFFFF','#f7fcf0','#e0f3db','#ccebc5','#a8ddb5','#7bccc4','#4eb3d3','#2b8cbe','#0868ac','#084081','#000000')
corparams = data.frame(floor=(-10:10)/10, color=color_ramp)

cgram$floor = floor(cgram$rho*10)/10
cgram$floor[cgram$rho==1] = 1
cgram$color = corparams$color[match(cgram$floor, corparams$floor)]
cgram$border_lty = ifelse(cgram$p < 0.05, 1, 3)
cgram$border_lwd = ifelse(cgram$p < 0.05, 1.5, 0.5)

cgram %>%
  mutate(cgram_id_left = gsub('aso7','ason',cgram_id_left)) %>%
  mutate(cgram_id_right = gsub('aso7','ason',cgram_id_right)) -> cgram_out

write_supp_table(cgram_out, 'Correlation matrix of residual target by cell type in all conditions.')


par(mar=c(4,8,1,1))
this_region = 'cortex'
subs = cgram %>% filter(region==this_region)
xlims = range(cgram$x) + c(-0.5, 0.5)
ylims = range(cgram$y) + c(-0.5, 0.5)
radius = 0.4
plot(x=NA, y=NA, xlim = xlims, ylim=ylims, axes=F, ann=F)
rect(xleft=subs$x-radius, xright=subs$x+radius, ybottom=subs$y-radius, ytop=subs$y+radius,
     col=subs$color, border=NA)
rect(xleft=subs$x-radius, xright=subs$x+radius, ybottom=subs$y-radius, ytop=subs$y+radius,
     col=NA, lty=subs$border_lty, lwd=subs$border_lwd, border='#000000')
leg = cgram_meta %>% 
  filter(region==this_region) %>%
  inner_join(tx_meta, by='tx') %>%
  mutate(longdisp = paste(species, disp, '-', weeks, 'weeks'))
par(xpd=T)
text(x=rep(0.4, nrow(leg)), y=leg$y, labels=leg$longdisp, adj=1, cex=0.7)
text(x=leg$x, y=rep(0.4, nrow(leg)), labels=leg$longdisp, srt=30, adj=1, cex=0.7)
par(xpd=F)
mtext(side=3, at=3.5, line=-1.2, text=this_region, cex=0.7)
mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.1, line = -0.4)
panel = panel + 1

par(mar=c(4,8,1,1))
this_region = 'cerebellum'
subs = cgram %>% filter(region==this_region)
xlims = range(cgram$x) + c(-0.5, 0.5)
ylims = range(cgram$y) + c(-0.5, 0.5)
radius = 0.4
plot(x=NA, y=NA, xlim = xlims, ylim=ylims, axes=F, ann=F)
rect(xleft=subs$x-radius, xright=subs$x+radius, ybottom=subs$y-radius, ytop=subs$y+radius,
     col=subs$color, border=NA)
rect(xleft=subs$x-radius, xright=subs$x+radius, ybottom=subs$y-radius, ytop=subs$y+radius,
     col=NA, lty=subs$border_lty, lwd=subs$border_lwd, border='#000000')
leg = cgram_meta %>% 
  filter(region==this_region) %>%
  inner_join(tx_meta, by='tx') %>%
  mutate(longdisp = paste(species, disp, '-', weeks, 'weeks'))
par(xpd=T)
text(x=rep(0.4, nrow(leg)), y=leg$y, labels=leg$longdisp, adj=1, cex=0.7)
text(x=leg$x, y=rep(0.4, nrow(leg)), labels=leg$longdisp, srt=30, adj=1, cex=0.7)
par(xpd=F)
mtext(side=3, at=2.5, line=-1.2, text=this_region, cex=0.7)
mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.1, line = -0.4)
panel = panel + 1

scale_xleft = 5.0
scale_xright = 6.6
scale_factor = (scale_xright - scale_xleft) / 21
scale_ybot = 5.1
scale_ytop = 5.4
y_offset = 0.1
par(xpd=T)
rect(xleft=scale_xleft,xright=scale_xright,ybottom=scale_ybot,ytop=scale_ytop,lwd=2.5)
rect(xleft=scale_xleft+(0:20)*scale_factor, xright=scale_xleft+(1:21)*scale_factor, ybottom=rep(scale_ybot,5), ytop=rep(scale_ytop,5), border=NA, col=corparams$color)
text(x=(scale_xright + scale_xleft)/2, y=scale_ytop-y_offset, pos=3, labels="correlation coefficient", cex=0.6)
text(x=c(scale_xleft,mean(c(scale_xleft,scale_xright)),scale_xright),y=c(scale_ybot,scale_ybot,scale_ybot)+y_offset,pos=1,labels=c('-1','0','1'),cex=0.6,font=3)
# p value scale
pleg_xleft = 5.0
pleg_ybot = 2.5
pleg_side = 0.5
pleg_space = 0.6
text(x=pleg_xleft + pleg_side + pleg_space/2, y=pleg_ybot+pleg_side-y_offset, pos=3, labels="P value", cex=0.6)
rect(xleft=pleg_xleft + c(0,1)*(pleg_space+pleg_side),
     xright=pleg_xleft + c(0,1)*(pleg_space+pleg_side) + pleg_side,
     ybottom=pleg_ybot + c(0,0),
     ytop=pleg_ybot + c(0,0) + pleg_side,
     col=NA,
     lty=c(3,1),
     lwd=c(0.5,1.5))
text(x=pleg_xleft + c(0,1)*(pleg_space+pleg_side) + pleg_side/2,
     y=pleg_ybot + c(0,0) + y_offset,
     pos=1,
     labels=c('≥0.05','<0.05'),
     cex=0.6, font=3)
par(xpd=F)




all_ana = rbind_files('data/analytical/','ana_.*full') %>%
  mutate(aggr = gsub('ana_','',file)) %>%
  select(-file) %>%
  relocate(aggr)
all_ana %>%
  group_by(aggr, tx, weeks, animal_short) %>%
  summarize(.groups='keep',
            prnp_rpm = 1e6*mean(prnp_umis/total_umis),
         malat1_rpm = 1e6*mean(malat1_umis/total_umis)) %>%
  ungroup() %>%
  group_by(aggr, weeks) %>%
  mutate(prnp_rel = prnp_rpm / mean(prnp_rpm[tx %in% c('acsf','pbs0')]),
         malat1_rel = malat1_rpm / mean(malat1_rpm[tx %in% c('acsf','pbs0')])) %>%
  ungroup() %>%
  mutate(target_rel = case_when(tx=='asom' ~ malat1_rel,
                                TRUE ~ prnp_rel)) -> sc_indiv_all
sc_indiv_all %>%  
  filter(!(tx %in% c('acsf','pbs0'))) %>%
  mutate(model = paste0(weeks,'-',tx)) %>%
  group_by(aggr, model) %>%
  summarize(.groups='keep',
            residual = mean(target_rel)) %>%
  ungroup() -> sc_overall_all

  
acos %>%
  filter(txc != '_ctl') %>%
  rename(tx=txc) %>%
  filter(celltype=='neuron') %>%
  mutate(weeks=as.numeric(gsub('-.*','',model))) %>%
  select(aggr, species, region, weeks, tx, model, subtype, normed, l95, u95, color) %>%
  inner_join(sc_overall_all, by=c('aggr','model')) %>%
  rename(subtype_residual=normed, overall_residual=residual) %>%
  mutate(subtype_sorter = case_when(subtype %in% c('excitatory','inhibitory') ~ paste0('_',subtype),
                                    TRUE ~ subtype)) %>%
  arrange(subtype_sorter, region, desc(species), desc(tx), weeks) %>%
  mutate(x = row_number(),
         y = max(row_number()) - row_number() + 1) %>%
  mutate(differential = subtype_residual - overall_residual,
         dl = l95 - overall_residual,
         du = u95 - overall_residual) -> nvo

nvo %>%
  rename(difference = differential, difference_l95 = dl, difference_u95 = du) %>%
  inner_join(select(tx_meta, tx, disp), by='tx') %>%
  select(aggr, species, region, weeks, treatment=disp, subtype, subtype_residual, overall_residual, difference, difference_l95, difference_u95) -> nvo_out
write_supp_table(nvo_out, 'Difference of subtype residual and overall residual for all neurons in all conditions.')


nvo %>%
  group_by(subtype) %>%
  summarize(.groups='keep',
            miny = min(y), maxy = max(y), midy = mean(y),
            minx = min(x), maxx = max(x), midx = mean(x)) -> s_leg

nvo %>%
  group_by(subtype, region) %>%
  summarize(.groups='keep',
            miny = min(y), maxy = max(y), midy = mean(y),
            minx = min(x), maxx = max(x), midx = mean(x)) %>%
  ungroup() %>%
  mutate(rshort = case_when(region=='cerebellum' ~ 'cb',
                            region=='cortex' ~ 'cx',
                            region=='thalamus' ~ 'th')) -> sr_leg

nvo %>%
  group_by(subtype, region, weeks) %>%
  summarize(.groups='keep',
            miny = min(y), maxy = max(y), midy = mean(y),
            minx = min(x), maxx = max(x), midx = mean(x)) -> srw_leg

par(mar=c(3.5,4,2,1))
xlims = range(nvo$x) + c(-0.5, 0.5)
ylims = c(-.5, .5) # range(c(nvo$dl, nvo$du))  + c(-1,1)*0.05
#plot(x=nvo$x, y=nvo$subtype_residual - nvo$overall_residual)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
axis(side=2, at=-10:10/10, labels=NA, tck=-0.02)
axis(side=2, at=-10:10/10, labels=percent(-10:10/10, signed=T), line=-0.75, lwd=0, las=2, cex.axis=0.6)
abline(h=0, lty=3)
barwidth=0.4
segments(x0=nvo$x-barwidth, x1=nvo$x+barwidth, y0=nvo$differential, col=nvo$color, lwd=2)
arrows(x0=nvo$x, y0=nvo$dl, y1=nvo$du, col=nvo$color, code=3, angle=90, length=0.05, lwd=2)
x_offset = -0.05
mtext(side=1, at=nvo$x, text=paste0(nvo$weeks,'w'), line=-0.2, cex=0.45)
sr_line = 1
overhang = 0.4
for (i in 1:nrow(sr_leg)) {
  axis(side=1, line=sr_line, at=c(sr_leg$minx[i], sr_leg$maxx[i]) + c(-1,1)*overhang, tck=0.02, labels=NA)
  mtext(side=1, line=sr_line, at=sr_leg$midx[i], text=sr_leg$rshort[i], cex=0.6)
}
s_line = 2.2
overhang = 0.4
for (i in 1:nrow(s_leg)) {
  axis(side=1, line=s_line, at=c(s_leg$minx[i], s_leg$maxx[i]) + c(-1,1)*overhang, tck=0.02, labels=NA)
  mtext(side=1, line=s_line, at=s_leg$midx[i], text=s_leg$subtype[i], cex=0.7)
}
mtext(side=2, line=2.5, text='difference from overall residual', cex=0.6)
mtext(side=2, line=1.6, at=c(-.33, .33), text=c('\u2190deeper KD','weaker KD\u2192'), cex=0.55)
tx_meta %>% filter(txc!='_ctl') -> active_tx
par(xpd=T)
legend(x=min(xlims),y=max(ylims)+0.05, active_tx$disp, col=active_tx$color, pch=15, bty='n', cex=0.8)
par(xpd=F)
mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.05, line = 0.4)
panel = panel + 1

unnecessary_message_end_of_figure_6 = dev.off()




#####
# FIGURE S1
#####

cat(file=stderr(), 'done.\nCreating Figure S1...')

resx=300
png(paste0('display_items/figure-s1.png'),width=6.5*resx,height=7*resx,res=resx)

layout_matrix = matrix(c(1,6,
                         2,7,
                         3,8,
                         4,9,
                         5,10), nrow=5, byrow=T)
layout(layout_matrix)
panel = 1


weeks_meta = tibble(weeks=c(2,3,12,13),
                    color=c('#12FE00','#01FEFE','#FE34FE','#F5785A'))

tx_meta %>%
  mutate(act = case_when(txc=='_ctl' ~ 'inactive', TRUE ~ 'active')) %>%
  mutate(acol = case_when(txc=='_ctl' ~ '#A9A9A9', TRUE ~ '#FFB00F')) -> act_meta

act_meta %>%
  distinct(act, acol) -> act_leg

aggrs = c("mouse_cerebellum_full", 
          "mouse_cortex_full", 
          "mouse_thalamus_full",
          "cyno_cerebellum_full", 
          "cyno_cortex_full")

par(mar=c(2,2,2,5))

for (aggr in aggrs) {
  aggr_disp = gsub('cyno','cynomolgus',gsub('_',' ',gsub('_full','',aggr)))
  this_aggr = aggr
  amt = read_csv(paste0('data/analytical/cla_',this_aggr,'.csv'), col_types=cols()) %>% clean_names()
  all_ana %>% 
    filter(aggr==this_aggr) %>%
    inner_join(amt, by='barcode') %>%
    filter(!grepl('exclude',assignments)) -> ana
  ana$act_color = act_meta$acol[match(ana$tx, act_meta$tx)]
  set.seed(1)
  ana %>%
    mutate(random_order = runif(nrow(ana),0,1)) %>%
    arrange(random_order) -> shuf
  xlims = range(shuf$umap_1) * 1.1
  ylims = range(shuf$umap_2) * 1.1
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
  mtext(side=1, line=0.25, text='UMAP 1', cex=0.8)
  axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
  mtext(side=2, line=0.25, text='UMAP 2', cex=0.8)
  points(x=shuf$umap_1, y=shuf$umap_2, pch=20, cex=0.05, col=alpha(shuf$act_color,ci_alpha))
  mtext(side=3, line=0, text=aggr_disp)
  par(xpd=T)
  legend(x=max(xlims),y=max(ylims),legend=act_leg$act,col=act_leg$acol,title = 'treatment group',title.col='#000000',pch=15,bty='n',cex=0.8)
  par(xpd=F)
}

for (aggr in aggrs) {
  aggr_disp = gsub('_',' ',gsub('_full','',aggr))
  this_aggr = aggr
  all_ana %>% filter(aggr==this_aggr) -> ana
  ana$weeks_color = weeks_meta$color[match(ana$weeks, weeks_meta$weeks)]
  set.seed(1)
  ana %>%
    mutate(random_order = runif(nrow(ana),0,1)) %>%
    arrange(random_order) -> shuf
  xlims = range(shuf$umap_1) * 1.1
  ylims = range(shuf$umap_2) * 1.1
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
  mtext(side=1, line=0.25, text='UMAP 1', cex=0.8)
  axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
  mtext(side=2, line=0.25, text='UMAP 2', cex=0.8)
  points(x=shuf$umap_1, y=shuf$umap_2, pch=20, cex=0.05, col=alpha(shuf$weeks_color,ci_alpha))
  mtext(side=3, line=0, text=aggr_disp)
  weeks_meta %>% filter(weeks %in% ana$weeks) -> weeks_leg
  par(xpd=T)
  legend(x=max(xlims),y=max(ylims),legend=weeks_leg$weeks,col=weeks_leg$color,title = 'weeks',title.col='#000000',pch=15,bty='n',cex=0.8)
  par(xpd=F)
}

unnecessary_message_end_of_figure_s1 = dev.off()




#####
# SUPPLEMENT
#####

cat(file=stderr(), 'done.\nFinalizing supplementary tables...')

# write the supplement directory / table of contents
supplement_directory %>% rename(table_number = name, description=title) -> contents
addWorksheet(supplement,'contents')
bold_style = createStyle(textDecoration = "Bold")
writeData(supplement,'contents',contents,headerStyle=bold_style,withFilter=T)
freezePane(supplement,'contents',firstRow=T)
# move directory to the front
original_order = worksheetOrder(supplement)
n_sheets = length(original_order)
new_order = c(n_sheets, 1:(n_sheets-1))
worksheetOrder(supplement) = new_order
activeSheet(supplement) = 'contents'
# now save
saveWorkbook(supplement,supplement_path,overwrite = TRUE)


#######
# DATA FOR SINGLE CELL PORTAL DEPOSITION
#######

cat(file=stderr(), 'done.\nGenerating files for Single Cell Portal deposition...')

all_amt = rbind_files('data/analytical/','cla_') %>%
  mutate(aggr=gsub('cla_','',file))

species_onto = tibble(species_short=c('mouse','cyno'),
                      species_id=c('NCBITaxon_10090','NCBITaxon_9541'),
                      species__ontology_label=c('Mus musculus','Macaca fascicularis'))

organ_onto = tibble(region_long=c('cerebellum','cortex','thalamus'),
                    uberon_id=c('UBERON_0002037','UBERON_0000956','UBERON_0004703'),
                    uberon_label=c('cerebellum','cerebral cortex','dorsal thalamus'))

cell_onto = read_tsv('data/other/cell_ontology_mapping.tsv', col_types=cols())

all_ana %>%
  inner_join(select(all_amt, barcode, assignments, aggr), by=c('barcode','aggr')) %>%
  mutate(NAME=paste0(aggr,'_',barcode)) %>%
  rename(biosample_id=sample) %>%
  rename(donor_id=animal_short) %>%
  inner_join(species_onto, by=c('species'='species_short')) %>%
  select(-species) %>%
  rename(species=species_id) %>% 
  mutate(disease='PATO_0000461', disease__ontology_label='normal') %>%
  inner_join(organ_onto, by=c('region_long')) %>%
  rename(organ=uberon_id, organ__ontology_label=uberon_label) %>% 
  mutate(library_preparation_protocol='EFO_0009922', library_preparation_protocol__ontology_label="10x 3' v3") %>% 
  mutate(sex=case_when(sex=='m' ~ 'male', sex=='f' ~ 'female')) %>% 
  inner_join(tx_meta, by='tx') %>%
  rename(treatment=disp) %>%
  rename(weeks_post_dose=weeks) %>%
  inner_join(cell_onto, by='assignments') %>%
  mutate(cell_type_assigned = case_when(grepl('exclude',assignments) ~ 'exclude', TRUE ~ assignments)) %>%
  select(NAME, 
         biosample_id, 
         donor_id, 
         species, 
         species__ontology_label,
         disease,
         disease__ontology_label,
         organ,
         organ__ontology_label,
         library_preparation_protocol,
         library_preparation_protocol__ontology_label,
         sex,
         cell_type,
         cell_type__ontology_label,
         cell_type_assigned,
         total_umis,
         malat1_umis,
         prnp_umis,
         rnaseh1_umis,
         treatment,
         weeks_post_dose,
         aggr,
         umap_1,
         umap_2) -> scp_metadata

scp_metadata %>% 
  summarize_all(class) %>%
  mutate(across(NAME:umap_2, ~ gsub('character','group',gsub('integer','numeric',.)))) %>%
  mutate(NAME='TYPE') -> double_header
  
write_tsv(double_header, 'display_items/scp_metadata.tsv')
write_tsv(scp_metadata, 'display_items/scp_metadata.tsv', col_names=F, append = T)

for (this_aggr in unique(scp_metadata$aggr)) {
  
  scp_metadata %>%
    filter(aggr == this_aggr) %>%
    select(NAME, 
           X=umap_1, 
           Y=umap_2, 
           aggr,
           Category=cell_type_assigned) -> scp_clustering
  
  scp_clustering %>% 
    summarize_all(class) %>%
    mutate(across(NAME:Category, ~ gsub('character','group',gsub('integer','numeric',.)))) %>%
    mutate(NAME='TYPE') -> double_header
  
  clust_out_path = paste0('display_items/scp_clustering_',gsub('_full','',this_aggr),'.tsv')
  write_tsv(double_header, clust_out_path)
  write_tsv(scp_clustering, clust_out_path, col_names=F, append = T)
  
}




elapsed_time = Sys.time() - overall_start_time
cat(file=stderr(), paste0('done.\nAll tasks complete in ',round(as.numeric(elapsed_time),1),' ',units(elapsed_time),'.\n'))






