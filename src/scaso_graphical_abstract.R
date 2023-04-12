options(stringsAsFactors=F)
if(interactive()) setwd('~/d/sci/src/scaso')
suppressMessages(library(tidyverse))
suppressMessages(library(janitor))
suppressMessages(library(weights))
suppressMessages(library(openxlsx))
suppressMessages(library(Hmisc))
suppressMessages(library(MASS)); select = dplyr::select; summarize = dplyr::summarize

alpha = function(rgb_hexcolor, proportion) {
  hex_proportion = sprintf("%02x",round(proportion*255))
  rgba = paste(rgb_hexcolor,hex_proportion,sep='')
  return (rgba)
}

ci_alpha = 0.35 # degree of transparency for shading confidence intervals in plot

resx=300
png('display_items/graphical_abstract_umap.png',width=5*resx,height=4.5*resx,res=resx)

met = read_tsv('data/analytical/met_mouse_cerebellum_full.tsv', col_types=cols())
ana = read_tsv('data/analytical/ana_mouse_cerebellum_full.tsv', col_types=cols())
amt = read_csv('data/analytical/cla_mouse_cerebellum_full.csv', col_types=cols()) %>% clean_names()
aco = read_tsv('data/analytical/aco_mouse_cerebellum_full.tsv', col_types=cols())

ana %>%
  inner_join(amt, by='barcode') %>%
  filter(!grepl('^excl',assignments)) %>%
  inner_join(met, by=c('assignments')) %>%
  select(umap_1, umap_2, color) -> umap

par(mar=c(1,1,1,1))
layout_matrix = matrix(c(1:7, rep(15,7), 8:14), nrow=3, byrow=T)
layout(layout_matrix, heights=c(1,2.5,1))

for(this_subtype in met$assignments) {
  aco %>%
    select(-color) %>%
    filter(subtype==this_subtype & model=='12-asom' & txc=='asom') %>%
    inner_join(met, by=c('subtype'='assignments')) -> this_aco
  xwidth = 0.5
  xcenter = 0
  xleft = xcenter - xwidth
  xright = xcenter + xwidth
  height = 1
  ybottom = 0
  ytop = ybottom + height
  yfill = ybottom + height * this_aco$normed
  plot(NA, NA, xlim=c(xleft, xright)*1.1, ylim=c(ybottom, ytop)*1.1, axes=F, ann=F, xaxs='i', yaxs='i')
  rect(xleft=xleft, xright=xright, ybottom=ybottom, ytop=ytop, border=this_aco$color, col=NA, lwd=2)
  rect(xleft=xleft, xright=xright, ybottom=ybottom, ytop=yfill, border=NA, col=this_aco$color, lwd=2)
}
xlims = range(umap$umap_1) * 1.1 
ylims = range(umap$umap_2) * 1.1
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
points(umap$umap_1, umap$umap_2, col=alpha(umap$color, 0.2), pch=20, cex=0.5)
dev.off()



n_nuclei = 50

resx=300
png('display_items/graphical_abstract_sc1.png',width=1*resx,height=1*resx,res=resx)
par(mfrow=c(1,1), mar=c(0,0,0,0))
set.seed(1)
sc_diagram = tibble(color = sample(met$color, replace=T, size=n_nuclei),
                    x = runif(n=n_nuclei, min=0, max=1),
                    y = runif(n=n_nuclei, min=0, max=1))
plot(sc_diagram$x, sc_diagram$y, col=sc_diagram$color, pch=20, axes=F, ann=F)
dev.off()

resx=300
png('display_items/graphical_abstract_sc2.png',width=1*resx,height=1*resx,res=resx)
par(mfrow=c(1,1), mar=c(0,0,0,0))
set.seed(2)
sc_diagram = tibble(color = sample(met$color, replace=T, size=n_nuclei),
                    x = runif(n=n_nuclei, min=0, max=1),
                    y = runif(n=n_nuclei, min=0, max=1))
plot(sc_diagram$x, sc_diagram$y, col=sc_diagram$color, pch=20, axes=F, ann=F)
dev.off()







resx=300
png('display_items/graphical_abstract_schist.png',width=5*resx,height=2*resx,res=resx)

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
binsize = 10
this_subtype = 'granule'
model_data %>%
  filter(tx %in% c('asom','pbs0')) %>%
  filter(weeks==12) %>%
  filter(this_subtype=='all' | subtype==this_subtype) -> this_subtype_data
asom_hist = hist(this_subtype_data$malat1_umis[this_subtype_data$tx=='asom'], breaks=seq(0,1e5,binsize), plot=F)
pbs0_hist = hist(this_subtype_data$malat1_umis[this_subtype_data$tx=='pbs0'], breaks=seq(0,1e5,binsize), plot=F)


par(mar=c(1,1,1,1))
ctl_color = '#A7A7A7'
mal_color = met$color[met$assignments=='granule']
xmax = as.numeric(ceiling(quantile(this_subtype_data$malat1_umis[this_subtype_data$tx=='pbs0'], .9)/binsize)*binsize)
ymax = max(asom_hist$counts)
xlims = c(0,xmax*1.2)
xbigs = seq(0,xmax,binsize*5)
xats = seq(0,xmax,binsize)
ylims = c(0,ymax)
plot(NA, NA, xlim=xlims, ylim = ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
rect(xleft=pbs0_hist$mids-binsize/2, xright=pbs0_hist$mids+binsize/2, ybottom=rep(0,length(pbs0_hist)), ytop=pbs0_hist$counts, col=alpha(ctl_color, .8), border=NA)
rect(xleft=asom_hist$mids-binsize/2, xright=asom_hist$mids+binsize/2, ybottom=rep(0,length(asom_hist)), ytop=asom_hist$counts, col=alpha(mal_color, .8), border=NA)
dev.off()