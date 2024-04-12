# alterations done to the current freeze after download
# these only exist to avoid re-running the whole thing and keeping a record
# changes required are implemented upstream simultaneously. these changes are applied
# to light data.holder directly to keep the original copy intact

library(dplyr)
library(magrittr)
library(gemma.R)
library(pbapply)
devtools::load_all()
dh = load_dif_exp()

# exposing annotation sources
data.holder.light = readRDS(file.path(DATADIR, 'DATA.HOLDER.light.rds'))


re_processAnnotations = function(x){
    attr = attributes(x)
    out = attr$env$response %>% gemma.R:::processAnnotations()
    attributes(out) = c(attributes(out),attr[!names(attr) %in% names(attributes(out))])
    return(out)
}


data.holder.light$human$contrast_metadata$experiment.annotations %<>% pblapply(re_processAnnotations)
data.holder.light$mouse$contrast_metadata$experiment.annotations %<>% pblapply(re_processAnnotations)
data.holder.light$rat$contrast_metadata$experiment.annotations %<>% pblapply(re_processAnnotations)



get_n_e = function(x,mask){
    pbsapply(seq_len(ncol(x)),function(i){
        x[,i][mask] %>% {!is.na(.)} %>% sum
    })
}

data.holder.light$human$gene_metadata$n.E = get_n_e(dh$human$data$adj.pv, dh$human$contrast_metadata$geeq.rawData ==1)
data.holder.light$mouse$gene_metadata$n.E = get_n_e(dh$mouse$data$adj.pv, dh$mouse$contrast_metadata$geeq.rawData ==1)
data.holder.light$rat$gene_metadata$n.E = get_n_e(dh$rat$data$adj.pv, dh$rat$contrast_metadata$geeq.rawData ==1)


data.holder.light$human$gene_metadata$DE.prior = data.holder.light$human$gene_metadata$n.DE/data.holder.light$human$gene_metadata$n.E
data.holder.light$mouse$gene_metadata$DE.prior = data.holder.light$mouse$gene_metadata$n.DE/data.holder.light$mouse$gene_metadata$n.E
data.holder.light$rat$gene_metadata$DE.prior = data.holder.light$rat$gene_metadata$n.DE/data.holder.light$rat$gene_metadata$n.E


saveRDS(data.holder.light,file.path(DATADIR, 'DATA.HOLDER.light.rds'))

