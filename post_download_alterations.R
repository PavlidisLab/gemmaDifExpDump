# alterations done to the current freeze after download
# these only exist to avoid re-running the whole thing and keeping a record
# changes required are implemented upstream simultaneously. these changes are applied
# to light data.holder directly to keep the original copy intact

library(dplyr)
library(magrittr)
library(gemma.R)
library(pbapply)


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

saveRDS(data.holder.light,file.path(DATADIR, 'DATA.HOLDER.light.rds'))

