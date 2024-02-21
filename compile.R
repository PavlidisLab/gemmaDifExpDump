print('data compilation')
devtools::load_all()
library(gemma.R)
library(magrittr)
library(dplyr)
library(pbapply)
library(data.table)

# options ---------

# paths. defaults are defined in package, can be overwritten
# DATADIR <- '/cosmos/data/project-data/GemmaDifExp'
# FREEZEDIR <- '/cosmos/data/project-data/GemmaDifExp/gemma_freeze'
# RAWDIR <- "/cosmos/data/project-data/GemmaDifExp/raw_data"
# CACHEDIR <- '/space/scratch/Ogdata/GemmaDifExp/gemma_cache'

# backups are saved in the middle of the process in case
# we wish to update the data partially
load_from_backups <- TRUE

# if FALSE, cache will be deleted
trust_cache <- TRUE
use_cache <- TRUE

# set up cache --------
gemma_memoise(use_cache,CACHEDIR)

if(!trust_cache){
    # just in case...
    stop("You are trying to delete the cache. It's currently disabled for your own good")
    # gemma.R::forget_gemma_memoised()
}


# objects to create for each species
lapply(c('human','mouse','rat'),function(taxon){
    dir.create(file.path(RAWDIR,'midway_backups',taxon),showWarnings =FALSE,recursive = TRUE)
    print(taxon)

    print('Getting all datasets')

    bc_path <- file.path(RAWDIR,'midway_backups',taxon,'all_datasets.rds')
    if(load_from_backups && file.exists(bc_path)){
        all_datasets <- readRDS(bc_path)
    } else {
        all_datasets <- gemma.R::get_datasets(taxa = taxon) %>%
            gemma.R::get_all_pages(file = bc_path,overwrite = TRUE)
    }

    print('Getting platforms')
    bc_path = file.path(RAWDIR,'midway_backups',taxon,'platform_ids.rds')
    if(load_from_backups && file.exists(bc_path)){
        all_platform_ids <- readRDS(bc_path)
    } else {
        all_platforms <- get_platforms_by_ids(taxa = taxon) %>%
            get_all_pages()
        all_platforms %<>% filter(platform.experimentCount != 0)
        all_platform_ids <- all_platforms$platform.ID
        saveRDS(all_platform_ids, bc_path)
    }

    print('Getting platform annotations')
    bc_path = file.path(RAWDIR,'midway_backups',taxon,'all_platform_annotations.rds')
    if(load_from_backups && file.exists(bc_path)){
        all_platform_annotations <- readRDS(bc_path)
    } else {
        all_platform_annotations <- all_platform_ids %>% pblapply(function(x){
            gemma.R::get_platform_annotations(x)
        })
        names(all_platform_annotations) <- all_platform_ids
        saveRDS(all_platform_annotations,bc_path)
    }

    null_platforms <- all_platform_annotations %>% purrr::map_lgl(is.null)
    all_platform_annotations <- all_platform_annotations[!null_platforms]


    print('Compiling the datasets table')
    bc_path = file.path(RAWDIR,'midway_backups',taxon,'metadata_first_pass.rds')

    if(load_from_backups && file.exists(bc_path)){
        contrast_metaData <- readRDS(bc_path)
    } else {
        seq_len(nrow(all_datasets)) %>% pblapply(function(i){
            dataset <- all_datasets[i,]

            dataset_annotations <- gemma.R::get_dataset_annotations(dataset$experiment.ID)


            dataset$experiment.annotations = list(dataset_annotations)



            differential <- gemma.R::get_dataset_differential_expression_analyses(dataset$experiment.ID)

            q_types <- gemma.R::get_dataset_quantitation_types(dataset$experiment.ID)

            if(length(differential)==0){
                return(NULL)
            }

            samples <- tryCatch(gemma.R::get_dataset_samples(dataset$experiment.ID),
                                error = function(e){
                                    if(e$message == "500: Internal server error."){
                                        return(list())
                                    }else{
                                        stop(e$message)
                                    }
                                })


            # resultset based sample counts
            resultSet_sampleCount <- seq_len(nrow(differential)) %>% sapply(function(i){
                resultSet_sampleCount <- gemma.R:::subset_factorValues(samples$sample.factorValues,
                                                                       differential_expressions = differential,
                                                                       resultSet = differential$result.ID[i]) %>% sum
            })

            contrast_sampleCount <- seq_len(nrow(differential)) %>% sapply(function(i){
                resultSet_sampleCount <- gemma.R:::subset_factorValues(samples$sample.factorValues,
                                                                       differential_expressions = differential,
                                                                       resultSet = differential$result.ID[i],
                                                                       contrast = differential$contrast.ID[i]) %>% sum
            })



            platform <- gemma.R::get_dataset_platforms(dataset$experiment.ID)
            platform_annots <- platform$platform.ID %>% lapply(function(x){
                all_platform_annotations[[as.character(x)]]
            }) %>% do.call(rbind,.)



            out <- data.table(result_contrast.ID = paste0(differential$result.ID,'.',differential$contrast.ID),
                              dataset,
                              experiment.scale = q_types %>% dplyr::filter(type == 'processed') %$% scale,
                              platform.numGenes = platform_annots$NCBIids %>% {.[.!='' | grepl('|',.,fixed = TRUE)]} %>% unique %>% length, # this counts genes a bit differently than what gemma displays as I think it should be closer to the intended purpose. gemma counts by splitting genes that are aligned to the same probeset which presumably shouldn't be used as differential expression results. need to reconsider
                              platform.ID = platform$platform.ID %>% paste0(collapse = ','),
                              differential %>% dplyr::select(-experiment.ID),
                              resultSet.sampleCount = resultSet_sampleCount,
                              contrast.sampleCount = contrast_sampleCount)


        })  %>% do.call(rbind,.) -> contrast_metaData

        saveRDS(contrast_metaData, bc_path)
    }


    differentials <- contrast_metaData$result.ID %>% unique
    print("Getting differential expression data")
    bc_path =  file.path(RAWDIR,'midway_backups',taxon,'differentials.rds')
    bc_path2 = file.path(RAWDIR,'midway_backups',taxon,'differentials_filtered.rds')

    if(load_from_backups && file.exists(bc_path) && !file.exists(bc_path2)){
        all_differential_values <- readRDS(bc_path)
    } else if(!(file.exists(bc_path2) && load_from_backups)){
        all_differential_values <- differentials %>% pblapply(function(x){
            tryCatch(
                gemma.R::get_differential_expression_values(resultSet = x)[[1]],
                error = function(e){
                    if(grepl('(502)|(404)',e$message)){
                        return(data.frame())
                    } else{
                        stop(e$message)
                    }
                }
            )

        })
        names(all_differential_values)= differentials
        saveRDS(all_differential_values, bc_path)
    }


    print('Processing differential expression data')
    if(load_from_backups && file.exists(bc_path2)){
        all_differential_values <- readRDS(bc_path2)
    } else {
        bad_differentials = all_differential_values %>% sapply(nrow) %>% {.==0}
        all_differential_values = all_differential_values[!bad_differentials]


        # remove duplicates, calculate additional statistics
        # looking at nathaniel's old freeze, it appears that
        # probeset with the lowest p value was selected to
        # to represent a gene. need a deeper look to make sure
        # but assuming that is the case for now -Ogan
        all_differential_values %<>% pblapply(function(x){
            x %<>% dplyr::filter(NCBIid!='' & !grepl('|',NCBIid,fixed = TRUE))
            p_value_cols = names(x)[grepl('contrast.*?pvalue',names(x))]
            logfc_cols = names(x)[grepl('contrast.*?log2fc',names(x))]
            stat_cols =names(x)[grepl('contrast.*?tstat',names(x))]
            dup_genes = x$NCBIid[duplicated(x$NCBIid)]
            out = x %>% dplyr::filter(!NCBIid %in% dup_genes) %>% dplyr::select(-Probe,-pvalue,-corrected_pvalue,-rank)
            # add dups with lowest p value
            if (length(p_value_cols)==0){
                return(data.frame())
            }

            lapply(seq_along(p_value_cols),function(i){
                x %>% dplyr::filter(NCBIid %in% dup_genes) %>% dplyr::arrange(!!sym(p_value_cols[i])) %>%
                    dplyr::filter(!duplicated(NCBIid)) %>%
                    dplyr::select(NCBIid,GeneSymbol,GeneName,!!sym(logfc_cols[i]),!!sym(stat_cols[i]),!!sym(p_value_cols[i]))
            }) %>% {
                if(length(.)>1){
                    out = .[[1]]
                    for (j in seq_len(length(.)-1)){
                        out = merge(out,.[[j+1]])
                    }
                    out
                } else if (length(.) == 1){
                    .[[1]]
                }
            } -> to_append

            out = rbind(out,to_append)

            # calculate additional stats
            for(i in seq_along(p_value_cols)){
                out[[paste0(p_value_cols[[i]],'_adjusted')]] = p.adjust(out[[p_value_cols[i]]],'BH')
            }
            return(out)
        })
        bad_differentials = all_differential_values %>% sapply(nrow) %>% {.==0}
        all_differential_values = all_differential_values[!bad_differentials]
        saveRDS(all_differential_values, bc_path2)
    }


    # construct the dataHolder object-----------
    dh = list()

    ncbi_ids = all_differential_values %>% lapply(function(x){
        as.integer(x$NCBIid)
    }) %>% unlist %>% unique %>% sort

    print('calculating adj.p values')
    bc_path = file.path(RAWDIR,'midway_backups',taxon,'adj.pv.rds')
    if(load_from_backups && file.exists(bc_path)){
        adj.pv = readRDS(bc_path)
    } else{
        pblapply(seq_len(nrow(contrast_metaData)),function(i){
            result_id = contrast_metaData$result.ID[i]
            contrast_id = contrast_metaData$contrast.ID[i]
            differential = all_differential_values[[as.character(result_id)]]
            differential[[
                glue::glue('contrast_{contrast_id}_pvalue_adjusted')
            ]][match(ncbi_ids,differential$NCBIid)]
        }) %>% {names(.) = contrast_metaData$result_contrast.ID;.} %>% do.call(cbind,.) -> adj.pv
        rownames(adj.pv) = ncbi_ids
        full_na_genes =  adj.pv %>% apply(1,function(x){all(is.na(x))})
        full_na_contrasts = adj.pv %>% apply(2,function(x){all(is.na(x))})
        adj.pv = adj.pv[!full_na_genes,!full_na_contrasts]

        saveRDS(adj.pv, file.path(bc_path))
    }

    # remove any filtered gene from the ncbi_ids
    ncbi_ids = rownames(adj.pv)
    dh$adj.pv = adj.pv


    contrast_metaData %<>% dplyr::filter(result_contrast.ID %in% colnames(dh$adj.pv))


    dh$fc = lapply(seq_len(nrow(contrast_metaData)),function(i){
        result_id = contrast_metaData$result.ID[i]
        contrast_id = contrast_metaData$contrast.ID[i]
        differential = all_differential_values[[as.character(result_id)]]
        differential[[
            glue::glue('contrast_{contrast_id}_log2fc')
        ]][match(ncbi_ids,differential$NCBIid)]
    }) %>% {names(.) = contrast_metaData$result_contrast.ID;.}%>% do.call(cbind,.) %>%
        {rownames(.) = ncbi_ids;.}


    assertthat::assert_that(all(contrast_metaData$result_contrast.ID == colnames(adj.pv)))

    contrast_metaData$contrast.difExpCount = matrixStats::colSums2(adj.pv <= 0.05, na.rm = T)
    contrast_metaData$contrast.meanFoldChange = matrixStats::colMeans2(abs(dh$fc), na.rm = T)
    contrast_metaData$contrast.meanUpFoldChange = matrixStats::colMeans2(dh$fc * ifelse(dh$fc > 0, 1, NA), na.rm = T)
    contrast_metaData$contrast.meanDownFoldChange = matrixStats::colMeans2(dh$fc * ifelse(dh$fc < 0, 1, NA), na.rm = T)

    saveRDS(contrast_metaData, file.path(RAWDIR,'midway_backups',taxon,'metadata_final.rds'))



    print('compiling gene metadata')


    seq(1,length(ncbi_ids),50) %>% pblapply(function(i){
        gemma.R::get_genes(ncbi_ids[seq(i,min(i+49,length(ncbi_ids)))],raw= FALSE, memoised = FALSE)
    }) %>% do.call(rbind,.) -> gene_metaData


    gene_metaData = gene_metaData[match(ncbi_ids,gene_metaData$gene.NCBI),]
    gene_metaData$n.DE = matrixStats::rowSums2(dh$adj.pv[, contrast_metaData$geeq.rawData == 1] <= 0.05, na.rm = T)
    gene_metaData$dist.mean =  rowMeans(dh$fc[, contrast_metaData$geeq.rawData == 1], na.rm = T)
    gene_metaData$dist.SD = Rfast::rowVars(dh$fc[,  contrast_metaData$geeq.rawData == 1], na.rm = T, std = T)

    # identify genes from the wrong species that snuck inside the metadata
    bad_genes = !gene_metaData$taxon.name %in% taxon
    # assertthat::assert_that(all(!bad_genes))
    assertthat::assert_that(all(gene_metaData$entrez.ID == rownames(dh$adj.pv)))
    assertthat::assert_that(all(gene_metaData$entrez.ID == rownames(dh$fc)))

    # identify contrasts that includes those genes
    bad_contrasts = dh$adj.pv[bad_genes,] %>% apply(2,function(x){!all(is.na(x))})

    # remove bad genes and bad contrasts from the results
    gene_metaData = gene_metaData[!bad_genes,]
    dh$adj.pv = dh$adj.pv[!bad_genes,]
    dh$fc = dh$fc[!bad_genes,]
    dh$adj.pv = dh$adj.pv[,!bad_contrasts]
    dh$fc = dh$fc[,!bad_contrasts]
    contrast_metaData = contrast_metaData[!bad_contrasts,]
    assertthat::assert_that(all(contrast_metaData$rsc.ID == colnames(dh$adj.pv)))


    print('calculating z scores')
    dh$zscore <- (dh$fc - gene_metaData$dist.mean) / gene_metaData$dist.SD

    print('getting go terms')
    go_terms <- mygene::queryMany(gene_metaData$gene.NCBI, scopes = 'entrezgene', fields = 'go', species = taxon)

    go_terms <- data.table::rbindlist(lapply(c('CC', 'BP', 'MF'), function(cat) {
        gocat <- paste0('go.', cat)
        data.table::rbindlist(lapply(1:nrow(go_terms), function(indx) {
            row <- go_terms@listData[[gocat]][[indx]]
            if(!is.null(row))
                data.frame(entrez.ID = go_terms@listData$query[indx], category = cat, id = row$id, term = row$term)
        }), fill = T)
    }), fill = T)

    contrast_metaData$nonNa.numGenes = dh$fc %>% apply(2,\(x){
        sum(!is.na(x))
    })

    out = list(
        taxon = taxon,
        data = dh,
        contrast_metadata = contrast_metaData,
        gene_metadata = gene_metaData,
        go = unique(go_terms)
    )

    dir.create(file.path(RAWDIR,taxon),showWarnings = FALSE)
    saveRDS(out, file.path(RAWDIR,taxon,'full_data.rds'))
    return(out)
}) -> data.holder

names(data.holder) = c('human','mouse','rat')
saveRDS(data.holder, file.path(RAWDIR, 'DATA.HOLDER.rds'))


fbm_path = file.path(DATADIR,'data_fbm')


for(taxon in names(data.holder)) {
    create_fbm(data.holder[[taxon]]$data$fc,
               file.path(fbm_path,taxon,'fc'))

    create_fbm(data.holder[[taxon]]$data$zscore,
               file.path(fbm_path,taxon,'zscore'))


    create_fbm(data.holder[[taxon]]$data$adj.pv,
               file.path(fbm_path,taxon,'adj.pv'))
}


for (taxon in names(data.holder)) {
    data.holder[[taxon]]$data$zscore <- NULL
    data.holder[[taxon]]$data$adj.pv <- NULL
    data.holder[[taxon]]$data$fc <- NULL
}
saveRDS(data.holder, paste(DATADIR, 'DATA.HOLDER.light.rds', sep='/'))

