
load_dif_exp = function(){
    dh_path = file.path(DATADIR, 'DATA.HOLDER.light.rds')
    fbm_path = file.path(DATADIR,'data_fbm')


    dif_exp_data =  readRDS(dh_path)

    for (t in names(dif_exp_data)) {
        dif_exp_data[[t]]$data <- load_fbms(file.path(fbm_path,t),suffix = '_contrast')
    }


    gc()

    return(dif_exp_data)
}

