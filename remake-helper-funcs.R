fetch_and_setup_zenodo_data <- function(){

    check_folders <- file.path("experiments",
        c("data", "results", "toy-seq-generation"))

    invisible(
        lapply(check_folders, function(x){
            if(!dir.exists(x)){

                message("Downloading supporting data ", basename(x),
                    " from Zenodo. This could take some time.")

                fbasenames <- paste0(basename(x), ".tar.gz")
                url_prefix <- "https://zenodo.org/record/7692742/files/"
                url_suffix <- "?download=1"
                src_urls <- paste0(url_prefix, fbasenames, url_suffix)
                dest_fnames <- file.path("experiments", fbasenames)
                utils::download.file(url = src_urls, destfile = dest_fnames,
                    method = "wget")
                utils::untar(dest_fnames, exdir = "experiments")

            }else{
                message("Required supporting data exists: ", x)
            }
        })
    )
}

get_hs_rmd <- function(){
    list("archR_human_figures.Rmd")
}

get_hs_add_rmd <- function(){
    list("archR_human_inr_KO_figures.Rmd",
        "archR_human_inr_KO_5_5_figures.Rmd",
        "archR_human_figures_150_downstream.Rmd")
}

get_dm_rmd <- function(){
    list("archR_dm_schor_figures.Rmd")
}

get_zf_rmd <- function(){
    list("archR_zf_nepal_figures.Rmd")
}


run_org_analysis <- function(fname){
    message("Building: ", fname)

    if(!is.null(fname)){
        for(fn in fname){
            rmarkdown::render(input = file.path("experiments", fn),
                output_format = "html_document")
        }
    }else{
        message("No filename provided.")
    }

}


write_article <- function(fname){
    message("Building: ", fname)

    if(!is.null(fname)){
        for(fn in fname){
            rmarkdown::render(input = file.path(fn))
        }
    }else{
        message("No filename provided.")
    }

}



