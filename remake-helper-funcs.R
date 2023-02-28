fetch_and_setup_zenodo_data <- function(){

    message("Downloading data from Zenodo. This could take some time.")
    fname <- "experiments/seqArchR-manuscript-zenodo-archive.zip"
    utils::download.file(url = "https://zenodo.org/record/5055408/files/snikumbh/archR-v0.1.8.zip?download=1",
        destfile = fname, method = "wget")

    utils::unzip(fname, exdir = "experiments")

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



