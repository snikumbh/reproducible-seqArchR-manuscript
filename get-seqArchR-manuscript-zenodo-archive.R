library(utils)

fname <- "experiments/seqArchR-manuscript-zenodo-archive.zip"
utils::download.file(url = "https://zenodo.org/record/5055408/files/snikumbh/archR-v0.1.8.zip?download=1",
    destfile = fname, method = "wget")

utils::unzip(fname, exdir = "experiments")
