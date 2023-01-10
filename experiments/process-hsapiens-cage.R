## Process CAGE data in human cell lines/tissues
##
## Data from ENCODEProjectCAGE package is used. This requires an earlier
## version of CAGEr (1.20).
##
library(CAGEr, lib.loc = "/usr/local/lib/R/site-library")
library(BSgenome.Hsapiens.UCSC.hg19)

data_path <- file.path("experiments", "data", "human")
# cage_data_path <- file.path("mnt", "storage", "cage_datasets", "homo_sapiens",
#                             "hg19", "cager_objects")

## Fetch data from package ENCODEprojectCAGE
if(!requireNamespace("ENCODEprojectCAGE", quietly = TRUE)){
    install.packages(file.path(data_path, "ENCODEprojectCAGE_1.0.1.tar.gz"), type = "source", repos = NULL)
}
##
library(ENCODEprojectCAGE)
data(ENCODEhumanCellLinesSamples)
# ENCODEhumanCellLinesSamples
#
###################
# Grouping/COllating CTSSs from all human cell lines with group "cell"
###################
idx <- which(ENCODEhumanCellLinesSamples$group == "cell")
idx <- idx[1:5]
encShort <- ENCODEhumanCellLinesSamples
human_cell_cage <- CAGEr::importPublicData(source="ENCODE",
                                    dataset = encShort$dataset[idx],
                                    group = rep("cell", length(idx)),
                                    sample = encShort$sample[idx])
## Plotting correlation -- this doesn't matter since we are anyway going to
## merge them
##
# pdf(file.path(data_path, "human_cellGroup_reps_corr.pdf"), width = 11, height = 11)
# human_cell_reps_corr <- CAGEr::plotCorrelation2(human_cell_cage, samples = "all",
#                                               tagCountThreshold = 1, applyThresholdBoth = FALSE,
#                                               method = "pearson")
# dev.off()

## Individual library sizes

# firstNames <- lapply(human_cell_cage@sampleLabels, function(x){
#     strsplit(x, split = "_")[[1]][1]
# })

pdf(file.path(data_path, "human_cellGroup_cage_libSizes_only5.pdf"), width = 11, height = 10)
human_cell_cage_libsizes <- tibble::tibble(CAGEr::librarySizes(human_cell_cage)/10^6,
                                         row.names = human_cell_cage@sampleLabels)
colnames(human_cell_cage_libsizes) <- c("Library_Size_(in_millions)", "Sample")
pl <- ggpubr::ggbarplot(human_cell_cage_libsizes, y = "Library_Size_(in_millions)", x = "Sample",
                        fill = names(human_cell_cage@sampleLabels),
                        x.text.angle = 90,
                        #format.scale = TRUE, yscale = "log10",
                        orientation = "horizontal", font.xtickslab = 14, font.ytickslab = 8,
                        #fill = RColorBrewer::brewer.pal(3, "Paired")[1:2],
                        width = 0.5)
pl
dev.off()

# merge?
CAGEr::mergeSamples(human_cell_cage, mergeIndex = rep(1,length(idx)),
             mergedSampleLabels = c("human_cellGroup_merged"))


CAGEr::plotReverseCumulatives(human_cell_cage, fitInRange = c(10, 10^6),
                              onePlot = TRUE,
                              xlim = c(1,10^8), ylim = c(1,10^8))

CAGEr::normalizeTagCount(human_cell_cage, method = "powerLaw",
                         fitInRange = c(5, 10^6), alpha = 1.14, T = 10^6)

CAGEr::clusterCTSS(object = human_cell_cage,
                   threshold = 1,
                   thresholdIsTpm = TRUE,
                   nrPassThreshold = 1,
                   method = "distclu",
                   maxDist = 20,
                   removeSingletons = TRUE,
                   keepSingletonsAbove = 5)

CAGEr::cumulativeCTSSdistribution(human_cell_cage, clusters = "tagClusters")
CAGEr::quantilePositions(human_cell_cage, clusters = "tagClusters",
                         qLow = 0.1, qUp = 0.9)

sample_names <- human_cell_cage@sampleLabels
flank_size_up <- c(500)
flank_size_down <- c(500)
ok_chr_names <- paste0("chr", c(1:21, "X", "Y"))
minTPM <- 10


for (sname in sample_names){
    for (flank_size_itr in seq_along(flank_size_up)){
        message("\nSample:", sname)

        bedFilename <- file.path(data_path, paste0("samarth_hsapiens_TC_sample_", sname,
                                                   "_minTPM", minTPM, "only5.bed"))
        fastaFilename <- file.path(data_path,
                                   paste0("samarth_hsapiens_TC_sample_", sname, "_minTPM", minTPM,
                                          "_flank_up", flank_size_up[flank_size_itr],
                                          "_flank_down", flank_size_down[flank_size_itr], "only5.fa"))

        cage_tc_granges_fname <- file.path(data_path,
                                        paste0("samarth_hsapiens_TC_sample_", sname, "_minTPM", minTPM,
                                               "_flank_up", flank_size_up[flank_size_itr],
                                               "_flank_down", flank_size_down[flank_size_itr],
                                               "only5.rds"))
        myCAGEobject_thisSample <-
            CAGEr::tagClusters(human_cell_cage, samples = sname,
                               returnInterquantileWidth = TRUE, qLow = 0.1, qUp = 0.9)
        myCAGEobject_thisSample$tpm <- as.numeric(myCAGEobject_thisSample$tpm)


        ## Keeping only standard chromosomes
        message("Keeping only standard chromosomes")
        ok_idx <- which(myCAGEobject_thisSample$chr %in% ok_chr_names)
        myCAGEobject_thisSample <- myCAGEobject_thisSample[ok_idx,]

        ## Filter >= 1TPM
        message("Filtering >= ", minTPM,"TPM")
        myCAGEobject_thisSample <-
            myCAGEobject_thisSample[which(myCAGEobject_thisSample[,c("tpm.dominant_ctss")] >= minTPM),]



        message("Turning into GRanges object...")
        gr_myCAGEobject_thisSample <-
            GenomicRanges::GRanges(seqnames = myCAGEobject_thisSample[,"chr"],
                                   ranges = IRanges::IRanges(start=myCAGEobject_thisSample[,"dominant_ctss"], width=1),
                                   myCAGEobject_thisSample[,"strand"])

        # ## Keeping only standard chromosomes
        # message("Keeping only standard chromosomes")
        # gr_myCAGEobject_thisSample <- keepStandardChromosomes(gr_myCAGEobject_thisSample, species="Homo_sapiens")

        seqlengths(gr_myCAGEobject_thisSample) <-
            seqlengths(BSgenome.Hsapiens.UCSC.hg19)[names(seqlengths(gr_myCAGEobject_thisSample))]
        gr_myCAGEobject_thisSample <- trim(gr_myCAGEobject_thisSample)


        message("Making promoters...")
        prom_myCAGEobject_thisSample <-
            promoters(gr_myCAGEobject_thisSample, upstream = flank_size_up[flank_size_itr],
                      downstream = flank_size_down[flank_size_itr]+1)


        prom_myCAGEobject_thisSample <- trim(prom_myCAGEobject_thisSample, use.names = TRUE)

        ## Exclude indices where the sequences have been trimmed
        omit_ids <- which(width(prom_myCAGEobject_thisSample) != ((flank_size_up[flank_size_itr]  +
                                                                       flank_size_down[flank_size_itr])+1)
        )
        message("Omitted IDs: ", length(omit_ids))#, "- ", width(ds_myCAGEobject_thisSample))
        if(length(omit_ids) > 0){
            myCAGEobject_thisSample <- myCAGEobject_thisSample[-omit_ids, ]
            prom_myCAGEobject_thisSample <- prom_myCAGEobject_thisSample[-omit_ids]
            gr_myCAGEobject_thisSample <- gr_myCAGEobject_thisSample[-omit_ids]
        }


        message("Attempting to create a DNAStringSet...")
        message("Getting sequences")

        fasta_names <- paste0("domCTSS=", seqnames(prom_myCAGEobject_thisSample), ":",
                              start(prom_myCAGEobject_thisSample), ";strand=",
                              strand(prom_myCAGEobject_thisSample), ";",
                              "up=", flank_size_up[flank_size_itr], ";",
                              "down=", flank_size_down[flank_size_itr])

        ds_myCAGEobject_thisSample <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19,
                                                       names= seqnames(gr_myCAGEobject_thisSample),
                                                       start = start(prom_myCAGEobject_thisSample),
                                                       end = end(prom_myCAGEobject_thisSample),
                                                       strand = strand(prom_myCAGEobject_thisSample))
        names(ds_myCAGEobject_thisSample) <- fasta_names

        # message("Putting together as DNAStringSet")
        # ds_myCAGEobject_thisSample <-
        #   DNAStringSet(ds_myCAGEobject_thisSample, use.names = FALSE)



        ## Write regions to BED file
        message("Writing tag cluster coordinates in BED file at: ", bedFilename)
        write.table(
            myCAGEobject_thisSample[ , c("chr", "start", "end", "interquantile_width", "tpm.dominant_ctss", "strand")],
            file = bedFilename,
            sep = "\t",
            col.names = TRUE,
            quote = FALSE,
            row.names = FALSE)



        ##
        message("Writing FASTA file at: ", fastaFilename)
        Biostrings::writeXStringSet(ds_myCAGEobject_thisSample,
                        filepath = fastaFilename,
                        format = "FASTA",
        )


        message("Saving CAGE TC for ", sname)
        saveRDS(myCAGEobject_thisSample, file = cage_tc_granges_fname)


    }
}
