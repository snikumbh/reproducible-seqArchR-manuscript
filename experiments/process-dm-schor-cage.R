

library(CAGEr)
library(BSgenome.Dmelanogaster.UCSC.dm6)

data_path <- file.path("experiments", "data", "drosophila-schor2017", "newly-processed")
cage_data_loc <- file.path("/mnt", "storage", "cage_datasets")
dm_cage <- file.path(cage_data_loc, "drosophila_melanogaster", "dm6")
raw_cage_bam_loc <- file.path(dm_cage, "BAM_bowtie")

## Lines from Leonie's code to pick the samples that we have been working with.
meta_data_file <- file.path(dm_cage, "meta_data",
                            "Metadata_dmelanogaster_Schor_sample_ctssAdded_post_qc.tsv")

info <- read.table(meta_data_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
info <- info[order(info$order_samples),]
sample.names <- c(info$uniqueLabel[1:9], "RAL908_6_to_8_sample12", "RAL28_10_to_12_sample12",info$uniqueLabel[14:17])

##
inputFiles <- file.path(raw_cage_bam_loc, basename(info$file_path_bam))




ce <- CAGEr::CAGEexp( metadata = list(genomeName = "BSgenome.Dmelanogaster.UCSC.dm6")
                    , colData  = DataFrame( inputFiles     = inputFiles
                                     , sampleLabels   = info$uniqueLabel
                                     , inputFilesType = "bam"
                                     , row.names      = info$uniqueLabel))


ce <- CAGEr::getCTSS(ce, removeFirstG = FALSE, useMulticore = TRUE, nrCores = 16)

# Merge samples
uniq.names <- unique(info$mergedSampleID)
sample.names <- info$mergedSampleID
mIndex <- unlist(lapply(sample.names, function(x){which(x == uniq.names)}))
ce <- CAGEr::mergeSamples(ce, mergeIndex = mIndex, mergedSampleLabels = uniq.names)


sam_pl <- ggpubr::ggbarplot(data.frame("samples" = uniq.names,
                    "LibrarySizes" = ce$librarySizes/10^6),
                    x = "samples", y = "LibrarySizes")
sam_pl <- ggpubr::ggpar(sam_pl, ylab = "Library sizes (in millions)",
                        rotate = TRUE)

pdf(file = file.path(data_path, "dm6_schor_reverse_cumulatives_5_10000.pdf"))
CAGEr::plotReverseCumulatives(ce, onePlot = TRUE, fitInRange = c(5, 10000))
dev.off()

ce <- CAGEr::normalizeTagCount(ce, method = "powerLaw",
                fitInRange = c(5, 10000), alpha = 1.08, T = 1*10^6)


# CLUSTER CTSS
ce <- CAGEr::clusterCTSS(object = ce, threshold = 0.5, thresholdIsTpm = TRUE,
                    nrPassThreshold = 1, method = "distclu", maxDist = 20,
                    removeSingletons = TRUE, keepSingletonsAbove = 5,
                    useMulticore = TRUE, nrCores = 16)
# TC interquantile widths
ce <- CAGEr::cumulativeCTSSdistribution(ce, clusters = "tagClusters",
                    useMulticore = TRUE, nrCores = 16)
ce <- CAGEr::quantilePositions(ce, clusters = "tagClusters",
                    qLow = 0.1, qUp = 0.9, useMulticore = TRUE, nrCores = 16)
# Save object
saveRDS(ce, file.path(data_path, "cager_objects/cageexp_Normalised_TagC.rds"))
# new sample labels
sample.names <- as.character(sampleLabels(ce))

########
## consensus clusters -- trial
ce <- CAGEr::aggregateTagClusters(ce, tpmThreshold = 0.5, qLow = 0.1,
                    qUp = 0.9, useMulticore = TRUE, nrCores = 16)

ce <- CAGEr::cumulativeCTSSdistribution(ce, "consensusClusters",
                    useMulticore = TRUE, nrCores = 16)
ce <- CAGEr::quantilePositions(ce, clusters = "consensusClusters",
                    qLow = 0.1, qUp = 0.9, useMulticore = TRUE, nrCores = 16)


###### create tracks ####
# Old CAGEr commands/functions
# CAGEr::exportCTSStoBedGraph(ce, values = "normalized", format = "bedGraph", oneFile = FALSE)
# CAGEr::exportToBed(object = ce, what = "tagClusters", qLow = 0.1, qUp = 0.9, oneFile = FALSE)

# Tagclusters go in BED format
# CTSSs go in BEDGRAPH format

# with new CAGEr
bedtrk <- CAGEr::exportToTrack(ce, what = "tagClusters",
                               qLow = 0.1, qUp = 0.9, oneTrack = FALSE)
for(sn in sample.names){
    pre_fname <- file.path(data_path, "tracks", paste0(sn, ".CTSS.normalized."))
    suffix <- ".bedGraph"
    trk <- CAGEr::exportToTrack(CTSSnormalizedTpmGR(ce, samples = sn), oneTrack = FALSE)

    splTrk <- split(trk, strand(trk), drop = TRUE)
    lapply(splTrk, function(x){
        p_fname <- paste0(pre_fname, "plus", suffix)
        rtracklayer::export.bedGraph(splTrk$`+`, con = p_fname, format = "bedGraph")
        m_fname <- paste0(pre_fname, "minus", suffix)
        splTrk$`-`$score <- -1*splTrk$`-`$score
        rtracklayer::export.bedGraph(splTrk$`-`, con = m_fname, format = "bedGraph", )
    })

    fname <- file.path(data_path, "tracks", paste0(sn, ".tagClusters.qLow0.1_qUp0.9.bed"))

    rtracklayer::export.bed(bedtrk[[sn]], con = fname, format = "bed")

}


flank_size_up <- c(500)
flank_size_down <- c(500)
# ok_chr_names <- paste0("chr", c(1:21, "X", "Y"))
minTPM <- 1


for (sname in sample.names){
    for (flank_size_itr in seq_along(flank_size_up)){
        message("\nSample:", sname)

        tcGR <- CAGEr::tagClustersGR(ce, sample = sname,
                                    returnInterquantileWidth = TRUE,
                                    qLow = 0.1, qUp = 0.9)



        ## Keeping only standard chromosomes
        message("Keeping only standard chromosomes")
        tcGR <- GenomeInfoDb::keepStandardChromosomes(tcGR,
                                        species = names(genomeStyles())[5],
                                        pruning.mode = "coarse")


        ## Filter >= 1TPM
        message("Filtering >= ", minTPM,"TPM")
        tcGR <- tcGR[which(tcGR$tpm.dominant_ctss >= minTPM),]


        # #This was required with previous versions of CAGEr.
        # #With version 2.0.0 onwards, tag clusters are returned as a GRanges
        # #already with the function CAGEr::tagClustersGR
        #
        message("Turning into GRanges object around the dominant tss...")
        domTcGR <- GenomicRanges::GRanges(seqnames = seqnames(tcGR),
                                   ranges = IRanges::IRanges(start=tcGR$dominant_ctss, width=1),
                                   strand = tcGR@strand)

        seqlengths(domTcGR) <- seqlengths(BSgenome.Dmelanogaster.UCSC.dm6)[names(seqlengths(tcGR))]
        domTcGR <- trim(domTcGR)


        message("Making promoters...")
        prom_tcGR <- promoters(domTcGR, upstream = flank_size_up[flank_size_itr],
                      downstream = flank_size_down[flank_size_itr]+1)
        seqlengths(prom_tcGR) <- seqlengths(BSgenome.Dmelanogaster.UCSC.dm6)[names(seqlengths(prom_tcGR))]


        prom_tcGR <- trim(prom_tcGR, use.names = TRUE)

        ## Exclude indices where the sequences have been trimmed
        omit_ids <- which(width(prom_tcGR) != ((flank_size_up[flank_size_itr]  +
                                            flank_size_down[flank_size_itr])+1)
        )
        message("Omitted IDs: ", length(omit_ids))
        if(length(omit_ids) > 0){
            tcGR <- tcGR[-omit_ids, ]
            prom_tcGR <- prom_tcGR[-omit_ids,]
            domTcGR <- domTcGR[-omit_ids,]
        }


        message("Attempting to create a DNAStringSet...")
        message("Getting sequences")

        fasta_names <- paste0("domCTSS=", seqnames(prom_tcGR), ":",
                              start(prom_tcGR), ";strand=",
                              strand(prom_tcGR), ";",
                              "up=", flank_size_up[flank_size_itr], ";",
                              "down=", flank_size_down[flank_size_itr])

        ds_tcGR <- BSgenome::getSeq(BSgenome.Dmelanogaster.UCSC.dm6,
                                       names= seqnames(prom_tcGR),
                                       start = start(prom_tcGR),
                                       end = end(prom_tcGR),
                                       strand = strand(prom_tcGR))
        names(ds_tcGR) <- fasta_names

        # message("Putting together as DNAStringSet")
        # ds_myCAGEobject_thisSample <-
        #   DNAStringSet(ds_myCAGEobject_thisSample, use.names = FALSE)

        ## Write regions to BED file
        tcGRBed <- as.data.frame(tcGR)
        colnames(tcGRBed)[1] <- "chr"
        bedFilename <- file.path(data_path, paste0("dm6_samarth_schor_et_al_TC_sample_", sname, "_minTPM", minTPM, ".bed"))
        message("Writing tag cluster coordinates in BED file at: ", bedFilename)
        write.table(
            tcGRBed[ , c("chr", "start", "end", "interquantile_width", "tpm.dominant_ctss", "strand")],
            file = bedFilename,
            sep = "\t",
            col.names = TRUE,
            quote = FALSE,
            row.names = FALSE)


        fastaFilename <- file.path(data_path,
                                   paste0("dm6_samarth_schor_et_al_TC_sample_", sname, "_minTPM", minTPM,
                                          "_flank_up", flank_size_up[flank_size_itr],
                                          "_flank_down", flank_size_down[flank_size_itr], ".fa"))
        ##
        message("Writing FASTA file at: ", fastaFilename)
        writeXStringSet(ds_tcGR,
                        filepath=fastaFilename,
                        format = "FASTA",
        )

        message("Saving CAGE TC for ", sname)
        saveRDS(tcGR,
                file = file.path(data_path,
                                 paste0("dm6_samarth_schor_et_al_TC_sample_", sname, "_minTPM", minTPM,
                                        # "_flank_up", flank_size_up[flank_size_itr],
                                        # "_flank_down", flank_size_down[flank_size_itr],
                                        ".rds")
                )
        )

    }
}
