
# Repository of supporting information for seqArchR manuscript

**Citation information:** 
<Insert here>


**seqArchR software links**
- *Bioconductor:* http://www.bioconductor.org/packages/seqArchR
- *Github link:* https://github.com/snikumbh/seqArchR
- *Zenodo link:* https://doi.org/10.5281/zenodo.5003648

**seqArchR manuscript supporting data links**
- *Zenodo link:* https://doi.org/10.5281/zenodo.7692742



1. `git clone` this repository

2. This is a `remake`-based reproducible pipeline
(originally available [here](https://github.com/richfitz/remake); a forked 
version [here](https://github.com/snikumbh/remake)).
To run the analyses, start R and when you have remake installed, just type 
```
remake::make()
```
at the R prompt. This will (re-)generate the figures by knitting/running the 
organism-specific Rmd files.

3. All script files (R scripts, Rmarkdown files, shell scripts) are available
  in this repository

Below is the detailed explanation of the folder structure and files.




==================  TOP-LEVEL ==================
```
- experiments
- manuscript
- remake-helper-funcs.R
- remake.yml

  (Rmd files: analysis/figures-related)
- archR_dm_schor_figures.Rmd
- archR_zf_nepal_figures.Rmd
- archR_human_figures.Rmd
- archR_human_figures_150_downstream.Rmd
- archR_human_figures_inr_KO_5_5_figures.Rmd
- archR_human_figures_inr_KO_figures.Rmd

  (Md files: Corresponding output of Rmd)
- archR_dm_schor_figures.md
- archR_zf_nepal_figures.md
- archR_human_figures.md
- archR_human_figures_150_downstream.md
- archR_human_figures_inr_KO_5_5_figures.md
- archR_human_figures_inr_KO_figures.md
```

All manuscript related files are available under folder `manuscript`, and all
experiments-related files (including generated figures) are available under
`experiments`. See additional details.


==================  EXPERIMENTS ==================
```
- experiments
    - data
        - simulated-data
        - drosophila-schor2017
        - drosophila-chen2014
        - zebrafish-nepal2013
        - human
        - dual-initiation-promoters
        - tissueSpecificity-tau
    
    - results
       - simulated-data
       - comparison-approaches
       - drosophila-chen2014
       - drosophila-schor2017
       - zebrafish-nepal2013
       - human
    
    - toy-seq-generation
       - All files obtained using
         git clone of `https://github.com/snikumbh/toy-seq-generation`
         (required for generating simulated data)
         
      (R scripts)
    - archR-on-simulated-data.R
    - archR-on-simulated-data-memory-usage.R
    - archR-on-human-encode.R
    - archR-on-zebrafish-nepal2013.R
    - archR-on-drosophila-schor2017.R
    - archR-on-drosophila-chen2014.R
    - run-NPLB-script.R
    - get-simulated-fasta.R
    - process-dm-schor-cage.R (data-processing)
    - process-hsapiens-cage.R (data-processing)
    - generate-simulated-fasta.R (data-processing)
    - helper-funcs.R
    
      (Shell scripts)
    - run_archR_drosophila_chen2014.sh
    - run_archR_drosophila_schor2017.sh
    - run_archR_hsapiens.sh
    - run_archR_zebrafish_nepal2012.sh
    - run_archR_simdata.sh
    - NPLB.sh
```

All figures for the paper can be generated using the
Rmarkdown (Rmd) documents in this folder. 
These Rmd files for generating figures for the manuscript can be run as is at 
your end, but some R/shell scripts can't be run as is. 
The shell scripts are used for submitting long-running parallel seqArchR runs 
(corresponding R scripts) as jobs on Slurm at our end. 
They also refer to the relevant conda environment (required for seqArchR at 
our end) that will not work at your end.

==================  MANUSCRIPT ==================
```
- manuscript
    - seqArchR_manuscript_GenomeBiology.Rmd
    - seqArchR_supplementary.Rmd
    - seqArchR-bibliography.bib
    - mystyles.sty
```

Some figures are generated directly in the manuscript Rmd file.
Where they are generated outside the manuscript Rmd files, the
code for it can be found in the corresponding organism-specific Rmd file.


