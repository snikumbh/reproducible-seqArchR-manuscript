# Generate toy (collection of) sequences

This is a very simple python script to generate toy (collection of) sequences with specific set of motifs with the desired mutation rates and starting positions of choice. Below is an example of how to do this. The output is a FASTA file (.fasta) with the generated sequences.

Note: This is a repurposed/rewritten version of datagen.py/datafuncs.py (just the relevant function) from the _[EasySVM](https://github.com/ratschlab/easysvm) toolbox based on Shogun_ for personal use. I rewrote it when I faced issues with installation of easysvm on a machine I was then working on.


## Usage
```
python generateData.py -T 1 -F <file_w_list_of_motifs> -N 500 -I 1000 -X 1500 -M 0.1 -O <output_filename>
```

The various command line arguments are explained in the help/usage information; fetch it via `python generateData.py -h`:

```
python generateData.py -h
Usage: generateData.py [options]

Options:
  -h, --help            show this help message and exit
  -T IFMOTIFSGIVENINFILE, --ifMotifsGivenInFile=IFMOTIFSGIVENINFILE
                        Specify 1 if the motifs are given in a file, 0 if a
                        single motif is to be planted. [default: 0]
  -F MOTIFSLISTFILENAME, --motifsListFilename=MOTIFSLISTFILENAME
                        Name of the file providing the list of
                        motifs.[default: motifsListFilename.list]. If single
                        motif (-T 0), this a motif.
  -N NUMOFSEQS, --numOfSeqs=NUMOFSEQS
                        Number of sequences to generate.[default: 1000]
  -I SEQLENMIN, --seqLenMin=SEQLENMIN
                        Minimum length of a sequence.[default: 1000]
  -X SEQLENMAX, --seqLenMax=SEQLENMAX
                        Maximum length of a sequence.[default: 1000]
  -S POSSTART, --posStart=POSSTART
                        Position-range start in the sequence to plant the
                        motifs.[default: 100]
  -E POSEND, --posEnd=POSEND
                        Position-range end in the sequence to plant the
                        motifs.[default: 100]
  -M MUTRATE, --mutationRate=MUTRATE
                        Mutation rate for the motifs.[default: 0.1]
  -O OUTPUTFASTAFILENAME, --outputFastaFilename=OUTPUTFASTAFILENAME
                        Name of the file to write FASTA output.[default:
                        output.fasta]
  -A APPENDSTR, --appendString=APPENDSTR
                        Append string to FASTA Id. use to specify
                        positive/negative.[default: positive]

```

Few examples of motifs given in files are provided (see `motif.list*`).


Any suggestions, comments, feature requests as either [#issues] or [#pull] are welcome!
