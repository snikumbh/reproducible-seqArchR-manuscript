import sys

import random
from optparse import OptionParser, OptionGroup, OptionValueError
from numpy.random import randn
from numpy import ones, concatenate, array, transpose


parser = OptionParser()

parser.add_option("-T", "--ifMotifsGivenInFile", type="str", action="store", default = "0",
                        dest="ifMotifsGivenInFile", help="Specify 1 if the motifs are given in a file, 0 if a single motif is to be planted. [default: %default]")
parser.add_option("-F", "--motifsListFilename", type="str", action="store", default = "motifsListFilename.list",
                        dest="motifsListFilename", help="Name of the file providing the list of motifs.[default: %default]. If single motif (-T 0), this a motif.")
			#All should be present in the current directory
parser.add_option("-N", "--numOfSeqs", type="int", action="store", default = 1000,
                        dest="numOfSeqs", help="Number of sequences to generate.[default: %default]")
parser.add_option("-I", "--seqLenMin", type="int", action="store", default = 1000,
                        dest="seqLenMin", help="Minimum length of a sequence.[default: %default]")
parser.add_option("-X", "--seqLenMax", type="int", action="store", default = 1000,
                        dest="seqLenMax", help="Maximum length of a sequence.[default: %default]")
parser.add_option("-S", "--posStart", type="int", action="store", default = 100,
                        dest="posStart", help="Position-range start in the sequence to plant the motifs.[default: %default]")
parser.add_option("-E", "--posEnd", type="int", action="store", default = 100,
                        dest="posEnd", help="Position-range end in the sequence to plant the motifs.[default: %default]")
parser.add_option("-M", "--mutationRate", type="float", action="store", default = 0.0,
                        dest="mutrate", help="Mutation rate for the motifs.[default: %default]")
parser.add_option("-P", "--nPositionsToMutate", type="int", action="store", default = 0,
                        dest="nPosMutate", help="Number of positions to be mutated in the motifs as per the given mutation rate.[default: %default]")
parser.add_option("-O", "--outputFastaFilename", type="str", action="store", default = "output.fasta",
                        dest="outputFastaFilename", help="Name of the file to write FASTA output.[default: %default]")
parser.add_option("-A", "--appendString", type="str", action="store", default = "positive",
 			dest="appendStr", help="Append string to FASTA Id. use to specify positive/negative.[default: %default]")

(options, args) = parser.parse_args()


def motifgen(nMotifs, motifs, numseq, seqlenmin, seqlenmax, posstart, posend, mutrate, nposmutate, dummyFlag=0):
    """Generate sequences with a particular motif at a particular location.
    Also allow a possible mutation rate of the motif.
    """
    if nMotifs == 1 and dummyFlag == 0:
        metadata = 'motifgen(%s,%d,%d,%d,%d,%d,%1.2f,%d)' % (motifs, numseq, seqlenmin, seqlenmax, posstart, posend, mutrate, nposmutate)
    else:
        metadata = 'motifgen(%s,%d,%d,%d,%1.2f,%d)' % (nMotifs, numseq, seqlenmin, seqlenmax, mutrate, nposmutate)
    acgt='acgt'
    seqlist = []
    for i in range(0,numseq):
        str=[] ;
        seqlen=random.randint(seqlenmin,seqlenmax);
        for l in range(0,seqlen):
            str.append(acgt[random.randint(0,3)])

        if nMotifs > 1 or dummyFlag == 1:
            for n in range(0,nMotifs):
                motif = motifs[n]
                if posend[n] == 0:
                    #place the motif throughout the sequence, separation is given by posstart[n] value
                    pos = posstart[n]
                    while pos < seqlen: 
                        for l in range(0,len(motif)):
                            if (pos+l<seqlen) and (pos+l>=0):
                                str[pos+l-1]=motif[l].upper()
                        pos = pos + posstart[n]
                else:
                    pos=random.randint(posstart[n],posend[n]);
                    for l in range(0,len(motif)):
                        if (random.random()>=mutrate) and (pos+l<seqlen) and (pos+l>=0):
                            str[pos+l-1]=motif[l].upper()
            seqlist.append(''.join(str))
        else:
            motif = motifs
            pos=random.randint(posstart,posend);
            # Select positions to mutate
            items = range(0,len(motif)-1)
            random.shuffle(items)
            mutate_this_pos = items[0:(nposmutate)]
            print(mutate_this_pos)
            for l in range(0,len(motif)):
                if (l in mutate_this_pos and random.random()<=mutrate):
                    print("mutate_samarth")
                else:
                    if (pos+l<seqlen and pos+l>=0):
                        str[pos+l-1]=motif[l].upper()
            seqlist.append(''.join(str))
    return metadata, seqlist


if options.ifMotifsGivenInFile == "1":
    with open(options.motifsListFilename, 'r') as f:
	    motifsLines = [l.strip() for l in f.readlines() if l.find('#') == -1]
    motifsList = []
    posStart = []
    posEnd = []
    for l in motifsLines:
        motif, start, end = l.split(' ')
        motifsList.append(motif)
        posStart.append(int(start))
        posEnd.append(int(end))
    if len(motifsList) == 1:
        a, b = motifgen(len(motifsList), motifsList, options.numOfSeqs, options.seqLenMin, options.seqLenMax, posStart, posEnd, options.mutrate, options.nPosMutate, 1)
    else:
        #print str(len(motifsList)) + ' motifs given\n'
        a, b = motifgen(len(motifsList), motifsList, options.numOfSeqs, options.seqLenMin, options.seqLenMax, posStart, posEnd, options.mutrate, options.nPosMutate)
        #print '\n'.join(b)
    with open(options.outputFastaFilename, 'w') as f:
        for l in range(0,len(b)):
            f.write('>' + a + '_sequence' + str(l+1)+ '_'+options.appendStr+ '\n')
            f.write(str(b[l]))
            f.write('\n')
else:
    #motifsListFilename now holds a motif
    a, b = motifgen(1, options.motifsListFilename, options.numOfSeqs, options.seqLenMin, options.seqLenMax, options.posStart, options.posEnd, options.mutrate, options.nPosMutate)
    #print a, b
    with open(options.outputFastaFilename, 'w') as f:
        for l in range(0,len(b)):
            f.write('>' + a + '_sequence' + str(l+1)+ '_'+options.appendStr + '\n')
            f.write(str(b[l]))
            f.write('\n')

