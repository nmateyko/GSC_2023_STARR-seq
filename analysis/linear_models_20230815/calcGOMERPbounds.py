# Calculates TF binding score using GOMER method.
# From Omar (modified from https://github.com/de-Boer-Lab/CRM2.0)

import argparse
import multiprocessing
import MYUTILS
import PWM
import sys
import warnings
from functools import partial
from tqdm import tqdm

def get_TF_params(line):
    if line is None or line == "" or line[0] == "#":
        return None
    curTF, curPFMFP = line.rstrip().split("\t")
    curPFM = PWM.loadPWM(curPFMFP)
    curPKDM = curPFM.to_PKdM()
    minKd, bestSeq = curPKDM.getMin()
    return (minKd, curPKDM, curTF)

def get_scores(line, allTFKds, ebound, verbose, allTFMotifs):
    if line is None or line == "" or line[0] == "#":
        return None
    curSeqName, curSeq = line.rstrip().split("\t")
    if verbose > 2:
        sys.stderr.write("Scanning %s...\n"%curSeqName)
    scores = []
    for i in range(0,len(allTFKds)):
        if ebound > 0:
            curScore = allTFMotifs[i].gomerScoreEBound(curSeq, allTFKds[i])
        else:
            curScore = allTFMotifs[i].gomerScore(curSeq, allTFKds[i])
        scores.append(curScore)
    return (curSeqName, scores)

def main():
    parser = argparse.ArgumentParser(description='Calculates the gomer P(bound) for a set of motifs using the minKD as [TF].')
    parser.add_argument('-if', dest='inFPTFs', metavar='<inFileTFs>', help='Input file of TFs, with two columns: TFID\tPFM_path', required=True)
    parser.add_argument('-is', dest='inFPSeqs', metavar='<inFileSeqs>', help='Input file of Sequences, tab separated (id\tseq), one per line [default = stdin]', required=False)
    parser.add_argument('-o', dest='outFP', metavar='<outFile>', help='Where to output results [default=stdout]', required=False)
    parser.add_argument('-l',dest='logFP', metavar='<logFile>',help='Where to output errors/warnings [default=stderr]', required=False)
    parser.add_argument('-t', '--cpus', type=int, dest='num_cpus', help="Number of cpus for multiprocessing", default=1)
    parser.add_argument('-v', dest='verbose', action='count', help='Verbose output?', required=False, default=0)
    parser.add_argument('-e', dest='ebound', action='count', help='Use expected binding instead of probability of binding?', required=False, default=0)
    parser.add_argument('-nh', dest='nohead', action='count', help='Don\'t print header?', required=False, default=0)

    args = parser.parse_args()
    verbose = args.verbose

    if (args.logFP is not None):
        logFile = MYUTILS.smartGZOpen(args.logFP, 'w')
        sys.stderr=logFile
    if (args.outFP is None):
        outFile= sys.stdout
    else:
        if args.verbose>0: warnings.warn("Outputting to file " + args.outFP)
        outFile = MYUTILS.smartGZOpen(args.outFP, 'w')

    inFileTFs = MYUTILS.smartGZOpen(args.inFPTFs, 'r')
    if args.inFPSeqs is None:
        inFileSeqs = sys.stdin
    else:
        inFileSeqs = MYUTILS.smartGZOpen(args.inFPSeqs, 'r')

    allTFKds = []
    allTFMotifs = []
    allTFNames = []
    if args.nohead == 0:
        outFile.write("ID")
    if verbose > 1:
        sys.stderr.write("Inputting TFs...\n")

    for line in inFileTFs:
        TF_params = get_TF_params(line)
        if not TF_params: continue
        else:
            minKd, curPKDM, curTF = TF_params
            allTFKds.append(minKd)
            allTFMotifs.append(curPKDM)
            allTFNames.append(curTF)
            if args.nohead==0:
                outFile.write("\t%s"%(curTF))
            if verbose>3:
                sys.stderr.write("Inputting %s...\n"%curTF)
    inFileTFs.close()

    if args.nohead == 0:
        outFile.write("\n")
    if verbose > 1:
        sys.stderr.write("Done.\n")
    if verbose > 1:
        sys.stderr.write("Inputting sequences...\n")

    p = multiprocessing.Pool(args.num_cpus)
    for res in tqdm(p.imap(partial(get_scores, allTFKds=allTFKds, ebound=args.ebound,
                                   verbose=verbose, allTFMotifs=allTFMotifs), inFileSeqs, chunksize=100)):
        if not res:
            continue
        else:
            curSeqName, scores = res
            outFile.write("%s"%(curSeqName))
            for score in scores:
                outFile.write("\t%g"%(score))
            outFile.write("\n")

    inFileSeqs.close()
    outFile.close()

    if verbose > 1:
        sys.stderr.write("All done.\n")
    if (args.logFP is not None):
        logFile.close()

if __name__ == "__main__":
    main()