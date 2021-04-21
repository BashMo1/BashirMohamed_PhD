#!/usr/bin/env python

import pysam
import sys
import argparse
import os.path
from collections import defaultdict
from trnasequtils import *
import itertools


def getdupes(namelist):
    allset = set()
    for currname in namelist:
        if currname in allset:
            yield currname
        else:
            allset.add(currname)
        


def main(**argdict):
    argdict = defaultdict(lambda: None, argdict)
    includebase = argdict["nofrag"]
    fullpretrnasonly = argdict["onlyfullpretrnas"]
    trnatable = argdict["trnatable"]
    removepseudo = argdict["removepseudo"]
    ensemblgtf = argdict["ensemblgtf"]
    nomultimap = argdict["nomultimap"]
    maxmismatches = argdict["maxmismatches"]
    typefile = None
    sampledata = samplefile(argdict["samplefile"])
    bedfiles = list()
    if "bedfile"  in argdict:
        bedfiles = argdict["bedfile"]
    trnalocifiles = list()
    if "trnaloci"  in argdict:
        trnalocifiles = argdict["trnaloci"]
    maturetrnas = list()
    if "maturetrnas" in argdict:
        maturetrnas = argdict["maturetrnas"]
        
    #trnalocifiles = argdict["trnaloci"]
    #maturetrnas=argdict["maturetrnas"]
    genetypefile = argdict["genetypefile"]
    trnacountfilename = argdict["trnacounts"]
    if "countfile" not in argdict or argdict["countfile"] == "stdout":
        countfile = sys.stdout
    else:
        countfile = open(argdict["countfile"], "w")
    
    trnacountfilename = argdict["trnacounts"]
    trnainfo = transcriptfile(trnatable)
    
    wholetrnas = dict()
    fivefrags = dict()
    threefrags = dict()
    trailerfrags = dict()
    otherfrags = dict()
    allfrags = dict()
    allowindels = False
    
    
    alltrnas = list()
    
    
    samplefiles = dict()
    
    
    samples = sampledata.getsamples()
    genetypes = dict()
    fullpretrnathreshold = 2
    #Grabbing all the features to count
    try:
        featurelist = list()
        trnaloci = list()
        for currfile in bedfiles:
            
            bedfeatures = list(readfeatures(currfile, removepseudo = removepseudo))
            for curr in bedfeatures:
                genetypes[curr.name] = os.path.basename(currfile)
                
            featurelist.extend(bedfeatures)
        trnalist = list()
        for currfile in trnalocifiles:
            trnaloci.extend(list(readbed(currfile)))
        for currfile in maturetrnas:
            trnalist.extend(list(readbed(currfile)))
        if ensemblgtf is not None:    
            embllist = list(readgtf(ensemblgtf, filterpsuedo = removepseudo))
        else:
            embllist = list()
    except IOError as e:
        print >>sys.stderr, e
        sys.exit()
    featcount = defaultdict(int)
    allfeats = featurelist+trnaloci+trnalist
    if len(set(curr.name for curr in allfeats)) < len(list(curr.name for curr in allfeats )):
        print >>sys.stderr, "Duplicate names in feature list"
    
    
    #featurelist = list(curr for curr in featurelist if curr.name == 'unknown20')
    alltrnas = list(curr.name for curr in featurelist)
    #print >>sys.stderr, "***"
    #setting up all the feature count dictionaries
    counts = defaultdict(lambda: defaultdict(int))
    trnacounts = defaultdict(lambda: defaultdict(int))
    trnawholecounts = defaultdict(lambda: defaultdict(int))
    trnafivecounts = defaultdict(lambda: defaultdict(int))
    trnathreecounts = defaultdict(lambda: defaultdict(int))
    trnalocuscounts = defaultdict(lambda: defaultdict(int))
    trnalocustrailercounts = defaultdict(lambda: defaultdict(int))
    partialtrnalocuscounts = defaultdict(lambda: defaultdict(int))
    fulltrnalocuscounts  = defaultdict(lambda: defaultdict(int))
    
    
    aminocounts  = defaultdict(lambda: defaultdict(int))
    anticodoncounts =  defaultdict(lambda: defaultdict(int))                                         
    
    #how much a pre-tRNA must extend off the end
    minpretrnaextend = 5
    #minimum mapq
    minmapq = 0
    if nomultimap:
        minmapq = 2
    #minimum number of reads for a feature to be reported
    minreads = 5
    for currsample in samples:
        
        currbam = sampledata.getbam(currsample)
        #print >>sys.stderr, currsample
        #doing this thing here why I only index the bamfile if the if the index file isn't there or is older than the map file
        try:
            if not os.path.isfile(currbam+".bai") or os.path.getmtime(currbam+".bai") < os.path.getmtime(currbam):
                pysam.index(""+currbam)
            bamfile = pysam.Samfile(""+currbam, "rb" )  
        except IOError as ( strerror):
            print >>sys.stderr, strerror
            sys.exit(1)
        except pysam.utils.SamtoolsError:
            print >>sys.stderr, "Can not index "+currbam
            print >>sys.stderr, "Exiting..."
            sys.exit(1)
            
        
        for currfeat in featurelist:
            #try catch is to account for weird chromosomes and the like that aren't in the genome
            #means that if I can't find a feature, I record no counts for it rather than bailing
            try:
                for currread in getbamrange(bamfile, currfeat, singleonly = nomultimap, maxmismatches = maxmismatches,allowindels = allowindels):
                    if currfeat.coverage(currread) > 10:
                        counts[currsample][currfeat.name] += 1
            except ValueError:
                pass
         
        for genename, featset in itertools.groupby(embllist,lambda x: x.data["genename"]):
            
            #pass 
            try:
                allreads =set()
                for currfeat in list(featset):
                    #print >>sys.stderr, currfeat.name
                    for currread in getbamrange(bamfile, currfeat, singleonly = nomultimap, maxmismatches = maxmismatches,allowindels = allowindels):
                        #print >>sys.stderr, "**"+currread.name
                        if currfeat.coverage(currread) > 10:
                            pass
                            counts[currsample][genename] += 1 
                            genetypes[genename] = currfeat.data["source"]
                            #print >>sys.stderr, currfeat.bedstring()
    
            except ValueError:
                pass
        for currfeat in trnaloci:
            #print >>sys.stderr,  currfeat.bedstring()
            #print >>sys.stderr,  currfeat.getdownstream(30).bedstring()
            for currread in getbamrange(bamfile, currfeat.addmargin(30), singleonly = nomultimap, maxmismatches = maxmismatches,allowindels = allowindels):
                #gotta be more than 5 bases off one end to be a true pre-tRNA
                #might want to shove these to the real tRNA at some point, but they are for now just ignored

                if currfeat.coverage(currread) > 10 and (currread.start + minpretrnaextend <= currfeat.start or currread.end - minpretrnaextend >= currfeat.end):
                    trnalocuscounts[currsample][currfeat.name] += 1
                    if currread.start + fullpretrnathreshold <  currfeat.start and currread.end - fullpretrnathreshold + 3 >  currfeat.end:
                        fulltrnalocuscounts[currsample][currfeat.name] += 1
                    else:
                        partialtrnalocuscounts[currsample][currfeat.name] += 1
                elif currfeat.getdownstream(30).coverage(currread) > 10:  #need the elif otherwise fragments that include trailer get in there
                    trnalocuscounts[currsample][currfeat.name] += 1
                    trnalocustrailercounts[currsample][currfeat.name] += 1
                else:
                    #print >>sys.stderr,  currfeat.getdownstream(30).coverage(currread)
                    pass
        
        for currfeat in trnalist:
            for currread in getbamrange(bamfile, currfeat, singleonly = nomultimap, maxmismatches = maxmismatches,allowindels = allowindels):

                if not currfeat.strand == currread.strand:
                    continue
                if not currfeat.coverage(currread) > 10:
                    continue
                curramino = trnainfo.getamino(currfeat.name)
                curranticodon = trnainfo.getanticodon(currfeat.name)
                trnacounts[currsample][currfeat.name] += 1
                    
                fragtype = getfragtype(currfeat, currread)

                if fragtype == "Whole":
                    trnawholecounts[currsample][currfeat.name] += 1
                elif fragtype == "Fiveprime":
                    trnafivecounts[currsample][currfeat.name] += 1
                elif fragtype == "Threeprime":
                    trnathreecounts[currsample][currfeat.name] += 1
                    
                if not isuniqueaminomapping(currread):
                    pass
                elif not isuniqueacmapping(currread):
                    aminocounts[currsample][curramino] += 1
                else:
                    aminocounts[currsample][curramino] += 1
                    anticodoncounts[currsample][curranticodon] += 1
                  
                            
    
    print >>countfile, "\t".join(samples)
    
    trnanames = set()
    for currfeat in trnalist:
        if max(itertools.chain((trnacounts[currsample][currfeat.name] for currsample in samples), [0])) < minreads:
            continue
        if includebase:
            print >>countfile, currfeat.name+"\t"+"\t".join(str(trnacounts[currsample][currfeat.name]) for currsample in samples)
        else:
            print >>countfile, currfeat.name+"_wholecounts\t"+"\t".join(str(trnawholecounts[currsample][currfeat.name]) for currsample in samples)
            print >>countfile, currfeat.name+"_fiveprime\t"+"\t".join(str(trnafivecounts[currsample][currfeat.name]) for currsample in samples)
            print >>countfile, currfeat.name+"_threeprime\t"+"\t".join(str(trnathreecounts[currsample][currfeat.name]) for currsample in samples)
            print >>countfile, currfeat.name+"_other\t"+"\t".join(str(trnacounts[currsample][currfeat.name] - (trnathreecounts[currsample][currfeat.name] + trnafivecounts[currsample][currfeat.name] + trnawholecounts[currsample][currfeat.name])) for currsample in samples)
        
        
    
    for currfeat in trnaloci:
        if max(itertools.chain((trnalocuscounts[currsample][currfeat.name] for currsample in samples),[0])) < minreads:
            continue
        if includebase:
            print >>countfile, currfeat.name+"\t"+"\t".join(str(trnalocuscounts[currsample][currfeat.name]) for currsample in samples)
        else:
            print >>countfile, currfeat.name+"_wholeprecounts\t"+"\t".join(str(fulltrnalocuscounts[currsample][currfeat.name]) for currsample in samples)
            print >>countfile, currfeat.name+"_partialprecounts\t"+"\t".join(str(partialtrnalocuscounts[currsample][currfeat.name]) for currsample in samples)
            print >>countfile, currfeat.name+"_trailercounts\t"+"\t".join(str(trnalocustrailercounts[currsample][currfeat.name]) for currsample in samples)
    
    #it's currently not used, but here is where I could count by amino acid or anticodon
    if typefile:
        for curramino in trnainfo.allaminos():
                print >>typefile, "AminoTotal_"+curramino+"\t"+"\t".join(str(aminocounts[currsample][curramino]) for currsample in samples)
        for currac in trnainfo.allanticodons():
                print >>typefile, "AnticodonTotal_"+currac+"\t"+"\t".join(str(anticodoncounts[currsample][currac]) for currsample in samples)
    if genetypefile is not None:
        genetypeout = open(genetypefile, "w")
    for currfeat in featurelist :
        if currfeat.name in trnanames:
            continue
        trnanames.add(currfeat.name)
        if max(counts[currsample][currfeat.name] for currsample in samples) > minreads:
            print >>countfile, currfeat.name+"\t"+"\t".join(str(counts[currsample][currfeat.name]) for currsample in samples)
            if genetypefile is not None:
                print >>genetypeout, currfeat.name+"\t"+genetypes[currfeat.name]   
    
    if genetypefile is not None:
        
        for currfeat in trnaloci:
            print >>genetypeout, currfeat.name+"_wholeprecounts"+"\t"+"trna_wholeprecounts"
            print >>genetypeout, currfeat.name+"_partialprecounts"+"\t"+"trna_partialprecounts"
            print >>genetypeout, currfeat.name+"_trailercounts"+"\t"+"trna_trailercounts"
            print >>genetypeout, currfeat.name+""+"\t"+"tRNA_locus"
        for currfeat in trnalist:
            print >>genetypeout, currfeat.name+"_wholecounts"+"\t"+"trna_wholecounts"
            print >>genetypeout, currfeat.name+"_fiveprime"+"\t"+"trna_fiveprime"
            print >>genetypeout, currfeat.name+"_threeprime"+"\t"+"trna_threeprime"
            print >>genetypeout, currfeat.name+"_other"+"\t"+"trna_other"
            print >>genetypeout, currfeat.name+""+"\t"+"tRNA"
    
    for currfeat in embllist:
        genename = currfeat.data['genename']
        if genename in trnanames:
            continue
        trnanames.add(genename)
        if genename is None:
            print >>sys.stderr, currfeat.name
            continue
        if max(counts[currsample][genename] for currsample in samples) > minreads:
            print >>countfile, genename+"\t"+"\t".join(str(counts[currsample][genename]) for currsample in samples)
            if genetypefile is not None:
                print >>genetypeout, genename+"\t"+genetypes[genename]          
            


                    
    if genetypefile is not None:
        genetypeout.close()    
            
            
            
    if trnacountfilename is not None:
        trnacountfile = open(trnacountfilename, "w")
        
        for currfeat in trnaloci:
            if max(trnalocuscounts[currsample][currfeat.name] for currsample in samples) < minreads:
                continue
            print >>trnacountfile, currfeat.name+"\t"+"\t".join(str(trnalocuscounts[currsample][currfeat.name]) for currsample in samples)
        for currfeat in trnalist:
            if max(trnacounts[currsample][currfeat.name] for currsample in samples) < minreads:
                continue
            print  >>trnacountfile, currfeat.name+"\t"+"\t".join(str(trnacounts[currsample][currfeat.name]) for currsample in samples)
        trnacountfile.close()            
        

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Generate fasta file containing mature tRNA sequences.')
    parser.add_argument('--samplefile',
                       help='Sample file in format')
    parser.add_argument('--bedfile',  nargs='+', default=list(),
                       help='bed file with non-tRNA features')
    parser.add_argument('--gtffile',  nargs='+', default=list(),
                       help='gtf file with non-tRNA features')
    parser.add_argument('--ensemblgtf',
                       help='ensembl gtf file with tRNA features')
    parser.add_argument('--trnaloci',  nargs='+', default=list(),
                       help='bed file with tRNA features')
    parser.add_argument('--maturetrnas',  nargs='+', default=list(),
                       help='bed file with mature tRNA features')
    parser.add_argument('--onlyfullpretrnas', action="store_true", default=False,
                       help='only include full pre-trnas')
    parser.add_argument('--trnatable',
                       help='table of tRNA features')
    parser.add_argument('--removepseudo', action="store_true", default=False,
                       help='remove psuedogenes from ensembl GTFs')
    parser.add_argument('--genetypefile',
                       help='Output file with gene types')
    parser.add_argument('--trnacounts',
                       help='Output file with just trna gene counts')
    parser.add_argument('--nofrag', action="store_true", default=False,
                       help='disable fragment determination')
    parser.add_argument('--nomultimap', action="store_true", default=True,
                       help='do not count multiply mapped reads')
    parser.add_argument('--maxmismatches', default=None,
                       help='Set maximum number of allowable mismatches')
    
    
    args = parser.parse_args()
        
    #main(samplefile=args.samplefile, bedfile=args.bedfile, gtffile=args.bedfile, ensemblgtf=args.ensemblgtf, trnaloci=args.trnaloci, onlyfullpretrnas=args.onlyfullpretrnas,removepseudo=args.removepseudo,genetypefile=args.genetypefile,trnacounts=args.trnacounts,maturetrnas=args.maturetrnas,nofrag=args.nofrag,nomultimap=args.nomultimap,maxmismatches=args.maxmismatches)
    argvars = vars(args)
    #argvars["countfile"] = "stdout"
    main(**argvars)
        
