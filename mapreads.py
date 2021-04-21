#!/usr/bin/env python

import sys
import subprocess
import argparse
from tempfile import NamedTemporaryFile
import os
import os.path
import re
from trnasequtils import *

MAXMAPS = 100


numcores = 35



def wrapbowtie2(bowtiedb, unpaired, outfile, scriptdir, trnafile, maxmaps = MAXMAPS,program = 'bowtie2', logfile = None, mapfile = None, expname = None):
    '''
    I think the quals are irrelevant, and N should be scored as only slightly better than a mismatch
    this does both
    ignoring quals forces the mismatches to have a single score, which is necessary for precise N weighting
    
     --ignore-quals --np 5
     Very sensitive is necessary due to mismatch caused by modified base misreads
    '''

    bowtiecommand = program+' -x '+bowtiedb+' -k '+str(maxmaps)+' --very-sensitive --ignore-quals --np 5 --reorder -p '+str(numcores)+' -U '+unpaired

    #print >>sys.stderr, bowtiecommand
    temploc = os.path.basename(outfile)
    #bowtiecommand = bowtiecommand + ' | '+scriptdir+'choosemappings.py '+trnafile+' | samtools sort - '+outfile
    bowtiecommand = bowtiecommand + ' | '+scriptdir+'choosemappings.py '+trnafile+' --progname='+"TRAX"+ ' --fqname=' +unpaired+' --expname='+expname + ' | samtools sort -T '+tempfile.gettempdir()+"/"+temploc+'temp - -o '+outfile+'.bam'
    #print >>sys.stderr,  bowtiecommand
    if logfile:
        print >>logfile,  bowtiecommand
    bowtierun = None
    logfile.flush()
    bowtierun = subprocess.Popen(bowtiecommand, shell = True, stderr = subprocess.PIPE)

    output = bowtierun.communicate()
    errinfo = output[1]
    if logfile is not None:
        print >>logfile, errinfo 
    logfile.flush()
    if bowtierun.returncode:
        print >>sys.stderr, "Failure to Bowtie2 map"
        print >>sys.stderr, "check logfile"
        logfile.close()
        sys.exit(1)

    rereadtotal = re.search(r'(\d+).*reads',errinfo )
    rereadunmap = re.search(r'\s*(\d+).*0 times',errinfo )
    rereadsingle = re.search(r'\s*(\d+).*exactly 1 time',errinfo )
    rereadmult = re.search(r'\s*(\d+).*>1 times',errinfo )
    if rereadtotal and rereadunmap and rereadsingle and rereadmult:
        totalreads = rereadtotal.group(1)
        unmappedreads = rereadunmap.group(1)
        singlemaps = rereadsingle.group(1)
        multmaps = rereadmult.group(1)
        return [unmappedreads,singlemaps,multmaps,totalreads]
        
    else:
        print >>sys.stderr, "Could not map "+unpaired +", check logfile"
        print >>sys.stderr, "Exiting..."
        sys.exit(1)
        return None
    

def checkheaders(bamname, fqname):
    try:
        bamfile = pysam.Samfile(bamname, "r" )
    except ValueError:
        return True
    newheader = bamfile.header
    if len(newheader["PG"]) > 1 and newheader["PG"][1]["PN"] == "TRAX":
        
        if newheader["RG"][0]["ID"] != fqname:
            return False
    return True

#print >>sys.stderr, os.path.dirname(os.path.realpath(sys.argv[0]))
#print >>sys.stderr, os.path.abspath(__file__)
def main(**argdict):
    argdict = defaultdict(lambda: None, argdict)
    scriptdir = os.path.dirname(os.path.realpath(sys.argv[0]))+"/"
    samplefilename = argdict["samplefile"]
    
    sampledata = samplefile(argdict["samplefile"])
    trnafile = argdict["trnafile"]
    logfile = argdict["logfile"]
    mapfile = argdict["mapfile"]
    bowtiedb = argdict["bowtiedb"]
    lazycreate = argdict["lazy"]
    #sys.exit()
    
    workingdir = './'
    #samplefile = open(args.samplefile)
    
    samples = sampledata.getsamples()
    
    trnafile = trnafile
        
    if logfile and lazycreate:
        logfile = open(logfile,'a')
        print >>logfile, "New mapping"
    elif logfile:
        logfile = open(logfile,'w')
    else:
        logfile = sys.stderr

    unmaps = defaultdict(int)
    singlemaps = defaultdict(int)
    multimaps = defaultdict(int)
    totalreads = defaultdict(int)
    
    if not os.path.isfile(bowtiedb+".fa"):
        print >>sys.stderr, "No bowtie2 database "+bowtiedb
        sys.exit(1)
    badsamples = list()
    for samplename in samples:
        bamfile = workingdir+samplename
        
        if lazycreate and os.path.isfile(bamfile+".bam"):   
            if not checkheaders(bamfile+".bam", sampledata.getfastq(samplename)):
                badsamples.append(bamfile+".bam")

                
            
        else:
            if os.path.isfile(bamfile+".bam"):

                if not checkheaders(bamfile+".bam", sampledata.getfastq(samplename)):
                    badsamples.append(bamfile+".bam")
    
    if len(badsamples) > 0:
        print >>sys.stderr, "Bam files "+",".join(badsamples)+" does not match fq files"
        print >>sys.stderr, "Aborting"
        sys.exit(1)               
    #'samtools sort -T '+tempfile.gettempdir()+"/"+outfile+'temp - -o '+outfile+'.bam'
    tempfilesover = list()
    missingfqfiles = list()
    for samplename in samples:
        #redundant but ensures compatibility
        bamfile = workingdir+samplename
        temploc = os.path.basename(bamfile)
        #print >>sys.stderr, "***"
        #print >>sys.stderr, samplename+'temp'
        
        for currfile in os.listdir(tempfile.gettempdir()):
            #
            if currfile.startswith(samplename+'temp'):
                tempfilesover.append(currfile)
        fqfile = sampledata.getfastq(samplename)
        if not os.path.isfile(fqfile):
            missingfqfiles.append(fqfile)
    if len(tempfilesover) > 0:
        for currfile in tempfilesover:
            print >>sys.stderr, tempfile.gettempdir() +"/"+ currfile + " temp bam files exists"
        print >>sys.stderr, "these files must be deleted to proceed"
        sys.exit(1)
    if len(missingfqfiles) > 0:
        print >>sys.stderr, ",".join(missingfqfiles) + " fastq files missing"
        sys.exit(1)
            
    for samplename in samples:
        bamfile = workingdir+samplename
        
        if lazycreate and os.path.isfile(bamfile+".bam"):
            pass
                
            print >>sys.stderr, "Skipping "+samplename
            
        else:
            if os.path.isfile(bamfile+".bam"):
                pass

            print >>logfile, "Mapping "+samplename
            print >>sys.stderr, "Mapping "+samplename
            logfile.flush()
            mapresults = wrapbowtie2(bowtiedb, sampledata.getfastq(samplename),bamfile,scriptdir, trnafile,  logfile=logfile, expname = samplefilename)

            
            if mapresults is not None:
                unmaps[samplename] = mapresults[0]
                singlemaps[samplename] = mapresults[1]
                multimaps[samplename] = mapresults[2]
                totalreads[samplename] = mapresults[3]
    
    if lazycreate:
        #here is where I might add stuff to read old files in lazy mode
        pass
    if mapfile is not None and not lazycreate:
        mapinfo = open(mapfile,'w')                
        print >>mapinfo, "\t".join(samples)
        print >>mapinfo, "unmap\t"+"\t".join(unmaps[currsample] for currsample in samples)
        print >>mapinfo, "single\t"+"\t".join(singlemaps[currsample] for currsample in samples)
        print >>mapinfo, "multi\t"+"\t".join(multimaps[currsample] for currsample in samples)
        #print >>mapinfo, "total\t"+"\t".join(totalreads[currsample] for currsample in samples)
        mapinfo.close()
    
        
        
        #print >>logfile, "Processing "+samplename +" mappings"
    logfile.close()
    
        
        
        #result = subprocess.call(scriptdir+'choosemappings.py '+trnafile+' <'+bamfile +' | samtools view -F 4 -b - | samtools sort - '+workingdir+samplename+'_sort', shell = True)
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Map reads with bowtie2 and process mappings')
    
    parser.add_argument('--samplefile',
                       help='Sample file in format')
    parser.add_argument('--trnafile',
                       help='tRNA file in format')
    parser.add_argument('--logfile',
                       help='log file for error messages and mapping stats')
    parser.add_argument('--mapfile',
                       help='output table with mapping stats')
    parser.add_argument('--bowtiedb',
                       help='Location of Bowtie 2 database')
    parser.add_argument('--lazy', action="store_true", default=False,
                       help='do not remap if mapping results exist')
    
    args = parser.parse_args()
    main(samplefile = args.samplefile, trnafile= args.trnafile, logfile = args.logfile, bowtiedb = args.bowtiedb, lazy = args.lazy, mapfile = args.mapfile)


