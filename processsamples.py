#!/usr/bin/env python

import sys
import argparse
import os.path
from collections import defaultdict
from trnasequtils import *
import itertools
import subprocess
import os
import time

import mapreads
import countreads
import getcoverage
import getends
import countreadtypes
import maketrackhub
from distutils.version import LooseVersion, StrictVersion




#expname is experiment name
#dbname is database name
#samplefile is sample file
#$4 is bed feature for other sRNAs




parser = argparse.ArgumentParser(description='Generate fasta file containing mature tRNA sequences.')
parser.add_argument('--experimentname',required=True,
                   help='experiment name to be used')
parser.add_argument('--databasename',required=True,
                   help='name of the tRNA database')
parser.add_argument('--samplefile',required=True,
                   help='sample file')
parser.add_argument('--ensemblgtf',
                   help='The ensembl gene list for that species')
parser.add_argument('--exppairs',
                   help='List of sample pairs to compare')
parser.add_argument('--bedfile',  nargs='+', default=list(),
                   help='Additional bed files for feature list')
parser.add_argument('--lazyremap', action="store_true", default=False,
                   help='Skip mapping reads if bam files exit')
parser.add_argument('--nofrag', action="store_true", default=False,
                   help='Omit fragment determination (Used for TGIRT mapping)')
parser.add_argument('--olddeseq', action="store_true", default=False,
                   help='Use old DESeq1 for analysis')
parser.add_argument('--nosizefactors', action="store_true", default=False,
                   help='Don\'t use Deseq size factors in plotting')
parser.add_argument('--maxmismatch',
                   help='Maximum allowed mismatches')
parser.add_argument('--mismatch',
                   help='Generate (experimental) mismatch coverage charts')

parser.add_argument('--paironly', action="store_true", default=False,
                   help='Generate only pair files (for adding a pair file after initial processing)')
parser.add_argument('--makehub', action="store_true", default=False,
                   help='make a track hub')
parser.add_argument('--maketdr', action="store_true", default=False,
                   help='create tdrs')
parser.add_argument('--makeall', action="store_true", default=False,
                   help='make both track hub and tdrs')


rlogname = "Rlog.txt"
rlogfile = open(rlogname, "w")

def runrscript(*script):
    print >>sys.stderr, "Rscript "+" ".join(script)
    retcode = subprocess.call("Rscript "+" ".join(script), shell=True, stdout = rlogfile, stderr = subprocess.STDOUT)
    print >>rlogfile, "*******************************************************" 
    print >>rlogfile, "Rscript "+" ".join(script)
    if retcode > 0:
        print >>rlogfile, script[0]+" failed"
        print >>sys.stderr, "R script "+script[0]+" failed"
        print >>sys.stderr, "Check "+rlogname+" for details"
        
        #sys.exit()
    return retcode
    

class trnadatabase:
    def __init__(self, dbname):
        self.trnatable = dbname+"-trnatable.txt"
        self.bowtiedb = dbname+"-tRNAgenome"
        self.locifile = dbname+"-trnaloci.bed"
        self.maturetrnas=dbname+"-maturetRNAs.bed"
        self.trnaalign = dbname+"-trnaalign.stk"
        self.locialign = dbname+"-trnaloci.stk"
        self.trnafasta = dbname+"-maturetRNAs.fa"
        self.modomics = dbname+"-modomics.txt"
    def test(self):
        bowtie2job = subprocess.Popen(["bowtie2","-x",self.bowtiedb, "-U", scriptdir+"test.fq"],stdout = subprocess.PIPE,stderr = subprocess.STDOUT )
        rstatsresults = rstatsjob.communicate()[0]
        if bowtie2job.returncode  != 0:
                print >>sys.stderr, "bowtie2 failed to run"

        
class expdatabase:
    def __init__(self, expname):
        self.mapinfo = expname+"/"+expname+"-mapinfo.txt"
        self.mapplot = expname+"/"+expname+"-mapinfo.pdf"
        
        self.maplog = expname+"/"+expname+"-mapstats.txt"
        self.genetypes = expname+"/"+expname+"-genetypes.txt"
        self.genecounts = expname+"/"+expname+"-counts.txt"
        self.normalizedcounts = expname+"/"+expname+"-normalized.txt"
        self.sizefactors = expname+"/"+expname+"-SizeFactors.txt"

        self.genetypecounts=expname+"/"+expname+"-typecounts.txt"
        self.genetypeplot=expname+"/"+expname+"-typecounts.pdf"

        self.genetyperealcounts=expname+"/"+expname+"-typerealcounts.txt"
        self.genetyperealplot=expname+"/"+expname+"-typerealcounts.pdf"
        
        self.trnaaminofile=expname+"/"+expname+"-aminocounts.txt"
        self.trnaaminoplot=expname+"/"+expname+"-aminocounts.pdf"
        
        self.trnalengthfile=expname+"/"+expname+"-readlengths.txt"
        self.trnalengthplot=expname+"/"+expname+"-readlengths.pdf"
        
        self.trnacoveragefile=expname+"/"+expname+"-coverage.txt"
        self.trnacoverageplot=expname+"/"+expname+"-coverage.pdf"
        self.trnacombinecoverageplot=expname+"/"+expname+"-combinecoverage.pdf"

        self.locicoveragefile=expname+"/"+expname+"-locicoverage.txt"
        self.locicoverageplot=expname+"/"+expname+"-locicoverage.pdf"
        self.locicombinecoverageplot=expname+"/"+expname+"-locicombinecoverage.pdf"
        
        self.trnamismatchfile = expname+"/"+expname+"-mismatchcoverage.txt"
        self.trnamismatchplot = expname+"/"+expname+"-mismatchcoverage.pdf"
        
        self.trnamismatchreport = expname+"/"+expname+"-mismatchreport.txt"
        
        

def mapsamples(samplefile, trnainfo,expinfo, lazyremap):
    mapreads.main(samplefile=samplefile, trnafile=trnainfo.trnatable,bowtiedb=trnainfo.bowtiedb,logfile=expinfo.maplog,mapfile=expinfo.mapinfo, lazy=lazyremap)
def countfeatures(samplefile, trnainfo,expinfo, ensgtf, bedfiles):
    countreads.main(samplefile=samplefile,ensemblgtf=ensgtf,maturetrnas=[trnainfo.maturetrnas],trnaloci=[trnainfo.locifile],removepseudo=True,genetypefile=expinfo.genetypes,trnatable=trnainfo.trnatable,countfile=expinfo.genecounts,bedfile=bedfiles, nofrag=nofrag)
def counttypes(samplefile, trnainfo,expinfo, ensgtf, bedfiles, ignoresizefactors = False):
    if not ignoresizefactors:
        
        countreadtypes.main(sizefactors=expinfo.sizefactors,combinereps= True ,samplefile=samplefile,maturetrnas=[trnainfo.maturetrnas],trnatable=trnainfo.trnatable,trnaaminofile=expinfo.trnaaminofile,ensemblgtf=ensgtf,trnaloci=[trnainfo.locifile],countfile=expinfo.genetypecounts,realcountfile=expinfo.genetyperealcounts,bedfile= bedfiles,readlengthfile =  expinfo.trnalengthfile )
        #Plot reads by gene type and tRNAs by amino acid
        runrscript(scriptdir+"/featuretypes.R",expinfo.genetypecounts,expinfo.genetypeplot)
        runrscript(scriptdir+"/featuretypes.R",expinfo.trnaaminofile,expinfo.trnaaminoplot)
        runrscript(scriptdir+"/featuretypesreal.R",expinfo.trnaaminofile,expinfo.trnaaminoplot)

        runrscript(scriptdir+"/readlengthhistogram.R",expinfo.trnalengthfile,expinfo.trnalengthplot)
    else:
        countreadtypes.main(combinereps= True ,samplefile=samplefile,maturetrnas=[trnainfo.maturetrnas],trnatable=trnainfo.trnatable,trnaaminofile=expinfo.trnaaminofile,ensemblgtf=ensgtf,trnaloci=[trnainfo.locifile],countfile=expinfo.genetypecounts,realcountfile=expinfo.genetyperealcounts,bedfile= bedfiles,readlengthfile =  expinfo.trnalengthfile)
        #Plot reads by gene type and tRNAs by amino acid
        runrscript(scriptdir+"/featuretypes.R",expinfo.genetypecounts,expinfo.genetypeplot)
        runrscript(scriptdir+"/featuretypes.R",expinfo.trnaaminofile,expinfo.trnaaminoplot)
        runrscript(scriptdir+"/featuretypesreal.R",expinfo.trnaaminofile,expinfo.trnaaminoplot)

        runrscript(scriptdir+"/readlengthhistogram.R",expinfo.trnalengthfile,expinfo.trnalengthplot)

def gettrnacoverage(samplefile, trnainfo,expinfo, ignoresizefactors = False):
    if not ignoresizefactors:
        getcoverage.main(samplefile=samplefile,bedfile=[trnainfo.maturetrnas],sizefactors=expinfo.sizefactors,stkfile=trnainfo.trnaalign,uniquename=expname+"/"+expname, allcoverage=expinfo.trnacoveragefile)
        runrscript(scriptdir+"/coverageplots.R","--cov="+expinfo.trnacoveragefile,"--trna="+trnainfo.trnatable,"--samples="+samplefile,"--allcov="+expinfo.trnacoverageplot,"--uniquename="+expname+"/"+expname,"--modomics="+trnainfo.modomics,"--combinecov="+expinfo.trnacombinecoverageplot,"--directory="+expname)
    else:
        getcoverage.main(samplefile=samplefile,bedfile=[trnainfo.maturetrnas],stkfile=trnainfo.trnaalign,uniquename=expname+"/"+expname, allcoverage=expinfo.trnacoveragefile)
        runrscript(scriptdir+"/coverageplots.R","--cov="+expinfo.trnacoveragefile,"--trna="+trnainfo.trnatable,"--samples="+samplefile,"--allcov="+expinfo.trnacoverageplot,"--uniquename="+expname+"/"+expname,"--modomics="+trnainfo.modomics,"--combinecov="+expinfo.trnacombinecoverageplot,"--directory="+expname)


def getendscoverage(samplefile, trnainfo,expinfo, ignoresizefactors = False):
    if not ignoresizefactors:
        getends.main(samplefile=samplefile,bedfile=[trnainfo.maturetrnas],sizefactors=expinfo.sizefactors,stkfile=trnainfo.trnaalign,uniquename=expname+"/"+expname, allmismatch=expinfo.trnamismatchfile,trnafasta = trnainfo.trnafasta,mismatchfile=expinfo.trnamismatchfile,mismatchreport=expinfo.trnamismatchreport )
        runrscript(scriptdir+"/endplots.R","--cov="+expinfo.trnamismatchfile,"--mismatchcov="+expinfo.trnamismatchfile,"--trna="+trnainfo.trnatable,"--samples="+samplefile,"--allcov="+expinfo.trnamismatchplot,"--uniquename="+expname+"/"+expname,"--modomics="+trnainfo.modomics,"--directory="+expname)
        runrscript(scriptdir+"/boxplotmismatches.R","--mismatch="+expinfo.trnamismatchreport,"--trna="+trnainfo.trnatable,"--samples="+samplefile,"--directory="+expname)
    else:
        getends.main(samplefile=samplefile,bedfile=[trnainfo.maturetrnas],stkfile=trnainfo.trnaalign,uniquename=expname+"/"+expname, allmismatch=expinfo.trnamismatchfile,trnafasta = trnainfo.trnafasta,mismatchfile=expinfo.trnamismatchfile,mismatchreport=expinfo.trnamismatchreport )
        runrscript(scriptdir+"/endplots.R","--cov="+expinfo.trnamismatchfile,"--mismatchcov="+expinfo.trnamismatchfile,"--trna="+trnainfo.trnatable,"--samples="+samplefile,"--allcov="+expinfo.trnamismatchplot,"--uniquename="+expname+"/"+expname,"--modomics="+trnainfo.modomics,"--directory="+expname)
        runrscript(scriptdir+"/boxplotmismatches.R","--mismatch="+expinfo.trnamismatchreport,"--trna="+trnainfo.trnatable,"--samples="+samplefile,"--directory="+expname)

def getlocuscoverage(samplefile, trnainfo,expinfo, ignoresizefactors = False):
    if not ignoresizefactors:
        getcoverage.main(samplefile=samplefile ,bedfile=[trnainfo.locifile],sizefactors=expinfo.sizefactors,stkfile=trnainfo.locialign,edgemargin=30, uniquegenome=expname+"/"+expname+"loci",allcoverage=expinfo.locicoveragefile) #removed minextend = 5
        runrscript(scriptdir+"/locuscoverage.R", "--cov="+expinfo.locicoveragefile,"--trna="+trnainfo.trnatable,"--samples="+samplefile,"--allcov="+expinfo.locicoverageplot,"--combinecov="+expinfo.locicombinecoverageplot,"--directory="+expname)
    else:
        getcoverage.main(samplefile=samplefile ,bedfile=[trnainfo.locifile],stkfile=trnainfo.locialign,edgemargin=30, uniquegenome=expname+"/"+expname+"loci",allcoverage=expinfo.locicoveragefile)
        runrscript(scriptdir+"/locuscoverage.R", "--cov="+expinfo.locicoveragefile,"--trna="+trnainfo.trnatable,"--samples="+samplefile,"--allcov="+expinfo.locicoverageplot,"--combinecov="+expinfo.locicombinecoverageplot,"--directory="+expname)
def gettdrinfo(samplefile, dbname,expname):
    
    tdrcommand = " ".join(["bash",scriptdir+"/"+"tdrtrax.bash", samplefile, dbname,expname+"/"+expname+"tdrs"])
    print >>sys.stderr, tdrcommand
    tdrjob = subprocess.Popen(tdrcommand,stdout = subprocess.PIPE,stderr = subprocess.STDOUT, shell=True )
    print >>sys.stderr, tdrjob.communicate()[0]

def createtrackhub(samplefile, trnainfo,expinfo):
    maketrackhub.main(genomedatabase=trnainfo, samplefile=samplefile,expname=expname)

        

args = parser.parse_args()
dbname = args.databasename
expname = args.experimentname
pairfile =  args.exppairs
ensgtf = args.ensemblgtf
samplefilename = args.samplefile
lazyremap = args.lazyremap
bedfiles= args.bedfile
nofrag= args.nofrag
nosizefactors = args.nosizefactors
olddeseq = args.olddeseq
mismatch = args.mismatch 
paironly= args.paironly

makehubs = args.makehub 
maketdrs= args.maketdr

if args.makeall:
    makehubs = True
    maketdrs = True


scriptdir = os.path.dirname(os.path.realpath(sys.argv[0]))+"/"



def testsamtools(): #Version: 1.6
    samversionre = re.compile(r"Version\:\s*([\.\d]+)")
    samtoolsloc = get_location("samtools")
    if samtoolsloc is None:
            print >>sys.stderr, "Cannot find samtools in path"
            print >>sys.stderr, "Make sure samtools is installed"
    samtoolsjob = subprocess.Popen([samtoolsloc,"--help"],stdout = subprocess.PIPE,stderr = subprocess.STDOUT )
    samtoolsresults = samtoolsjob.communicate()[0]
    if samtoolsjob.returncode  != 0:
            print >>sys.stderr, "Samtools failed to run"
            print >>sys.stderr, "Make sure samtools is functioning" 
    samtoolsres = samversionre.search(samtoolsresults)
    if samtoolsres:
        if LooseVersion(samtoolsres.group(1)) < LooseVersion("1.0.0"):
            print >>sys.stderr, "Old samtools version "+samtoolsres.group(1)+" found"
            print >>sys.stderr, "Upgrade to latest version"
            sys.exit(1)
    else:
        print >>sys.stderr, "Could not find samtools version number"
        
def testrstats():
    rstatsversionre = re.compile(r"R\s+version\s+((\d+)\.(\d+)\.(\d+))")
    rstatsloc = get_location("R")
    if rstatsloc is None:
            print >>sys.stderr, "Cannot find R in path"
            print >>sys.stderr, "Make sure R is installed"
            sys.exit(1)
    rstatsjob = subprocess.Popen([rstatsloc, "--version"],stdout = subprocess.PIPE,stderr = subprocess.STDOUT )
    rstatsresults = rstatsjob.communicate()[0]
    if rstatsjob.returncode  != 0:
            print >>sys.stderr, "R failed to run"
            print >>sys.stderr, "Make sure R is functioning" 
    rstatsres = rstatsversionre.search(rstatsresults)
    if rstatsres:
        if LooseVersion(rstatsres.group(1)) < LooseVersion("3.1.2"):
            print >>sys.stderr, "Old R version "+rstatsres.group(1)+" found"
            print >>sys.stderr, "Upgrade to latest version"
            sys.exit(1)
    else:
        print >>sys.stderr, "Could not find R version number"





        
        
testrstats()
get_location("Rscript")

testsamtools()
get_location("bowtie2")

gitversion, gitversionhash = getgithash(scriptdir)

#trnainfo.test(trnainfo)


sampledata = samplefile(samplefilename)
samples = sampledata.getsamples()
for currsample in samples:
    if '-' in currsample:
        print >>sys.stderr, "Sample names containing '-' character are not allowed"
        sys.exit(1)
    if currsample[0].isdigit():
        print >>sys.stderr, "Sample names starting with digits are not allowed"
        sys.exit(1)
replicates = sampledata.allreplicates()
for currsample in replicates:
    if '-' in currsample:
        print >>sys.stderr, "Sample names containing '-' character are not allowed"
        sys.exit(1)
    if currsample[0].isdigit():
        print >>sys.stderr, "Sample names starting with digits are not allowed"
        sys.exit(1)


if runrscript(scriptdir+"checkRmodules.R") > 0:
    print >>sys.stderr, "Not all R modules needed are installed"
    print >>sys.stderr, "check README for needed R modules"
    sys.exit(1)
    





#mkdir -p expname
if not os.path.exists(expname):
    os.makedirs(expname)
if not os.path.exists(expname+"/indiv"):
    os.makedirs(expname+"/indiv")

    
    
dbname = os.path.expanduser(dbname)
if ensgtf is not None:
    ensgtf = os.path.expanduser(ensgtf)
bedfiles = list(os.path.expanduser(curr) for curr in bedfiles)


trnainfo = trnadatabase(dbname)
expinfo = expdatabase(expname)
getsamples = samplefile(samplefilename)
if len(getsamples.getsamples()) == 1:
    nosizefactors = True
    

#getendscoverage(samplefilename, trnainfo,expinfo, nosizefactors)

#if only pairfile
if pairfile and paironly:
    if olddeseq:
        deseqret = runrscript(scriptdir+"/deseq1.R",expname,expinfo.genecounts,samplefilename)
        if deseqret == 2:
            print >>sys.stderr, "Deseq analysis failed, cannot continue"
            sys.exit(1)    
    else:
        deseqret = runrscript(scriptdir+"/analyzecounts.R",expname,expinfo.genecounts,samplefilename, pairfile)
        print >>sys.stderr, scriptdir+"/analyzecounts.R",expname,expinfo.genecounts,samplefilename

        if deseqret == 2:
            print >>sys.stderr, "Deseq analysis failed, cannot continue"
            sys.exit(1)
    
	runrscript(scriptdir+"/makescatter.R",expname,expinfo.normalizedcounts,trnainfo.trnatable,expinfo.genetypes,samplefilename,pairfile)
    sys.exit(0)
elif paironly:
    print >>sys.stderr, "pair only mode used but no --pairfile used"
    sys.exit(1)





#Map the reads
runtime = time.time()
loctime = time.localtime(runtime)
print >>sys.stderr, "Mapping Reads"
#need to check here for names with dashes
mapsamples(samplefilename, trnainfo,expinfo, lazyremap)

runinfoname = expname+"/"+expname+"-runinfo.txt"
dbinfo = None
if not lazyremap:
    dbinfo = open(runinfoname,"w")
    print >>dbinfo, "Starting"

else:
    dbinfo = open(runinfoname,"a")
    print >>dbinfo, "---------------------------------------------------------"
    print >>dbinfo, "redoing"
    
print >>dbinfo, "expname\t"+expname
print >>dbinfo, "time\t"+str(runtime)+" ("+str(loctime[1])+"/"+str(loctime[2])+"/"+str(loctime[0])+")"
print >>dbinfo, "samplefile\t"+os.path.realpath(samplefilename)
print >>dbinfo, "dbname\t"+os.path.realpath(dbname)
print >>dbinfo, "git version\t"+gitversion

print >>dbinfo, "git version hash\t"+gitversionhash
dbinfo.close()

runrscript(scriptdir+"/featuretypes.R",expinfo.mapinfo,expinfo.mapplot)

#Count the reads for DEseq2 and scatter plots
print >>sys.stderr, "Counting Reads"
countfeatures(samplefilename, trnainfo,expinfo, ensgtf, bedfiles)
#Create a plot of mapped reads                                
print >>sys.stderr, "Analyzing counts"



#Analyze counts and create scatter plots if pair file is provided
if pairfile:
    if olddeseq:
        deseqret = runrscript(scriptdir+"/deseq1.R",expname,expinfo.genecounts,samplefilename)
        if deseqret == 2:
            print >>sys.stderr, "Deseq analysis failed, cannot continue"
            sys.exit(1)    
    else:
        deseqret = runrscript(scriptdir+"/analyzecounts.R",expname,expinfo.genecounts,samplefilename, pairfile)
        print >>sys.stderr, scriptdir+"/analyzecounts.R",expname,expinfo.genecounts,samplefilename

        if deseqret == 2:
            print >>sys.stderr, "Deseq analysis failed, cannot continue"
            sys.exit(1)
    
	runrscript(scriptdir+"/makescatter.R",expname,expinfo.normalizedcounts,trnainfo.trnatable,expinfo.genetypes,samplefilename,pairfile)
elif not nosizefactors:
    if olddeseq:
        deseqret = runrscript(scriptdir+"/deseq1.R",expname,expinfo.genecounts,samplefilename)
        if deseqret == 2:
            print >>sys.stderr, "Deseq analysis failed, cannot continue"
            sys.exit(1)    
    else:
        deseqret = runrscript(scriptdir+"/analyzecounts.R",expname,expinfo.genecounts,samplefilename)
        #print >>sys.stderr, scriptdir+"/analyzecounts.R",expname,expinfo.genecounts,samplefilename
        if deseqret == 2:
            print >>sys.stderr, "Deseq analysis failed, cannot continue"
            sys.exit(1)
#Count the reads by gene type
print >>sys.stderr, "Counting Read Types"
counttypes(samplefilename, trnainfo,expinfo, ensgtf, bedfiles, ignoresizefactors = nosizefactors)


#coverage plot of tRNAs
print >>sys.stderr, "Generating Read Coverage plots"      
gettrnacoverage(samplefilename, trnainfo,expinfo, ignoresizefactors = nosizefactors)

#coverage plot of pre-tRNAs
getlocuscoverage(samplefilename, trnainfo,expinfo, nosizefactors)

#coverage plot of mismatches
#getmismatchcoverage(samplefilename, trnainfo,expinfo, nosizefactors)
print >>sys.stderr, "Counting mismatches"      

getendscoverage(samplefilename, trnainfo,expinfo, nosizefactors)

if makehubs:
    print >>sys.stderr, "Creating trackhub"      

    createtrackhub(samplefilename, dbname,expname)


if (os.path.isfile(scriptdir+"/"+"tdrtrax.bash") and maketdrs):
    print >>sys.stderr, "Creating tdrs"      

    gettdrinfo(samplefilename, dbname,expname)
        #tdrtrax.bash samplefile.txt traxdb outputname


