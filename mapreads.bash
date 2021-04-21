#!/usr/bin/env bash

#$1 is experiment name
#$2 is database name
#$3 is sample file
#$4 is bed feature for other sRNAs

function print_usage() {
  echo "USAGE: $0 experimentname databasename samplefile.txt otherfeatures.bed" >&2
  echo "    experimentname: Name of experiment that will be given to output files " >&2
  echo "    databasename: Name of database created by maketrnadb.bash " >&2
  echo "    samplefile.txt: tRNAscan-SE file containing tRNAs to be used " >&2
  echo "    otherfeatures.bed:  Bed file containing non-tRNA features" >&2
 
}



REALNAME=$(readlink -f $0)
SCRIPTDIR=$( cd "$( dirname "$REALNAME" )" && pwd )
mkdir -p $1


"$SCRIPTDIR/mapreads.py" --samplefile=$3 --trnafile=$2-trnatable.txt --bowtiedb=${2}-tRNAgenome --logfile=${1}/$1-mapstats.txt

#exit

"$SCRIPTDIR/countreads.py" --samplefile=$3 --ensemblgtf=${4} --maturetrnas=$2-maturetRNAs.bed --trnaloci=${2}-trnaloci.bed --removepseudo --genetypefile=${1}/${1}-genetypes.txt >${1}/$1-counts.txt
#exit /projects/lowelab/users/holmes/pythonsource/trnatest/tgirttest/hg19.gtf.gz
if [ $# -eq 5 ]
then
	Rscript "$SCRIPTDIR/analyzecounts.R" $1 ${1}/$1-counts.txt $3 $5
	Rscript "$SCRIPTDIR/makescatter.R" ${1} ${1}/$1-normalized.txt ${2}-trnatable.txt ${1}/${1}-genetypes.txt $3 $5
else
	Rscript "$SCRIPTDIR/analyzecounts.R" $1 ${1}/$1-counts.txt $3
fi
"$SCRIPTDIR/countreadtypes.py" --sizefactors=${1}/$1-SizeFactors.txt --combinereps --samplefile=$3  --maturetrnas=$2-maturetRNAs.bed --trnatable=${2}-trnatable.txt --trnaaminofile=${1}/${1}-aminocounts.txt --ensemblgtf $4 --trnaloci=${2}-trnaloci.bed   >${1}/${1}-typecounts.txt #--countfrags 
Rscript "$SCRIPTDIR/featuretypes.R" ${1}/${1}-typecounts.txt ${1}/${1}-typecounts.pdf
Rscript "$SCRIPTDIR/featuretypes.R" ${1}/${1}-aminocounts.txt ${1}/${1}-aminocounts.pdf

#exit
"$SCRIPTDIR/countfragsize.py"  --combinereps --samplefile=$3  --maturetrnas=$2-maturetRNAs.bed --countfrags --trnaloci=${2}-trnaloci.bed --bedfile $4 --trnanormfile=$1-tRNANormFactors.txt --allreadsnormfile=$1-AllreadNormFactors.txt >${1}/${1}-readlengths.txt
Rscript "$SCRIPTDIR/readlengthhistogram.R" ${1}/${1}-readlengths.txt ${1}/${1}-readlengths.pdf



#"$SCRIPTDIR/getcoverage.py" --samplefile=$3  --bedfile=$2-maturetRNAs.bed --sizefactors=${1}/$1-SizeFactors.txt --stkfile=$2-trnaalign.stk >${1}/${1}-coverage.txt
#Rscript "$SCRIPTDIR/coverageplots.R" --cov=${1}/${1}-coverage.txt --trna=${2}-trnatable.txt --samples=$3 --allcov=${1}/${1}-coverage.pdf  --modomics=${2}-modomics.txt --multicov=${1}/${1}-multipagecoverage.pdf --combinecov=${1}/${1}-combinecoverage.pdf --directory=${1}

"$SCRIPTDIR/getcoverage.py" --samplefile=$3  --bedfile=$2-maturetRNAs.bed --sizefactors=${1}/$1-SizeFactors.txt --stkfile=$2-trnaalign.stk --uniquename=${1}/${1} >${1}/${1}-coverage.txt
Rscript "$SCRIPTDIR/coverageplots.R" --cov=${1}/${1}-coverage.txt --trna=${2}-trnatable.txt --samples=$3 --allcov=${1}/${1}-coverage.pdf  --uniquename=${1}/${1} --modomics=${2}-modomics.txt --multicov=${1}/${1}-multipagecoverage.pdf --combinecov=${1}/${1}-combinecoverage.pdf --directory=${1}

#"YeastAging","featurecounts.txt","agingshort.txt"
"$SCRIPTDIR/getcoverage.py" --samplefile=$3  --bedfile=$2-trnaloci.bed --sizefactors=${1}/$1-SizeFactors.txt --stkfile=$2-trnaloci.stk --edgemargin=30 >${1}/${1}-locicoverage.txt
Rscript "$SCRIPTDIR/locuscoverage.R" --cov=${1}/${1}-locicoverage.txt --trna=${2}-trnatable.txt --samples=$3 --allcov=${1}/${1}-locicoverage.pdf --multicov=${1}/${1}-locimultipagecoverage.pdf --combinecov=${1}/${1}-locicombinecoverage.pdf 

