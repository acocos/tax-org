#!/usr/bin/env bash

export GRB_LICENSE_FILE=/path/to/gurobi.lic
RELATIONFILE=../taxdata/ppdb-local/relations/lexnet-distrib.db

######################
# DMST
######################
OUTDIR=../taxdata/ppdb-local/pred_taxo/dmst
mkdir -p $OUTDIR
LOGDIR=$OUTDIR/logs
mkdir -p $LOGDIR

for termfile in ../taxdata/ppdb-local/terms_test/*2.0.terms;
do
    outfile=${termfile%terms}taxo
    outfile=$OUTDIR/${outfile#../taxdata/ppdb-local/terms_test/}
    echo $outfile;
    python dmst.py -f $termfile -r $RELATIONFILE -o $OUTDIR -L $LOGDIR
done

######################
# DMST+clus
######################
OUTDIR=../taxdata/ppdb-local/pred_taxo/dmst-clus
mkdir -p $OUTDIR
LOGDIR=$OUTDIR/logs
mkdir -p $LOGDIR

for termfile in ../taxdata/ppdb-local/terms_test/*2.0.terms;
do
    outfile=${termfile%terms}taxo
    outfile=$OUTDIR/${outfile#../taxdata/ppdb-local/terms_test/}
    clusfile=${termfile#../taxdata/ppdb-local/terms_test/}
    clusfile=../taxdata/ppdb-local/clusters/${clusfile%terms}clusters
    echo $outfile;
    python dmst.py -f $termfile -r $RELATIONFILE -o $OUTDIR -L $LOGDIR -c $clusfile
done

######################
# NoCYC
######################
LAMBDA=0.7
OUTDIR=../taxdata/ppdb-local/pred_taxo/nocyc
mkdir -p $OUTDIR
LOGDIR=$OUTDIR/logs
mkdir -p $LOGDIR

for termfile in ../taxdata/ppdb-local/terms_test/*2.0.terms;
do
    outfile=${termfile%terms}taxo
    outfile=$OUTDIR/${outfile#../taxdata/ppdb-local/terms_test/}
    echo $outfile;
    python nocyc.py -f $termfile -r $RELATIONFILE -o $OUTDIR -L $LOGDIR -l $LAMBDA
done

######################
# NoCYC+clus
######################
LAMBDA=0.7
OUTDIR=../taxdata/ppdb-local/pred_taxo/nocyc-clus
mkdir -p $OUTDIR
LOGDIR=$OUTDIR/logs
mkdir -p $LOGDIR

for termfile in ../taxdata/ppdb-local/terms_test/*2.0.terms;
do
    outfile=${termfile%terms}taxo
    outfile=$OUTDIR/${outfile#../taxdata/ppdb-local/terms_test/}
    clusfile=${termfile#../taxdata/ppdb-local/terms_test/}
    clusfile=../taxdata/ppdb-local/clusters/${clusfile%terms}clusters
    echo $outfile;
    python nocyc.py -f $termfile -r $RELATIONFILE -o $OUTDIR -L $LOGDIR -c $clusfile -l $LAMBDA
done

######################
# MaxTransGraph
######################
LAMBDA=0.7
REL=entmax
OUTDIR=../taxdata/ppdb-local/pred_taxo/mtg
mkdir -p $OUTDIR
LOGDIR=$OUTDIR/logs
mkdir -p $LOGDIR

for termfile in ../taxdata/ppdb-local/terms_test/*2.0.terms;
do
    outfile=${termfile%terms}taxo
    outfile=$OUTDIR/${outfile#../taxdata/ppdb-local/terms_test/}
    echo $outfile;
    python maxtransgraph.py -t $REL -f $termfile -r $RELATIONFILE -o $OUTDIR -L $LOGDIR -l $LAMBDA
done

######################
# MaxTransForest
######################
LAMBDA=0.7
REL=entmax
OUTDIR=../taxdata/ppdb-local/pred_taxo/mtf
mkdir -p $OUTDIR
LOGDIR=$OUTDIR/logs
mkdir -p $LOGDIR

for termfile in ../taxdata/ppdb-local/terms_test/*2.0.terms;
do
    outfile=${termfile%terms}taxo
    outfile=$OUTDIR/${outfile#../taxdata/ppdb-local/terms_test/}
    echo $outfile;
    python maxtransforest.py -t $REL -f $termfile -r $RELATIONFILE -o $OUTDIR -L $LOGDIR -l $LAMBDA
done