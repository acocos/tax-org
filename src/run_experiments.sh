#!/usr/bin/env bash

export GRB_LICENSE_FILE=/path/to/gurobi.lic
RELATIONFILE=../taxdata/ppdb-local/relations/lexnet-distrib.db

######################
# DMST
######################
echo "DMST"
OUTDIR=../taxdata/ppdb-local/pred_taxo/dmst
mkdir -p $OUTDIR
LOGDIR=$OUTDIR/logs
mkdir -p $LOGDIR
# build taxonomies
for termfile in ../taxdata/ppdb-local/terms_test/*2.0.terms;
do
    outfile=${termfile%terms}taxo
    outfile=$OUTDIR/${outfile#../taxdata/ppdb-local/terms_test/}
    echo $outfile;
    python dmst.py -f $termfile -r $RELATIONFILE -o $OUTDIR -L $LOGDIR
done
# eval
for f in $(ls $OUTDIR | grep .graph);
do
    python convert_ti_output_to_formatted.py $OUTDIR/$f dmst $OUTDIR/${f/.graph/.converted.filt.taxo} ../taxdata/ppdb-local/gs_terms/wn.${f/.graph/.terms};
done
python eval_dir_simple.py $OUTDIR ../taxdata/ppdb-local/gs_taxo

######################
# DMST+clus
######################
echo "DMST+clus"
OUTDIR=../taxdata/ppdb-local/pred_taxo/dmst-clus
mkdir -p $OUTDIR
LOGDIR=$OUTDIR/logs
mkdir -p $LOGDIR
# build taxonomies
for termfile in ../taxdata/ppdb-local/terms_test/*2.0.terms;
do
    outfile=${termfile%terms}taxo
    outfile=$OUTDIR/${outfile#../taxdata/ppdb-local/terms_test/}
    clusfile=${termfile#../taxdata/ppdb-local/terms_test/}
    clusfile=../taxdata/ppdb-local/clusters/${clusfile%terms}clusters
    echo $outfile;
    python dmst.py -f $termfile -r $RELATIONFILE -o $OUTDIR -L $LOGDIR -c $clusfile
done
# eval
for f in $(ls $OUTDIR | grep .graph);
do
    python convert_ti_output_to_formatted.py $OUTDIR/$f dmst $OUTDIR/${f/.graph/.converted.filt.taxo} ../taxdata/ppdb-local/gs_terms/wn.${f/.graph/.terms};
done
python eval_dir_simple.py $OUTDIR ../taxdata/ppdb-local/gs_taxo

######################
# NoCYC
######################
echo "NoCYC"
LAMBDA=0.7
OUTDIR=../taxdata/ppdb-local/pred_taxo/nocyc
mkdir -p $OUTDIR
LOGDIR=$OUTDIR/logs
mkdir -p $LOGDIR
# build taxonomies
for termfile in ../taxdata/ppdb-local/terms_test/*2.0.terms;
do
    outfile=${termfile%terms}taxo
    outfile=$OUTDIR/${outfile#../taxdata/ppdb-local/terms_test/}
    echo $outfile;
    python nocyc.py -f $termfile -r $RELATIONFILE -o $OUTDIR -L $LOGDIR -l $LAMBDA
done
# eval
for f in $(ls $OUTDIR | grep .graph);
do
    python convert_ti_output_to_formatted.py $OUTDIR/$f nocyc $OUTDIR/${f/.graph/.converted.filt.taxo} ../taxdata/ppdb-local/gs_terms/wn.${f/.graph/.terms};
done
python eval_dir_simple.py $OUTDIR ../taxdata/ppdb-local/gs_taxo

######################
# NoCYC+clus
######################
echo "NoCYC+clus"
LAMBDA=0.7
OUTDIR=../taxdata/ppdb-local/pred_taxo/nocyc-clus
mkdir -p $OUTDIR
LOGDIR=$OUTDIR/logs
mkdir -p $LOGDIR
# build taxonomies
for termfile in ../taxdata/ppdb-local/terms_test/*2.0.terms;
do
    outfile=${termfile%terms}taxo
    outfile=$OUTDIR/${outfile#../taxdata/ppdb-local/terms_test/}
    clusfile=${termfile#../taxdata/ppdb-local/terms_test/}
    clusfile=../taxdata/ppdb-local/clusters/${clusfile%terms}clusters
    echo $outfile;
    python nocyc.py -f $termfile -r $RELATIONFILE -o $OUTDIR -L $LOGDIR -c $clusfile -l $LAMBDA
done
# eval
for f in $(ls $OUTDIR | grep .graph);
do
    python convert_ti_output_to_formatted.py $OUTDIR/$f nocyc $OUTDIR/${f/.graph/.converted.filt.taxo} ../taxdata/ppdb-local/gs_terms/wn.${f/.graph/.terms};
done
python eval_dir_simple.py $OUTDIR ../taxdata/ppdb-local/gs_taxo

######################
# MaxTransGraph
######################
echo "MaxTransGraph"
LAMBDA=0.7
REL=entmax
OUTDIR=../taxdata/ppdb-local/pred_taxo/mtg
mkdir -p $OUTDIR
LOGDIR=$OUTDIR/logs
mkdir -p $LOGDIR
# build taxonomies
for termfile in ../taxdata/ppdb-local/terms_test/*2.0.terms;
do
    outfile=${termfile%terms}taxo
    outfile=$OUTDIR/${outfile#../taxdata/ppdb-local/terms_test/}
    echo $outfile;
    python maxtransgraph.py -t $REL -f $termfile -r $RELATIONFILE -o $OUTDIR -L $LOGDIR -l $LAMBDA
done
# eval
for f in $(ls $OUTDIR | grep .graph);
do
    python convert_ti_output_to_formatted.py $OUTDIR/$f mtg $OUTDIR/${f/.graph/.converted.filt.taxo} ../taxdata/ppdb-local/gs_terms/wn.${f/.graph/.terms};
done
python eval_dir_simple.py $OUTDIR ../taxdata/ppdb-local/gs_taxo

######################
# MaxTransForest
######################
echo "MaxTransForest"
LAMBDA=0.7
REL=entmax
OUTDIR=../taxdata/ppdb-local/pred_taxo/mtf
mkdir -p $OUTDIR
LOGDIR=$OUTDIR/logs
mkdir -p $LOGDIR
# build taxonomies
for termfile in ../taxdata/ppdb-local/terms_test/*2.0.terms;
do
    outfile=${termfile%terms}taxo
    outfile=$OUTDIR/${outfile#../taxdata/ppdb-local/terms_test/}
    echo $outfile;
    python maxtransforest.py -t $REL -f $termfile -r $RELATIONFILE -o $OUTDIR -L $LOGDIR -l $LAMBDA
done
# eval
for f in $(ls $OUTDIR | grep .graph);
do
    python convert_ti_output_to_formatted.py $OUTDIR/$f mtf $OUTDIR/${f/.graph/.converted.filt.taxo} ../taxdata/ppdb-local/gs_terms/wn.${f/.graph/.terms};
done
python eval_dir_simple.py $OUTDIR ../taxdata/ppdb-local/gs_taxo