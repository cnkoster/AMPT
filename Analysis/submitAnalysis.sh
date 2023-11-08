#!/bin/bash
outdir=/data/alice/nkoster/AMPT_out/Run2_Energy/nEvents100

if [ ! -d ${outdir} ]
then
mkdir -p ${outdir}
mkdir -p ${outdir}/logs/
fi

if [ ! -d ${outdir}/logs ]
then
mkdir -p ${outdir}/logs/
fi

rm -rf ${outdir}/runStbcJob.sh
cp runStbcJob.sh ${outdir}

cd ${outdir}
qsub -o ${outdir}/logs/logOut1 -e ${outdir}/logs/logErr1 -q generic -v OUTPUTDIR=${outdir} runStbcJob.sh
cd -





