#!/bin/bash
addon=${1}
etaFlag=${2}
outdir=/data/alice/nkoster/AMPT_out/Run2_Energy/nEvent100/EPM
#Cent30_60
nGroups=1000

if [ ! -d ${outdir} ]
then
mkdir -p ${outdir}
mkdir -p ${outdir}/logs/
fi

if [ ! -d ${outdir}/logs ]
then
mkdir -p ${outdir}/logs/
fi

rm -rf ${outdir}/runStbcJob${addon}.sh
cp runStbcJob${addon}.sh ${outdir}

cd ${outdir}
qsub -o ${outdir}/logs/logOut${addon} -e ${outdir}/logs/logErr${addon} -q generic -v ETAFLAG=${etaflag},OUTPUTDIR=${outdir} runStbcJob${addon}.sh
cd -





