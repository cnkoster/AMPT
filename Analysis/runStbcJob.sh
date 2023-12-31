#!/bin/bash
source /cvmfs/alice.cern.ch/etc/login.sh
eval $(alienv printenv VO_ALICE@ROOT::v6-24-06-97)

cd ${TMPDIR}
echo ${TMPDIR}

export INPUT_FILES_DIR=/user/nkoster/PhD/AMPT/AnalyseTreeOutput

iStart=0
iEnd=200 

echo ${OUTPUTDIR}

cp ${INPUT_FILES_DIR}/CalculateFlowCME.h ${TMPDIR}
cp ${INPUT_FILES_DIR}/CalculateFlowCME.cxx ${TMPDIR}
cp ${INPUT_FILES_DIR}/run.C ${TMPDIR}
cp ${INPUT_FILES_DIR}/Event.h ${TMPDIR}
cp ${INPUT_FILES_DIR}/Particle.h ${TMPDIR}

pwd
ls

##run the analysis
echo "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
root.exe -b -q CalculateFlowCME.cxx run.C'(55.,'${iStart}','${iEnd}',1,500)';


if [ ! -f AnalysisResults*.root ]; then
        echo "Failed, trying again"
	#root.exe -b -q CalculateFlowCME.cxx run.C'(55.,${iStart},${iEnd},1,500)';
fi
if [ ! -f AnalysisResults*.root ]; then
        echo "Failed yet again. Second time's the charm"
        #root.exe -b -q CalculateFlowCME.cxx run.C'(55.,${iStart},${iEnd},1,500)';
fi
if [ ! -f AnalysisResults*.root ]; then
        echo "Failed yet again. Third time's the charm"
	#root.exe -b -q CalculateFlowCME.cxx run.C'(55.,${iStart},${iEnd},1,500)';
fi
if [ ! -f AnalysisResults*.root ]; then
        echo "Forget it. Clean tree AnalysisResults up, this is bogus, sorry"
        #rm *.root
fi
echo "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"

mv *.root ${OUTPUTDIR}

echo "Anything left?"
ls
echo "Done! Enjoy!"
