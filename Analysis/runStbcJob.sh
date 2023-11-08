#!/bin/bash
source /cvmfs/alice.cern.ch/etc/login.sh
eval $(alienv printenv VO_ALICE@ROOT::v6-24-06-97)

outdir=/dcache/alice/nkoster/PhD/AMPT_out/Run2_Energy_PbPb/nEvents100
addon='spec'

export INPUT_FILES_DIR=/user/nkoster/PhD/AMPT/createTTree
#!/bin/bash
echo ${OUTPUTDIR}

cd outdir

pwd
ls

##run the analysis
echo "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
root.exe -b -q ${INPUT_FILES_DIR}/amptTTree_v6${addon}.C

if [ ! -f TreeOutput*.root ]; then
    echo "Failed, trying again"
	root.exe -b -q ${INPUT_FILES_DIR}/amptTTree_v6${addon}.C
fi
if [ ! -f TreeOutput*.root ]; then
    echo "Failed yet again. Second time's the charm"
    root.exe -b -q ${INPUT_FILES_DIR}/amptTTree_v6${addon}.C
fi
if [ ! -f TreeOutput*.root ]; then
    echo "Failed yet again. Third time's the charm"
	root.exe -b -q ${INPUT_FILES_DIR}/amptTTree_v6${addon}.C
fi
if [ ! -f TreeOutput*.root ]; then
    echo "Forget it. Clean tree AnalysisResults up, this is bogus, sorry"
    #rm *.root
fi
echo "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"

echo "Anything left?"
ls
echo "Done! Enjoy!"
