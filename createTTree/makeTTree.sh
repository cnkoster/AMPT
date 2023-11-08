#!/bin/sh

#  makeTTree.sh
#  
#
#  Created by Noor Koster on 13/09/2023.
#  

# Input amptTTree_v4.C(Int_t Njobs = 10, Int_t Nsplit = 1, Int_t nEventsperJob = 500)
directory=${1}
addon=${2}

#cd /dcache/alice/nkoster/PhD/AMPT_out/standard_input/nEvent500

#cd /data/alice/nkoster/AMPTtest/nEvent1
cd ${directory}

root -q /user/nkoster/PhD/AMPT/createTTree/amptTTree_v6${addon}.C

cd -
