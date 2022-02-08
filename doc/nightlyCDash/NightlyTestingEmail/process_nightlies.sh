#!/bin/bash

rm -rf TriBITS
rm -rf test_history
rm -rf *html
rm -rf *out
rm -rf *json
rm -rf  clderaNightlyBuildsTwoif.csv 

git clone git@github.com:TriBITSPub/TriBITS.git


now=$(date +"%Y-%m-%d")

./cldera_cdash_status.sh --date=$now --email-from-address=ikalash@solo-login1.sandia.gov --send-email-to=ikalash@sandia.gov,bhillma@sandia.gov,lbertag@sandia.gov


