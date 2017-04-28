#!/bin/bash

# Argument 1: full paht directory with jobs output

inRootFileName="root/FastRotation.root"
#inRootFileName="root/frs.root"
#inRootFileName="frs_fromAnalytical.root"
rebin="false"
rebinWidth="5"
t0="74.09" #"47.25"
tS="0" # in ns
histoName="h_frs"
maxTime="700000" # in ns
samplingPeriod="5" # in ns

./bin/FFT -f ${inRootFileName} -rebin ${rebin} -rebinWidth ${rebinWidth} -t0 ${t0} -histoName ${histoName} -maxTime ${maxTime} -samplingPeriod ${samplingPeriod} -tS ${tS}
