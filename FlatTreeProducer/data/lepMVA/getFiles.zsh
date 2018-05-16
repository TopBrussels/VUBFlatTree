#!/bin/env zsh

wpath="https://raw.githubusercontent.com/peruzzim/cmgtools-lite/tree/94X_dev_M18_EGM17/TTHAnalysis/data/leptonMVA/tth"

wget "${wpath}/el_BDTG.weights.xml"
wget "${wpath}/mu_BDTG.weights.xml"
