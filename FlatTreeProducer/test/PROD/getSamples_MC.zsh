#!/bin/env zsh

tag=RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v

das_client --query="dataset dataset=/*/${tag}*/MINIAODSIM" --format=plain --limit=0 > samples_MC.txt
