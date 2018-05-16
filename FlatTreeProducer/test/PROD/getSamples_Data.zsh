#!/bin/env zsh

tag=Run2017*-31Mar2018-v

das_client --query="dataset dataset=/*/${tag}*/MINIAOD" --format=plain --limit=0 > samples_Data.txt
