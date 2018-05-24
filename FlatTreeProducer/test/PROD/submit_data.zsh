#!/bin/env zsh

slist="list.txt"
pset="crabConfigTemplate.py"
ver="Test2017Analysis_Data"
prodv="/store/user/smoortga/Analysis/FlatTree/${ver}/"

rm -f crabConfig.py*

samp=()
is=1
cat ${slist} | while read line
do
  samp[${is}]=${line}
  is=$[$is+1]
done

for i in ${samp}
do
  spl=($(echo $i | tr "/" "\n"))
  pubdn=$(echo "${spl[2]}_${spl[3]}" | sed 's%-%_%g')
  nam=$(echo "${spl[1]}" | sed 's%-%_%g')
  reqn=$(echo "${nam}_${pubdn}" | sed 's%_RunIIFall17MiniAODv2.*%%g')
  cat ${pset} | sed "s%INPUTDATASET%${i}%g" \
  | sed "s%OUTLFN%${prodv}%g" \
  | sed "s%REQUESTNAME%${reqn}%g" \
  | sed "s%PUBLISHDATANAME%${pubdn}%g" \
  > crabConfig.py
  
  echo "${reqn}"
  crab submit
  
done

rm -f crabConfig.py*
