import os

for i in os.listdir(os.getcwd()):
    if not "crab_" in i: continue
    os.system("crab resubmit -d ./%s/"%i)
