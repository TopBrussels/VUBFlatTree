import os


for i in os.listdir(os.getcwd()):
    if not "crab_TT" in i: continue
    #if not "trial2" in i: continue
    print "\033[94mcrab resubmit -d ./%s/\033[0m"%i
    os.system("crab resubmit -d ./%s/"%i)# --siteblacklist='T2_TR_METU'"%i)
