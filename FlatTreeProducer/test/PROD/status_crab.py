import os
import subprocess
import sys
from numpy import mean


def progress(finished, running, failed, total, status=''):
    bar_len = 100
    finished_len = int(round(bar_len * finished / float(total)))
    running_len = int(round(bar_len * running / float(total)))
    failed_len = int(round(bar_len * failed / float(total)))

    percents_finished = round(100.0 * finished / float(total), 1)
    percents_running = round(100.0 * running / float(total), 1)
    percents_failed = round(100.0 * failed / float(total), 1)
    percents_other = round(100.0 * (total-finished-running-failed) / float(total), 1)
    bar = '\033[92m=\033[0m' * finished_len + '\033[94m=\033[0m' * running_len + '-' * (bar_len - finished_len - running_len - failed_len) + '\033[91m=\033[0m' * failed_len
    textbar = bool(percents_finished)*(' ' * max(0,(finished_len-len(str(percents_finished))-1)) + '\033[92m' + str(percents_finished) + '%\033[0m') + bool(percents_running)*(' ' * max(0,(running_len-len(str(percents_running))-1)) + '\033[94m' + str(percents_running) + '%\033[0m') + bool(percents_other)*(' ' * max(0,(bar_len - finished_len - running_len - failed_len - len(str(percents_finished))-1)) + str(percents_other) + '%') +  bool(percents_failed)*( ' ' * max(0,(failed_len-len(str(percents_failed))-1)) + '\033[91m' + str(percents_failed) + '%\033[0m' )
    
    sys.stdout.write('[%s] \033[92m%s%s fin\033[0m \033[94m%s%s run\033[0m \033[91m%s%s fail\033[0m %s%s o %s\n' % (bar, percents_finished, '%', percents_running, '%', percents_failed, '%', percents_other, '%', status))
    sys.stdout.write('[%s]'%(textbar))
    #sys.stdout.flush()  # As suggested by Rom Ruben (see: http://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console/27871113#comment50529068_27871113)


os.system("voms-proxy-init --voms cms:/cms/becms --valid 192:0")
os.system("cp $X509_USER_PROXY /user/$USER/")

Mean_finished_percentage_array = []
for i in os.listdir(os.getcwd()):
    if not "crab_TT" in i: continue
    #if not "trial2" in i: continue
    print "\033[94mcrab status -d ./%s/\033[0m"%i

    #output=os.popen("crab status -d ./%s/"%i).read()
    #getCrabStatus(output)
    
    
    finished = 0
    running = 0
    failed = 0
        
    proc = subprocess.Popen(['crab','status', '-d', './%s/'%i],stdout=subprocess.PIPE)
    while True:
      line = proc.stdout.readline()
      if not line:
        break
      #the real code does filtering here
      if "finished" in line and "% (" in line: finished = float(line.split("%")[-2].split(" ")[-1])
      elif "running" in line and "% (" in line: running = float(line.split("%")[-2].split(" ")[-1])
      elif "failed" in line and "% (" in line: failed = float(line.split("%")[-2].split(" ")[-1])
        
    progress(finished, running, failed, 100)
    Mean_finished_percentage_array.append(finished)
    
    #os.system("crab status -d ./%s/"%i)# --siteblacklist='T2_TR_METU'"%i)
    print ""
    print ""
    print ""


print "\033[92m --> Average percentage of finished events: ", mean(Mean_finished_percentage_array), "%\033[0m"