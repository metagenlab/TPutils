
"""
job submission on vital-it using python
- creation of temporary submission files: class BSUB_Script
- possibility to execute job arrays
- possibility to wait for the completion of one or several jobs: functions waitForJobCompletion, waitForMultipleJobCompletion

Example of submission script:

#!/bin/bash
#BSUB -q normal
#BSUB -J velvet
#BSUB -R rusage[mem=4000]
#BSUB -M 4000000
#BSUB -J array[41-96:4]
#BSUB -o %I_output.txt
#BSUB -e %I_error.txt
#BSUB -n 4 -R span[ptile=4]
./my_script.py

Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
Date: 2014
# ---------------------------------------------------------------------------
"""

from shell_command import shell_command
import sys
from time import sleep
import time
from tempfile import NamedTemporaryFile
import re
import string
import random


#def is_job_completed(job_id):
#    job_info=shell_command("bjobs %d" % job_id )[0].split('\n')[1]
#    job_status=re.findall("EXIT|DONE|PEND|RUN",job_info)[0]
#    if job_status == "RUN" or job_status== "PEND":
#        return False
#    elif job_status == "DONE":
#        return True
#    else:

#        raise(Exception('bsub command failed with status: ' + str(job_status)))

def id_generator(size=6, chars=string.ascii_uppercase + string.ascii_lowercase + string.digits):
   return ''.join(random.choice(chars) for _ in range(size))


def unique(seq):
   # not order preserving                                                                                  
   set = {}                                                                                                
   map(set.__setitem__, seq, [])                                                                           
   return set.keys()                                                                                       


def get_job_status(job_id):
   try:
      merged_text= " ".join(shell_command("bjobs %d" % job_id )[0].split('\n'))
      stat= re.findall("EXIT|DONE|PEND|RUN", merged_text)
      all_status = list(set(stat))
   except:
      raise(Exception("unknown job ID"))
   # return exit only if all jobs (in case of job array) were exited
   if len(all_status) == 1 and all_status[0] == "EXIT":
      return "EXIT"
   elif "RUN" in all_status and "EXIT" in all_status:
      return "partial EXIT"
   elif "RUN" in all_status or "PEND" in all_status:
      return "RUN"
   elif "EXIT" in all_status:
      return "partial DONE"
   else:
      return "DONE"


def detailed_job_status(job_id):
   detailed_job_status = shell_command("bjobs %d" % job_id )[0].split('\n')
   detailed_exit_job_status = []
   for i in detailed_job_status:
      if "EXIT" in i:
         detailed_exit_job_status.append(i)
   return detailed_exit_job_status
   
                                                                                                           
def is_job_completed(job_id):
    try:
       job_status = get_job_status(job_id)
    except:
       print ("Error with jobID", job_id)
    print("cheking job status:", job_status)
    if job_status == "EXIT":
       #raise(Exception('all jobs with id %s failed with status: EXIT' % job_id))
       return detailed_job_status(job_id)
    elif job_status == "RUN" or job_status == "partial EXIT":       
       return False
    elif job_status == "DONE":
       return "DONE"
    elif job_status == "partial DONE":
       print ("Partial DONE")
       return detailed_job_status(job_id)
    else:
       raise(Exception("Problem with job completion for jobID: %s" % job_id))


def is_job_running(job_id):
    try:
       merged_text= " ".join(shell_command("bjobs %d" % job_id )[0].split('\n'))
       stat= re.findall("EXIT|DONE|PEND|RUN",merged_text)
       print (stat)
       all_status = unique(stat)
    except:
       print ("unknown job ID")
       return True
    if "EXIT" in all_status:
       raise(Exception('bsub command failed with status: EXIT'))
    if "PEND" in all_status:
       return False
    else:
       return True

def check_pending_jobs(job_id_list):
    job_list=list(job_id_list)
    """ check if all listed jobs are running
    if jobs are arrays, check if all jobs of the array are running"""
    while len(job_list) != 0:
        for job_id in job_list:
            if is_job_running(job_id):
                job_list.remove(job_id)
                print (job_id, "running!")
        sleep(10)
        print ("wait!")

def wait_for_job_completion(job_id):
    while(not is_job_completed(job_id)):
        print ("wait!")
        sleep(10)


def run_job_and_wait(script):
    #print jobScript
    job_id = script.launch()
    #print jobID
    wait_for_job_completion(job_id)


def run_job(script):
    job_id = script.launch()
    return(job_id)


def wait_multi_jobs(job_id_list):
    print ('job list:', job_id_list)
    job_list=list(job_id_list)
    job_exited = []
    while len(job_list) != 0:
        for job_id in job_list:
           status = is_job_completed(job_id)
           if status == "DONE":
              job_list.remove(job_id)
              print (job_id, "fertig!")
           elif status == False:
              continue
           else:
              job_exited += status
              job_list.remove(job_id)
        sleep(200)
        print ("Waiting for job completion : %s" % time.ctime())
    return job_exited


class BSUB_script(object):
    def __init__(self, command,
                 name=None,
                 mem_in_GB=2,
                 queue='normal',
                 module_list=None,
                 log_file="log.txt",
                 error_file="error.txt",
                 array_start=None,
                 array_stop=None,
                 array_step=None,
                 n_cores=None):
        self.command = command
        if queue in ['normal', 'long']:
            self.queue = queue
        else:
            self.queue = 'normal'
        self.name = name
        self.mem_in_GB = mem_in_GB
        self.module_list = module_list
        self.log_file = log_file
        self.error_file = error_file
        self.array_start = array_start
        self.array_stop = array_stop
        self.array_step = array_step
        self.n_cores = n_cores

    # lsf file for job submission 
    def __str__(self):
        script = ['#!/bin/bash']
        # XXX fixme
        # should include job id in the output name.
        # should use the proper log directory.
        script.append('#BSUB -q %s' % self.queue)
        script.append('#BSUB -q %s' % self.queue)

        if self.name:
            script.append('#BSUB -J %s' % self.name)
        if self.n_cores:
           script.append('#BSUB -n %s' % self.n_cores)
           script.append('#BSUB -R span[hosts=1]')
        if self.mem_in_GB:
            script.append('#BSUB -R rusage[mem=%d]' % (self.mem_in_GB*1000))
            script.append('#BSUB -M %d' % (self.mem_in_GB*1000000))
        if self.array_start and self.array_stop and self.array_step:
            script.append('#BSUB -J array[%d-%d:%d]' % (self.array_start,self.array_stop,self.array_step))
            script.append('#BSUB -o %I_output.txt')
            script.append('#BSUB -e %I_error.txt')
        else:
          if self.log_file:
           script.append('#BSUB -o %s' % self.log_file)
        if not self.array_start and not self.array_stop and not self.array_step and self.error_file:
           script.append('#BSUB -e %s' % self.error_file)          

        # load modules
        if type(self.module_list) == list and len(self.module_list) > 0:
            script.append('module add %s' % ' '.join(self.module_list))

        script.append(self.command)

        return '\n'.join(script) + '\n'

    def launch(self):
        #import os
        # generate temporary file
        file_name = id_generator(24) + ".sub"
        #with NamedTemporaryFile() as temp_file:
           
        # add content to temporary file
        #temp_file.write(str(self))
        with open(file_name, "w") as f:
           f.write(str(self))
        #f.close()
        # write file content to the disk
        #temp_file.flush()
        #file_name.close()
        #import time
        #time.sleep(0.5)
        # define command line (file.name contain complete path)
        command = 'bsub < ' + file_name # file_name
        # execute command line
        print('command:', command)
        (stdout, stderr, return_code) = shell_command(command)
        # close temp file
        #file.close()
        #time.sleep(15)
        print ("out:", stdout, "err:", stderr, "code:", return_code)
        if return_code == 0:
           # return job id
           print ("Job submitted:", stdout)
           job_id=re.search("\d+", stdout).group(0)
           print (job_id)
           return int(job_id)
        else:
          print (command)
          raise(Exception('bsub submission command failed with exit status: ' + str(return_code)))
