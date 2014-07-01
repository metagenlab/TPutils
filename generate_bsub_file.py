
"""
job submission on vital-it using python
- creation of temporary submission files: class BSUB_Script
- possibility to execute job arrays
- possibility to wait for the completion of one or several jobs: functions waitForJobCompletion, waitForMultipleJobCompletion

"""

from shell_command import shell_command
import sys
from time import sleep
from tempfile import NamedTemporaryFile
import re


#def is_job_completed(job_id):
#    job_info=shell_command("bjobs %d" % job_id )[0].split('\n')[1]
#    job_status=re.findall("EXIT|DONE|PEND|RUN",job_info)[0]
#    if job_status == "RUN" or job_status== "PEND":
#        return False
#    elif job_status == "DONE":
#        return True
#    else:

#        raise(Exception('bsub command failed with status: ' + str(job_status)))


def unique(seq):                                                                                           
   # not order preserving                                                                                  
   set = {}                                                                                                
   map(set.__setitem__, seq, [])                                                                           
   return set.keys()                                                                                       


                                                                                                           
def is_job_completed(job_id):
    try:
       #job_info=shell_command("bjobs %d" % job_id )[0].split('\n')[1]
       merged_text= " ".join(shell_command("bjobs %d" % job_id )[0].split('\n'))
       stat= re.findall("EXIT|DONE|PEND|RUN",merged_text)
       all_status = unique(stat)
    except:
       print "unknown job ID"
       return True   
    if "EXIT" in all_status:
       raise(Exception('bsub command failed with status: EXIT'))                                          

    elif "RUN" in all_status or "PEND" in all_status:       
       return False                                                                                       
    else:                                                                                                  
       return True                                                                                         

def is_job_running(job_id):
    try:
       merged_text= " ".join(shell_command("bjobs %d" % job_id )[0].split('\n'))
       stat= re.findall("EXIT|DONE|PEND|RUN",merged_text)
       print stat
       all_status = unique(stat)
    except:
       print "unknown job ID"
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
                print job_id, "running!"
        sleep(10)
        print "wait!"








def wait_for_job_completion(job_id):
    while(not is_job_completed(job_id)):
        print "wait!"
        sleep(10)

def run_job_and_wait(script):
    #print jobScript
    job_id = script.launch()
    #print jobID
    wait_for_job_completion(job_id)

def run_job(script):
    job_id = script.launch()
    return(job_id)

#def waitForMultipleJobCompletion(jobIDlist):
#    t=0
#    for jobID in jobIDlist:
#        print t/60, "minutes"
#        while(not isJobCompleted(jobID)):                                                          
#            print "wait!"                                                                         
#            sleep(10)
#            t=t+10


def wait_multi_jobs(job_id_list):
    job_list=list(job_id_list)
    while len(job_list) != 0:
        for job_id in job_list:
            if is_job_completed(job_id):
                job_list.remove(job_id)
                print job_id, "fertig!"
        sleep(10)
        print "wait!"


class BSUB_script(object):
    def __init__(self, command, name=None, mem_in_GB=16, queue='normal', module_list=None, log_file="log.txt",error_file="error.txt",array_start=None,array_stop=None,array_step=None):
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


    # lsf file for job submission 
    def __str__(self):
        script = ['#!/bin/bash']
        # XXX fixme
        # should include job id in the output name.
        # should use the proper log directory.
        script.append('#BSUB -q %s' % self.queue)

        if self.name:
            script.append('#BSUB -J %s' % self.name)
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
        # generate temporary file
        file = NamedTemporaryFile()
        # add content to temporary file
        file.write(str(self))
        # write file content to the disk
        file.flush()
        # define command line (file.name contain complete path)
        command = 'bsub <' + file.name
        # execute command line
        (stdout, stderr, return_code) = shell_command(command)
        # close temp file
        file.close()

        if return_code == 0:
            # return job id
            job_id=re.search("\d+", stdout).group(0)
            print job_id
            return int(job_id)
        else:
            raise(Exception('bsub command failed with exit status: ' + str(return_code)))

