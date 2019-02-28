#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

"""
job resubmission on vital-it using python

Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
Date: 2016
# ---------------------------------------------------------------------------
"""

import generate_bsub_file


def check_vital_it_output(output_file):

    import re

    mem_kill = re.compile("TERM_MEMLIMIT: job killed after reaching LSF memory usage limit")
    limit_kill = re.compile("TERM_RUNLIMIT: job killed after reaching LSF run time limit")
    success = re.compile("^Successfully completed")
    script_match = re.compile("^\#!/bin/bash")
    arg_match = re.compile("\#")

    with open(output_file, 'r') as f:
        script= {}
        keep=False
        for line in f:
            if keep:
                if re.match(arg_match, line):
                    data = line.rstrip().split(" ")
                    script[data[1]] = data[2]
                #print len(line.rstrip())
                else:
                    script["cmd"] = line.rstrip()
                    keep = False


            if re.match(success, line):
                return "SUCCESS", script
            if re.match(mem_kill, line):
                return "MEMKILL", script
            if re.match(limit_kill, line):
                return "LIMITKILL", script
            if re.match(script_match, line):
                keep=True

    return "UNKNOWN", script

def relaunch_vital_it_job(status, cmd_data):
    import generate_bsub_file
    import shell_command
    import sys
    if status == "MEMKILL":
        mem = int(cmd_data["-M"])/1000000
        mem+=2
        script = generate_bsub_file.BSUB_script(command=cmd_data["cmd"],
                                        mem_in_GB=mem,
                                        name=cmd_data["-J"],
                                        log_file=cmd_data["-o"],
                                        error_file=cmd_data["-e"])
        generate_bsub_file.run_job(script)
        sys.stdout.write("Job %s relaunched with increased memory limit: %s GB\n" % (cmd_data["-J"], mem))
        shell_command.shell_command("rm %s" % cmd_data["-o"])
        shell_command.shell_command("rm %s" % cmd_data["-e"])

    elif status == "LIMITKILL":
        script = generate_bsub_file.BSUB_script(command=cmd_data["cmd"],
                                        mem_in_GB=int(cmd_data["-M"])/1000000,
                                        name=cmd_data["-J"],
                                        log_file=cmd_data["-o"],
                                        error_file=cmd_data["-e"],
                                        queue="long")
        generate_bsub_file.run_job(script)
        sys.stdout.write("Job %s relaunched with queue long" % cmd_data["-J"])
        shell_command.shell_command("rm %s" % cmd_data["-o"])
        shell_command.shell_command("rm %s" % cmd_data["-e"])

    elif status == "SUCCESS":
        pass
    elif status == "UNKNOWN":
        print ('status unknown for %s' % cmd_data)
    else:
        raise IOError("Uknwon LFS error status %s\n" % status)

if __name__ == '__main__':
    import argparse
    import sys

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_files', help="output LFS files", nargs="+")

    args = parser.parse_args()

    for out in args.input_files:
        status, job_data = check_vital_it_output(out)
        relaunch_vital_it_job(status, job_data)
