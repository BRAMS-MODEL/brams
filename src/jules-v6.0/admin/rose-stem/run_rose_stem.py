#!/usr/bin/env python
# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************

# Owner: JULES System Manager
# Purpose: Checkout the trunk and run the rose-stem suite for nightly JULES
# testing


from optparse import OptionParser
import os
import re
import subprocess
import sys

HOME = os.environ['HOME']
USER = os.environ['USER']
LOCALDATA = "/data/local/%s" % (USER)
LOGDIR = "%s/automated_testing_jules/rose" % (HOME)
FCM = "/opt/ukmo/utils/bin/fcm"
ROSE = "/opt/ukmo/utils/bin/rose"


def parse_options():
    parser = OptionParser()
    parser.add_option("-g", "--group", dest="group", action="append",
                      help="Task group to run")
    parser.add_option("-k", "--keep", dest="keep", action="store_true",
                      help="Run with housekeeping OFF")
    parser.add_option("-n", "--name", dest="name", action="store",
                      help="Name of suite")
    parser.add_option("-s", "--source", dest="source", action="append",
                      help="Source tree to run")
    parser.add_option("-t", "--task", dest="group", action="append",
                      help="Task group to run")
    parser.add_option("-L", "--linuxmachine", dest="linuxmachine",
                      action="store",
                      help="Machine or group of linux machines to run on")
    parser.add_option("-N", "--next", dest="version", action="store_true",
                      help="Use Rose/Cylc/FCM next")
    parser.add_option("-X", "--setxcs", dest="xcs", action="store_true",
                      help="Use xcs via METO_HPC_GROUP=\'xcs\'")

    (opts, args) = parser.parse_args()
    return opts, args


def fcm_version_next(host, fname, version):
    '''Use site versions unless opt provided, change in rose-suite.conf.'''
    lines = read_remote_file(host, fname)
    newlines = []
    for line in lines:
        if re.search(r'FCM_VERSION', line):
            line = "FCM_VERSION='%s'" % (version)
        newlines.append(line)
    write_remote_file(host, fname, newlines, newline=True)


def meto_hpc(host, fname, hpc):
    '''Use site versions unless opt provided, change in rose-suite.conf.'''
    lines = read_remote_file(host, fname)
    newlines = []
    for line in lines:
        if re.search(r'METO_HPC_GROUP', line):
            line = "METO_HPC_GROUP='%s'" % (hpc)
        newlines.append(line)
    write_remote_file(host, fname, newlines, newline=True)

def read_remote_file(host, fname):
    '''Return contents of a file given filename on a host.'''
    rc, lines, stderr = run_remote_command(host, "cat %s"%(fname))

    # Reformat stdout
    lines = ''.join(lines)
    lines = lines.split("\n")
    return lines


def remote_host(opts):
    '''Define the remote host from the opt if provided.'''
    host = 'linux'
    if opts.linuxmachine:
        host = opts.linuxmachine
    rc, stdout, stderr = run_command(". /etc/profile; %s host-select %s"
                                     % (ROSE, host))
    return stdout


def run_command(command):
    '''Create the Popen object and connect out and err to pipes.'''
    p = subprocess.Popen(command, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         shell=True)

    # Do the communiate and wait to get the results of the command
    stdout, stderr = p.communicate()
    rc = p.wait()
    return rc, stdout, stderr


def run_remote_command(host, command):
    command = "ssh -Y %s '. /etc/profile; %s'" % (host, command)
    return run_command(command)


def write_file(fname, lines, newline=False):
    '''Write a file given names and contents. The optional newline argument
    adds a newline at the end of each element of the list.'''
    with open(fname, 'w') as fh:
        for line in lines:
            if newline:
                fh.write("%s\n" % (line))
            else:
                fh.write(line)


def write_remote_file(host, fname, lines, newline=False):
    '''Write a file given names and contents. The optional newline argument
    adds a newline at the end of each element of the list.'''
    tmpfile = 'run_rose_stem.tmp'
    write_file(tmpfile, lines, newline)
    cmd = "scp %s %s:%s"%(tmpfile, host, fname)
    rc, stdout, stderr, = run_command(cmd)
    if rc:
        print "Failed to run scp:", cmd
        print stderr
        sys.exit(1)
    else:
        os.unlink(tmpfile)


if __name__ == '__main__':
    opts, args = parse_options()

    host = remote_host(opts)
    host = re.sub(r'\n', r'', host)

    if not opts.group:
        print "ERROR: No task group specified."
        sys.exit(1)

    if not opts.source:
        print "ERROR: No source specified."
        sys.exit(1)

    if not opts.name:
        opts.name = opts.group[0]

    url = opts.source[0]

    workdir = "%s/rose_%s" % (LOCALDATA, opts.name)
    logfilename = "%s/rose_%s.log" % (LOGDIR, opts.name)

    # Create the dir for the log to be saved to.
    rc, stdout, stderr = run_remote_command(host, "mkdir -p %s" % (LOGDIR))
    with open(logfilename, 'w') as logfile:
        logfile.write("Running on host %s\n" % (host))
        vertxt = ''
        if opts.version:
            vertxt = 'export ROSE_VERSION=next; export CYLC_VERSION=next; '
            logfile.write("Using rose-next and cylc-next\n")

        logfile.write("Removing any existing directory %s\n" % (workdir))
        run_remote_command(host, "rm -rf %s" % (workdir))

        logfile.write("Creating directory %s\n" % (workdir))
        rc, stdout, stderr = run_remote_command(host, "mkdir -p %s"
                                                % (workdir))
        if rc is 0:
            logfile.write(stdout)
        else:
            logfile.write("mkdir failed: RC=%s\n" % (rc))
            logfile.write(stdout)
            logfile.write(stderr)
            sys.exit(1)

        logfile.write("Checking out %s to %s\n" % (url, workdir))
        rc, stdout, stderr = run_remote_command(host, "%s co %s %s"
                                                % (FCM, url, workdir))
        if rc is 0:
            logfile.write(stdout)
        else:
            logfile.write("Checkout failed: RC=%s\n" % (rc))
            logfile.write(stdout)
            logfile.write(stderr)
            sys.exit(1)

        # Shutdown any instance of the suite if it is still running.
        logfile.write("Shutdown any suite named %s\n" % (opts.name))
        rc, stdout, stderr = run_command(("%s suite-shutdown" +
                                         "--non-interactive --name=%s")
                                         % (ROSE, opts.name))

        # Set-up the rose stem command and add opts if provided.
        stemcmd = ((". /etc/profile; %s %s stem -v -v --new --source=%s " +
                   "--name=%s") % (vertxt, ROSE, workdir, opts.name))

        if opts.keep:
            logfile.write("Running with housekeeping OFF\n")
            stemcmd = "%s -S HOUSEKEEPING=false" % (stemcmd)

        if opts.xcs:
            logfile.write("Running on the xcs\n")
            logfile.write("Overriding METO_HPC_GROUP='xc' to xcs\n")
            rosesuiteconf = os.path.join(workdir, 'rose-stem',
                                         'rose-suite.conf')
            meto_hpc(host, rosesuiteconf, 'xcs')

        if opts.version:
            logfile.write("Overriding FCM_VERSION to next\n")
            rosesuiteconf = os.path.join(workdir, 'rose-stem',
                                         'rose-suite.conf')
            fcm_version_next(host, rosesuiteconf, 'next')

        for source in opts.source[1:]:
            stemcmd = "%s --source=%s" % (stemcmd, source)

        for group in opts.group:
            stemcmd = "%s --task=%s" % (stemcmd, group)

        logfile.write("Stem command is '%s'\n" % (stemcmd))

        # Run the suite with the --new option to remove old files/builds.
        rc, stdout, stderr = run_remote_command(host, stemcmd)
        logfile.write(stdout)
        if rc is not 0:
            logfile.write(stderr)
            sys.exit(1)
