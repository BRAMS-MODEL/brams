#!/usr/bin/env python
# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************
'''
This module is part of the new release process. It copies the HEAD metadata
into a vnX.Y folder and consistently renames the url fields from "latest".

Usage:
   create_jules_version_metadata.py  <current version> <new version>

The version numbers can be specified with or without "vn".
'''
import os
import sys
import re


def read_file(fname):
    '''Return contents of a file given filename.'''
    with open(fname, 'r') as fh:
        lines = fh.readlines()
    return lines


def write_file(fname, lines, newline=False):
    '''Write a file given names and contents. The optional newline argument
    adds a newline at the end of each element of the list.'''
    with open(fname, 'w') as fh:
        for line in lines:
            if newline:
                fh.write("%s\n"%(line))
            else:
                fh.write(line)


def run_command(command, shell=False):
    '''Given a command as a string, run it and return the exit code, standard
    out and standard error. The optional shell argument allows a shell to
    be spawned to allow multiple commands to be run.'''

    import subprocess

    if shell:
        # Create the Popen object and connect out and err to pipes using
        # the shell=True option.
        p = subprocess.Popen(command, stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE, shell=True)
    else:
        # Turn command into a list
        command_list = command.split()

        # Create the Popen object and connect out and err to pipes
        p = subprocess.Popen(command_list, stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE)

    # Do the communicate and wait to get the results of the command
    stdout, stderr = p.communicate()
    rc = p.wait()

    # Reformat stdout
    stdout = ''.join(stdout)
    stdout = stdout.split("\n")

    return rc, stdout, stderr

def copy_metadata(new_version):
    '''Create a vnX.Y metadata. This copies the HEAD metadata to vnX.Y for
     jules-fcm-make and jules-standalone.'''

    # First change any spurious umX.Y or vnX.Y in urls to "latest"
    check_for_spurious_url_tags("../rose-meta/jules-standalone/HEAD/rose-meta.conf")
    check_for_spurious_url_tags("../rose-meta/jules-fcm-make/HEAD/rose-meta.conf")

    # Then copy
    rc, stdout, stderr = run_command("fcm cp ../rose-meta/jules-standalone/HEAD ../rose-meta/jules-standalone/vn%s"%(new_version))
    rc, stdout, stderr = run_command("fcm cp ../rose-meta/jules-fcm-make/HEAD ../rose-meta/jules-fcm-make/vn%s"%(new_version))

def check_for_spurious_url_tags(fname):
    # Change any spurious umX.Y or vnX.Y in urls to "latest"
    lines = read_file(fname)
    newlines = []
    for line in lines:
        if re.search(r'url=', line):
            line = re.sub(r'/um\d+\.\d/', r"/latest/", line)
            line = re.sub(r'/vn\d+\.\d/', r"/latest/", line)
        newlines.append(line)
    write_file(fname, newlines)


def update_jules_version_metadata(fname, new_version):
    '''Change url tags from "latest" to vnX.Y'''

    verstr = "vn%s"%(float(new_version))
    lines = read_file(fname)
    newlines = []
    for line in lines:
        if re.search(r'url=', line):
            line = re.sub(r'/latest/', r"/%s/"%(verstr), line)
        newlines.append(line)
    write_file(fname, newlines)

def check_version_numbers(fname, current_version, new_version):
    '''Check that the command line version numbers are sensible'''

    lines = read_file(fname)
    for line in lines:
        if re.search(r'BEFORE_TAG', line):
            print line.strip()
            if not re.search(r"%s"%(current_version), line):
                print "%s does not match 'Current version'\n"%(fname)
                something = raw_input("Press CTRL-c to abort, or enter to continue: ")
            else:
                print "%s consistent with 'Current version'"%(fname)

    if abs(new_version - current_version - 0.1) > 0.01:
        print ("Increment between 'Current version' and 'New version' greater"
               " than normal value 0.1. Was this intentional?\n")
        something = raw_input("Press CTRL-c to abort, or enter to continue: ")

if __name__ == '__main__':

    # Check that $PWD is rose-stem directory
    cwd = os.getcwd()
    if not re.search(r'rose-stem$', cwd):
        sys.exit("Please run in the rose-stem subdirectory")

    # Get version number information from command line
    if len(sys.argv) < 2:
        sys.exit("Syntax: create_jules_version_metadata.py <current version> <new version>\ne.g. 'create_jules_version_metadata.py 5.7 5.8'")

    for i in range(1, 3):
      if 'vn' in sys.argv[i]:
          sys.argv[i] = re.sub(r'vn', r'', sys.argv[i])

    current_version = float(sys.argv[1])
    new_version = float(sys.argv[2])

    verstr = "%s_%s"%(int(current_version*10), int(new_version*10))

    print "Current version: ",current_version
    print "New version: ",new_version

    # We currently don't need the current version, but it might be used later
    # and serves as a useful check to make sure that the new version number is
    # also correct.
    print "Checking 'Current version' against versions.py BEFORE_TAG"
    check_version_numbers("../rose-meta/jules-standalone/versions.py",
                          current_version, new_version)

    # Copy metadata for jules-standalone and jules-fcm-make
    print "Copying HEAD metadata"
    copy_metadata(new_version)

    print "Changing the JULES url tag from 'latest' to vn%s in ../rose-meta/jules-*/vn%s"%(new_version, new_version)
    update_jules_version_metadata("../rose-meta/jules-standalone/vn%s/rose-meta.conf"%(new_version), new_version)
    update_jules_version_metadata("../rose-meta/jules-fcm-make/vn%s/rose-meta.conf"%(new_version), new_version)

