#!/usr/bin/env python
#
# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************
'''
Top level module for the UMDP3 fixer / code styling tool

Usage:
 To apply UMDP3 styling to a specific file or set of files:

   umdp3_fixer.py <filename1> [<filename2> <filename3> ...

 Or to apply UMDP3 styling to files showing differences to
 the JULES trunk so they match the UM trunk coding style:

   umdp3_fixer.py --branch-diff

 Fortran files must end in .f90, .F90 or .inc extension for the 
 branch-diff version to apply to them
'''

import os
import re
import sys
import subprocess
from optparse import OptionParser
from indentation import apply_indentation
from styling import apply_styling
from ampersands import apply_ampersand_shift, print_message

def get_branch_diff():
    '''If in a local working copy of an FCM branch, return a 
    list of files returned by running a branch-diff command'''

    # Use the bdiff command to extract a list of files that have changed
    # on the user's branch
    bdiff = subprocess.Popen("fcm bdiff --summarize",
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             shell=True)
    _ = bdiff.wait()
    if bdiff.returncode == 0:
        bdiff_stdout = bdiff.stdout.read() 
        bdiff_files = bdiff_stdout.strip().split("\n")
        bdiff_files = [bfile[8:] for bfile in bdiff_files
                       if bfile[0] == "M" or bfile[0] == "A"]

    if len(bdiff_files) == 1 and bdiff_files[0] == "":
        raise ValueError("Unable to run fcm bdiff command")

    # Use the fcm info command to extract the root name of the repository
    info = subprocess.Popen("fcm info",
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             shell=True)
    _ = info.wait()
    repos_root = ""
    if info.returncode == 0:
        info_out = info.stdout.read().strip().split("\n")
        for line in info_out:
            search = re.match("^Repository Root: (?P<url>.*)", line)
            if search:
                repos_root = search.group("url")
    if repos_root == "":
        raise ValueError("Unable to run fcm info command")

    # Deal with external/internal repository differences and strip the output
    # of the bdiff command to become a relative path to the file
    if re.match("svn://fcm\d/.*", repos_root):
        bdiff_files = [os.path.relpath(bfile, 
                                       os.path.join(repos_root, "JULES", "trunk"))
                       for bfile in bdiff_files]
    elif re.match("https://code.metoffice.gov.uk/.*", repos_root):
        bdiff_files = [os.path.relpath(bfile, 
                                       os.path.join(repos_root, "main", "trunk"))
                       for bfile in bdiff_files]
    else:
        raise ValueError("Unrecognised Repository root: {0:s}"
                         .format(repos_root))

    return bdiff_files

def main():
    '''Main toplevel function'''
    parser = OptionParser(usage="""
    %prog [--branch-diff] [file_1 [file_2] [file_3] ...]
    
    This script will attempt to apply UMDP3 conformant styling to a Fortran
    source file or a set of source files.  It accepts an unlimited number of
    filenames as arguments and will apply styling to each file in turn.

    NOTE: It will overwrite the contents of the files so it is recommended to
    run only on working copies which do not have uncommitted changes - making
    it easy to revert should the results be undesired.

    The optional --branch-diff flag will instead assume your current directory
    is within a working copy and apply the styling only to files listed by the
    \"fcm branch-diff\" command.
    """)
    parser.add_option("--branch-diff", dest="bdiff", 
                      action="store_true", help="Run on a branch diff")
    parser.add_option("--max-line-len", dest="maxlinelen", default=79,
                      type="int", help="Maximum length of lines")
    (opts, args) = parser.parse_args()

    if opts.bdiff:
        if len(args) > 0:
            sys.exit("ERROR: Cannot specify filenames and --branch-diff")
        files = get_branch_diff()
    else:
        files = args

    for input_file in files:
        print "Processing: {0:s}".format(input_file)
        if (input_file.split(".")[-1] != "F90" and
            input_file.split(".")[-1] != "f90" and
            input_file.split(".")[-1] != "inc"):
            print ("Input file {0:s} not a "
                   "Fortran file, skipping".format(input_file))
            continue

        with open(input_file, "r+") as file_in:
            lines = file_in.read().split("\n")
            styled_lines = apply_styling(lines)
            indented_lines = apply_indentation(styled_lines)
            file_in.seek(0)
            if indented_lines:
                amp_lines, not_parsed = apply_ampersand_shift(indented_lines,
                                                        opts.maxlinelen)
                if len(not_parsed) > 0:
                    print_message("WARNING",
                                  "Ampersand alignment failed for some lines "
                                  "due to parsing errors. Please check lines "
                                  "and make sure they are correct.",
                                  fname=input_file)

                file_in.write("\n".join(amp_lines))
            else:
                print_message("WARNING",
                              "Indentation failed. Ampersand alignment not "
                              "applied",
                              fname=input_file)
                file_in.write("\n".join(styled_lines))
            file_in.truncate()


if __name__ == "__main__":
    main()
