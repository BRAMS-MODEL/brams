#!/usr/bin/env python
#
# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************

'''
Rose-stem test for checking the branch working copy with the trunk.

Fail if any files are changed by the umdp3_checker.py script.

Return a list of any files changed by the umdp3_checker.py script.

Usage:
 Fortran files must end in .f90, .F90 or .inc extension.
 Copy files to work/tmp.
 Run the script on the whole tmp dir.
 Diff with the working copy.
 Fail if changes and return files that need changing.

'''

from optparse import OptionParser
import os
import shutil
import subprocess
import tempfile


def copy_working_branch(jules_source):
    '''Copy the working version of the branch to a tmp dir.'''
    # Make the tmp dir in the cwd which is the work/ of the task name.
    tmpdir = tempfile.mkdtemp()
    tmp_filename = 'jules_diff'

    # Ensure the file is read/write by the creator only
    saved_umask = os.umask(0077)

    tmp_path = os.path.join(tmpdir, tmp_filename)

    # Copy the src/ from the working copy to the tmp dir.

    src_wkcopy = os.path.join(jules_source, "src")
    try:
        shutil.copytree(src_wkcopy, tmp_path)
    # Directories are the same
    except shutil.Error as err:
        print('Directory not copied. Error: %s' % err)
    # Any error saying that the directory doesn't exist
    except OSError as err:
        print('Directory not copied. Error: %s' % err)
    return tmp_path, saved_umask, tmpdir


def diff_cwd_working(jules_source, path, saved_umask, tmpdir):
    '''Diff the tmp dir with the working branch and report diff.'''
    diff = subprocess.Popen("diff -qr " + jules_source + "/src " + path,
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            shell=True)
    _ = diff.wait()

    if diff.returncode == 0:
        print("[OK] No changes were made by the UMDP3 checker script and "
              "the working copy complies with the coding standards. "
              "No action required.")
    else:
        diff_stdout = diff.stdout.read()
        diff_files = diff_stdout.strip().split("\n")
        print("[FAIL] The following files were changed when the "
              "umdp3_fixer.py script was run:")
        for diff_filesname in diff_files:
            # diffs are of the form "Files <x> and <y> differ"
            # we select only <x>
            print("[FAIL] " + diff_filesname.split()[1])
        print("Please run rose-stem/bin/umdp3_fixer.py on each of the "
              "failed files in your working copy and check the changes. "
              "Then commit the changes to your branch and then re-run all "
              "rose-stem testing.")
        os.umask(saved_umask)
        shutil.rmtree(tmpdir)
        raise ValueError("Ran fcm diff command and changes were made by " +
                         "the umdp3_fixer.py script.")

    return


def run_umdp3checker(jules_source, path):
    '''Run the umdp3 fixer script in the tmp dir copy of the working branch.'''
    umdp3 = subprocess.Popen(jules_source + "/rose-stem/bin/umdp3_fixer.py " +
                             "$(find " + path +
                             " -name '*.[F|f]90' -o -name '*.inc' | xargs)",
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             shell=True)
    rtncode = umdp3.wait()
    if rtncode != 0:
        print "[FAIL] Problem while attempting to run umdp3_fixer.py"
        raise ValueError("Problem while attempting to run umdp3_fixer.py")
    return


def main():
    '''Take in the location of working branch location and run the tests.'''
    # Initialise the command line parser.
    description = 'Args for the source code...'
    parser = OptionParser(description=description)
    parser.add_option('--source',
                      dest='source',
                      action='store',
                      help='source of the JULES branch',
                      default='None')
    # e.g. "--source jules_source_branch"

    # Parse the command line.
    (opts, _) = parser.parse_args()
    jules_source = opts.source

    (path, saved_umask, tmpdir) = copy_working_branch(jules_source)
    run_umdp3checker(jules_source, path)
    diff_cwd_working(jules_source, path, saved_umask, tmpdir)

    os.umask(saved_umask)
    shutil.rmtree(tmpdir)

if __name__ == "__main__":
    main()
