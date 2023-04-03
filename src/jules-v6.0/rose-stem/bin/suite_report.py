#!/usr/bin/env python
# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************
"""Script to process the results of a suite and write a summary to file. The
   summary is in Trac wiki mark-up. Any projects that do not have a local
   mirror repository are assumed not to be used at that site and are
   excluded from the report.

   Owner: UM System Development Team
   Cylc Suite Syntax: shutdown handler = "suite_report.py -s"
   Command Line syntax:
       suite_report.py -S <suite_dir> [-H] [-s] [-L <log_dir>]
"""

CYLC_SUITE_ENV_FILE = 'cylc-suite-env'
PROCESSED_SUITE_RC = 'suite.rc.processed'
ROSE_SUITE_RUN_CONF = 'rose-suite-run.conf'
SUITE_DB_FILENAME = 'cylc-suite.db'
TRAC_LOG_FILE = 'trac.log'

FCM = {
    'meto': 'fcm',
    'nci': 'fcm',
    'bom' : 'fcm',
    'niwa': 'fcm',
    'vm': 'fcm',
    'jasmin': 'fcm',
    'cehwl1': 'fcm',
    'Unknown': 'true',
}
ROSE_BUSH_URL = {
    'meto': 'http://fcm1/cylc-review',
    'nci': 'http://accessdev.nci.org.au/cylc-review',
    'bom' : 'http://scs-watchdog-dev/rose-bush',
    'niwa':  'http://w-rose01.maui.niwa.co.nz/cylc-review',
    'vm': 'http://localhost/cylc-review',
    'jasmin': 'Unavailable',
    'cehwl1': 'Unavailable',
    'Unknown': 'Unavailable',
}

import glob
import os
import re
import sqlite3
import sys
import traceback
import time
import subprocess
from optparse import OptionParser


def _read_file(filename):
    '''Takes filename (str)
       Return contents of a file, as list of strings.'''
    if os.path.exists(filename):
        with open(filename, 'r') as filehandle:
            lines = filehandle.readlines()
    else:
        print("[ERROR] Unable to find file :\n    \"{0:s}\"".format(filename))
        raise IOError('_read_file got invalid filename : \"{0:s}\"'
                      .format(filename))
    return lines


def _write_file(filename, lines, newline=False):
    '''Takes filemname and list of strings and opt newline boolean.
       Writes array to file of given name.
       The optional newline argument adds a newline at the end of each
       element of the list.
       Returns None'''
    retn = "\n" if newline else ""
    with open(filename, 'w') as filehandle:
        for line in lines:
            filehandle.write("{0:s}{1:s}".format(line, retn))


def _run_command(command, ignore_fail=False):
    '''Takes command and command line options as a list.
       Runs the command with subprocess.Popen.
       Returns the exit code, standard out and standard error as list.
    '''
    pobj = subprocess.Popen(command, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    pobj.wait()
    retcode, stdout, stderr = (
        pobj.returncode, pobj.stdout.read(), pobj.stderr.read())
    if retcode is not 0 and not ignore_fail:
        print("[ERROR] running {0:s}".format(command))
        print("[INFO] RC: {0:}".format(retcode))
        print("[INFO] Stdout: {0:s}".format(stdout))
        print("[INFO] Stderr: {0:s}".format(stderr))
        raise IOError('run_command')
    # Reformat stdout into a list
    stdout = ''.join(stdout)
    stdout = stdout.split("\n")

    return retcode, stdout, stderr


def _remove_quotes(string):
    '''Takes, modifies and returns string.
       Removes all quotes from the string.
       None input results in None output'''
    if string is not None:
        string = re.sub(r'"', r'', string)
        string = re.sub(r"'", r"", string)
    return string


def _dict_merge(main_dict, addon_dict, force=False):
    '''Merge addon dictionary into main dictionary.
       Takes main_dict, addon_dict and optional bool 'force'
       Returns new merged dictionary.
       Optional argument force=True allows forced overwrite of existing
       value with None from the addon dictionary. Otherwise original
       value is preserved when value in addon dict is None.
       This preserving behaviour differentiates it from main.update(addon)'''
    merged_dict = main_dict.copy()
    for key, value in addon_dict.iteritems():
        if isinstance(value, dict):
            if key not in merged_dict:
                merged_dict[key] = {}
            merged_dict[key] = _dict_merge(merged_dict[key], value)
        else:
            if (force         # Switch to Force main to take whatever addon has
                or key not in merged_dict  # No matching key in main - take
                                           # whatever addon has including None
                or value is not None):  # Override main with contents of addon
                merged_dict[key] = value
    return merged_dict


def _select_preferred(option_list):
    '''Takes a list of strings, returns the fist one that is not None.
       If the strings are report text in preffered order it essentially
       ensures you get the preffered option from a list of choices.'''
    pref_opt = None
    for choice in option_list:
        if choice is not None:
            pref_opt = choice
            break
    return pref_opt


def _escape_svn(url):
    '''Takes and returns url as string.
       Escape 'svn:' urls as Trac tries to convert them to links.'''
    if not re.search(r'!svn://', url):  # Make sure it's not already escaped.
        url = re.sub(r'svn://', r'!svn://', url)
    return url


def _get_current_head_revision(mirror_url, fcm_exec):
    '''Given a mirror repository (local) url, uses fcm branch-info to
       retrieve and append the head revision number.
       Requires url and fcm exec path (strings)
       Returns revision number as string'''
    revision = ''
    _, stdout, _ = _run_command([fcm_exec, "branch-info", mirror_url])
    find_last_changed_rev = re.compile(r'Last Changed Rev:\s*(\d+)')
    for line in stdout:
        result = find_last_changed_rev.search(line)
        if result:
            revision = str(result.group(1))
            break
    return revision


def _url_to_trac_link(url):
    '''Takes a URL as string, edits text to resemble a Trac link for code
       on the SRS.
       Returns Trac link form of URL or None if 'svn' was absent from the url.
    '''
    if re.search(r'/svn/', url):
        link_2_url = re.sub(r'svn', r'trac', url)
        elements = link_2_url.split('/')
        elements.insert(elements.index('trac')+2, 'browser')
        link_2_url = '/'.join(elements)
        link_2_url = re.sub(r'@', r'?rev=', link_2_url)
    else:
        link_2_url = None
    return link_2_url


def _parse_string(varname, lines, remove_quotes=True, split_on_comma=False,
                 default_unknown=False):
    """Given a variable name in the rose-suite-run.conf file, return its
    value."""
    find_var = re.compile(r'{0}\s*=\s*(.*)'.format(varname))
    if split_on_comma:
        value = [None]
    elif default_unknown:
        value = 'Unknown'
    else:
        value = None
    for line in lines:
        result = find_var.search(line)
        if result:
            value = result.group(1)
            if remove_quotes:
                value = _remove_quotes(value.rstrip())
            if split_on_comma:
                # Remove brackets and split on comma
                value = re.sub(r'\[', r'', value)
                value = re.sub(r'\]', r'', value)
                value = value.split(',')
    return value       
        

class SuiteReport(object):
    '''Object to hold data and methods required to produce a suite report
       from a rose-stem suite output.'''
    def __init__(self, suite_path, log_path=None, omit_housekeeping=True,
                 sort_by_status=False):
        '''Requires a path to the suite output directory.
           Takes optional arguments for log_path (output dir), and
           booleans omit_housekeeping and sort_by_status to switch on
           omission of houskeeping tasks and sort by status over task name
           when generating the task table in the report.'''
        suite_path = re.sub(r'/$', '', suite_path) # remove trailing '/'
        self.suite_path = suite_path
        self.log_path = log_path
        self.omit_housekeeping = omit_housekeeping
        self.sort_by_status = sort_by_status
        self.creation_time = time.strftime("%Y/%m/%d %X")
        self.uncommitted_changes = 0
        self.suitename = os.path.basename(self.suite_path)

        self.site = 'Unknown'
        self.rose_orig_host = None
        self.suite_owner = None
        self.groups = []
        self.job_sources = {}
        self.projects = {}

        self.suite_owner = self.ascertain_suite_owner()
        (self.site, self.groups, self.trustzone, self.rose, self.cylc, 
         self.fcm) = self.parse_rose_suite_run(self.suite_path)
        self.projects = self.initialise_projects(FCM[self.site])
        self.rose_orig_host, self.job_sources, self.multi_branches = (
                  self.parse_suite_rc_processed(self.suite_path))
        projects = self.check_versions_files()
        self.job_sources = _dict_merge(self.job_sources, projects)

        fcm_exec = FCM[self.site]
        invalid = []
        for project in self.job_sources:
            proj_dict = self.job_sources[project]
            if 'repo loc' in proj_dict:
                proj_dict["repo loc"] = self.convert_to_srs(
                                          proj_dict["repo loc"], self.projects)
            else:
                proj_dict["repo loc"] = self.convert_to_srs(
                                     proj_dict["tested source"], self.projects)

            proj_dict["repo mirror"] = self.convert_to_mirror(
                                          proj_dict["repo loc"], self.projects)
            # If the mirror doesn't exist, move on to the next project.
            if not self.check_repository(fcm_exec, proj_dict["repo mirror"]):
                invalid.append(project)
                continue

            proj_dict["parent mirror"] = self.set_parent(fcm_exec,
                                                     proj_dict["repo mirror"])

            proj_dict["parent loc"] = self.convert_to_srs(
                                     proj_dict["parent mirror"], self.projects)
            # Check "repo loc" and "parent loc" have revisions,
            # and if not, try to get a head of 'trunk' one for them.
            for location in ("repo",  "parent"):
                url = proj_dict[location + " loc"]
                mirror_url = proj_dict[location + " mirror"]
                if url is None or mirror_url is None:
                    continue
                if ":" in url and "@" not in url:
                    revision = _get_current_head_revision(mirror_url, fcm_exec)
                    proj_dict[location + " loc"] = url + '@' + revision
                    proj_dict[location + " mirror"] = (mirror_url + '@' +
                                                             revision)
            proj_dict["repo link"] = self.generate_link(proj_dict["repo loc"])
            proj_dict["parent link"] = self.generate_link(
                                                      proj_dict["parent loc"])
            # If those attempts to generate links didn't work, try the hope
            # and guess approach.
            if proj_dict["repo link"] is None:
                proj_dict["repo link"] = self.link_from_loc_layout(
                    proj_dict["repo link"], proj_dict["repo mirror"], fcm_exec)
            if proj_dict["parent link"] is None:
                proj_dict["parent link"] = self.link_from_loc_layout(
                                               proj_dict["parent loc"],
                                               proj_dict["parent mirror"],
                                               fcm_exec)
            # Final attempt to ensure the links have revision numbers and not
            # keywords which aren't evaluated in the browser.
            if (proj_dict["repo link"] is not None and
                re.search(r'rev=[a-zA-z]', proj_dict["repo link"])):
                revision = self.revision_from_loc_layout(
                                         proj_dict["repo mirror"], fcm_exec)
                proj_dict["repo link"] = re.sub(r'rev=[a-zA-z0-9.]+',
                                     'rev=' + revision, proj_dict["repo link"])
            proj_dict["human repo loc"] = self.convert_to_keyword(
                                      proj_dict["repo loc"], self.projects)
            proj_dict["human parent"] = self.convert_to_keyword(
                                      proj_dict["parent loc"], self.projects)
            proj_dict["ticket no"] = self.ascertain_ticket_number(
                                      proj_dict["repo mirror"], fcm_exec)

        # Finally, remove any projects which were deemed invalid.
        for project in invalid:
            del self.job_sources[project]

    def debug_print_obj(self):
        ''' Debug print method.
            Prints everything in the SuiteReport object.'''
        print("-"*80 + "\nSet up SuiteReport object\n" + "-"*80 + "\n\n")
        for key, value in self.__dict__.iteritems():
            if key == "projects":
                print("{0:s} contains \"{1:d}\" entries.".format(key,
                                                                 len(value)))
            elif key == "omit_housekeeping":
                if value:
                    print("{0:s} is :\"True\"".format(key))
                else:
                    print("{0:s} is :\"False\"".format(key))
            elif key == "job_sources":
                self.print_job_sources(value)
            else:
                print("{0:s} is :\"{1:}\"".format(key, value))
        print("\n" + "-"*80 + "\nEnd of SuiteReport object\n" + "-"*80 + "\n")

    @staticmethod
    def print_job_sources(job_srcs_dict):
        ''' Debug print method.
            Prints everything in projects dictionary.'''
        for key, value in job_srcs_dict.iteritems():
            print("    {0:s} :".format(key))
            for sub_key, sub_value in value.iteritems():
                if isinstance(sub_value, bool):
                    if sub_value:
                        print("        {0:s} is :\"True\"".format(sub_key))
                    else:
                        print("        {0:s} is :\"False\"".format(sub_key))
                else:
                    print("        {0:s} is :\"{1:}\""
                         .format(sub_key, sub_value))

    @staticmethod
    def parse_suite_rc_processed(suite_dir):
        '''Parse the suite.rc.processed file.
           Extract all projects present that begin with a "SOURCE_".
           Allow SOURCE_<project> to override any SOURCE_<project>_<extension>
           entries. Creating a dictionary of format {<project> : <URL>,...}
           Also Extract the host machine rose was launched on.
           Takes full path for suite dir.
           Returns original host as string and a dictionary of projects.'''
        rose_orig_host = 'Unknown rose_orig_host'
        srp_file = os.path.join(suite_dir, PROCESSED_SUITE_RC)
        find_orig_host = re.compile(r'ROSE_ORIG_HOST\s*=\s*(.*)')
        # in pattern below, need to include "_REV" after the project name and
        # before the " *=" and then exclude lines with "_REV" later as
        # otherwise the search will identify PROJ_REV as a unique project
        # name. The other option would be to have an alternate 3rd group match
        # of "_.*?" but that would exclude any future project names that might
        # have an underscore in them.
        find_sources = re.compile(
                     r'\s*(?:HOST_SOURCE_|SOURCE_)(.*?)(|_BASE|_MIRROR|_REV)\s*=\s*(.*)')
        sources = {}
        multiple_branches = {}
        for line in _read_file(srp_file):
            # check for ROSE_ORIG_HOST
            result = find_orig_host.search(line)
            if result:
                rose_orig_host = result.group(1).rstrip()
            # check for SOURCE_.*
            result = find_sources.match(line)
            # Discard the ones which were SOURCE_PROJ_REV
            if result and result.group(2) != "_REV":
                # Allow SOURCE_PROJ to override any existing entries
                # Otherwise only add new entries
                if result.group(1) not in sources or result.group(2) == "":
                    sources[result.group(1)] = {}
                    if ' ' in result.group(3):
                        multiple_branches[(result.group(1))] = result.group(3)
                        sources[result.group(1)]["tested source"] = result.group(3).split()[0]
                    else:
                        sources[result.group(1)]["tested source"] = result.group(3)

        return rose_orig_host, sources, multiple_branches

    @staticmethod
    def ascertain_suite_owner():
        """Find suite owner.
           To retrieve suite owner :
              Take default from ENV['USER'] or 'Unknown Suite Owner' if none.
              Use ENV['CYLC_SUITE_OWNER'] if present.
           Returns suite owner as string."""

        # Get a default
        suite_owner = os.environ.get('USER', 'Unknown Suite Owner')
        suite_owner = os.environ.get('CYLC_SUITE_OWNER', suite_owner)
        return suite_owner

    @staticmethod
    def parse_rose_suite_run(suite_dir):
        """Parse rose-suite-run.conf file.
           Takes full path for suite dir.
           Returns site as string and list of groups rose-stem was run on"""
        rsr_file = os.path.join(suite_dir, 'log', ROSE_SUITE_RUN_CONF)
        lines = _read_file(rsr_file)
        site = _parse_string('SITE', lines, default_unknown=True)
        groups = _parse_string('RUN_NAMES', lines, split_on_comma=True,
                                    remove_quotes=False)
        fcm = _parse_string('FCM_VERSION', lines)
        cylc = _parse_string('CYLC_VERSION', lines)
        rose = _parse_string('ROSE_VERSION', lines)
        trustzone = None
        if 'TRUSTZONE' in os.environ:
            trustzone = os.environ['TRUSTZONE']
        return (site, groups, trustzone, rose, cylc, fcm)

    @staticmethod
    def initialise_projects(fcm_exec):
        '''Uses fcm kp to initialise a directory containing project keywords
           linked to SVN URLS. Format {<project> : <URL>,...}
           Takes full path for suite dir.
           Returns project dictionary.'''

        projects = {}
        _, stdout, _ = _run_command([fcm_exec, "kp"])
        find_primary_loc = re.compile(r'location{primary}')
        find_projects = re.compile(r'\[(.*)\]\s*=\s*(.*)')
        find_x_keyword = re.compile(r'.x$')
        find_xm_keyword = re.compile(r'.xm$')
        find_srs_url = re.compile(r'https://code.metoffice')
        find_mirror_url = re.compile(r'svn:|https://')
        for line in stdout:
            if not find_primary_loc.search(line):
                continue
            result = find_projects.search(line)
            if result:
                project = result.group(1)
                url = result.group(2)
                # Check for keywords conforming to the meto prescribed pattern
                # of ending in '.x' for the external repo and '.xm' for the
                # local mirror.
                if ((find_x_keyword.search(project) and
                     find_srs_url.match(url)) or
                    (find_xm_keyword.search(project) and
                     find_mirror_url.match(url))):
                    projects[project] = url
        return projects

    def check_versions_files(self):
        '''Locate the log/*.version files.
           Call parse_versions_fileto parse the contents of each file.
           Recover which projects are being augmented by branch or WC
           Takes full path for suite dir.
           Returns dictionary of project dictionares and number of projects
           with uncommitted changes.'''
        projects = {}
        self.uncommitted_changes = 0
        find_proj_name = re.compile(r'/(\w+)-\d+.version')

        version_files = glob.glob("{0:s}/*.version".format
                                        (os.path.join(self.suite_path, 'log')))

        for vfile in version_files:
            if 'rose-suite-run.version' in vfile:
                continue
            result = find_proj_name.search(vfile)
            if result:
                project = (result.group(1).upper())
                projects[project] = {}
                (url, revision,
                 wc_changes) = self.parse_versions_file(vfile)
                projects[project]["last changed rev"] = revision
                projects[project]["working copy changes"] = wc_changes
                projects[project]["version file"] = os.path.basename(vfile)
                if wc_changes:
                    self.uncommitted_changes += 1
                if url is not None:
                    if revision is not None:
                        ending = '@' + revision
                    else:
                        ending = ''
                    projects[project]["repo loc"] = url + ending
        return projects

    @staticmethod
    def parse_versions_file(vfile):
        '''Parse a versions file to extract the url and revision for
           the branches behind any working copies, plus any uncommitted
           changes.
           Takes full path to a .version file.
           Returns url and revision as strings plus wc changes as boolean.'''
        url = None
        revision = None
        working_copy_changes = False
        find_svn_status = re.compile(r'SVN STATUS', re.IGNORECASE)
        find_url = re.compile(r'URL:\s*')
        find_last_changed_rev = re.compile(r'Last Changed Rev:\s*')
        for line in _read_file(vfile):
            if find_svn_status.search(line):
                working_copy_changes = True
            if find_url.match(line):
                url = find_url.sub(r'', line).rstrip()
            if find_last_changed_rev.match(line):
                revision = find_last_changed_rev.sub(r'', line).rstrip()
        return url, revision, working_copy_changes

    @staticmethod
    def set_parent(fcm_exec, mirror_url):
        '''For given URL, on the internal mirror repository, use
           'fcm branch-info' to try and ascertain the branch parent, if any.
           Takes fcm_exec path and mirror_url as strings.
           Returns parent URL or None'''

        parent = None
        stdout = ""
        command = [fcm_exec, "branch-info", mirror_url]
        _, stdout, _ = _run_command(command, ignore_fail=True)
        find_branch_parent = re.compile(r'Branch Parent:\s*(.*)')
        for line in stdout:
            result = find_branch_parent.search(line)
            if result:
                parent = result.group(1).rstrip()
        return parent

    @staticmethod
    def check_repository(fcm_exec, url):
        '''Checks whether a given repository is accessible or not.
           Takes fcm_exec path and a url (either SRS or mirror) as strings.
           Returns True if the repository exists, False otherwise.'''
        retcode = 0
        command = [fcm_exec, "info", url]
        retcode, _, _ = _run_command(command, ignore_fail=True)
        if retcode == 0:
            return True
        return False

    @staticmethod
    def generate_task_table(data, omit_housekeeping=True,
                            sort_by_status=False):
        '''Returns a trac-formatted table of the tasks run in this suite.
           Tasks are provided in a dictionary of format {"task" : "status",...}
           omit_houskeeping (bool) removes housekeeping tasks from the
           list when true
           sort_by_status (bool) sorts by status when true, otherwise sorting
           is done by task name'''
        def key_by_name_or_status(task_item):
            '''A key generating function for use by sorted.
               task_item is a tuple of (name, status).
               If sorting by status, return a tuple of (status, name),
               otherwise return (name, status) for use as the sorting key'''
            if sort_by_status:
                return (task_item[1], task_item[0])
            else:
                return task_item

        lines = [" || '''Task''' || '''State''' || "]
        find_housekeep = re.compile(r'housekeep')
        find_monitor = re.compile(r'monitor')
        for task, state in sorted(data.items(), key=key_by_name_or_status):
            if omit_housekeeping and find_housekeep.match(task):
                continue
            if find_monitor.match(task):
                continue
            highlight = "'''"
            if 'succeeded' in state:
                highlight = ""
            lines.append(" || {0:s} || {2:s}{1:s}{2:s} || ".format(task,
                                      state, highlight))
        return lines

    @staticmethod
    def convert_to_mirror(url, projects_dict):
        '''Take a URL as a string, and a dictionary of {project : url, ...}
           If url is a shared repository URL in the projects dictionary convert
           to an internal mirror URL if also available.
           Otherwise return the original URL
           Assumes mirror loc of proj with url svn:something/somewhere/project
           is given as svn:something/somewhere/projectm
        '''
        if url is None:
            return None
        mirror_url = url
        for proj, proj_url in projects_dict.iteritems():
            # checking given url against urls in the projects_dict
            if proj_url in url:
                new_proj = proj + 'm'
                if new_proj in projects_dict:
                    old_proj_url = proj_url
                    new_proj_url = projects_dict[new_proj]
                    mirror_url = re.sub(old_proj_url, new_proj_url,
                                        url, count=1)
                    break
            # checking given url against keywords in the projects_dict
            elif proj in url:
                new_proj = proj + 'm'
                if new_proj in projects_dict:
                    mirror_url = re.sub(proj, new_proj, url, count=1)
                    break
        return mirror_url

    @staticmethod
    def convert_to_srs(url, projects_dict):
        '''Take a URL as a string, and a dictionary of {project : url, ...}
           If url is a mirror repository URL in the projects dictionary convert
           to an SRS URL if also availble.
           Otherwise return the original URL
        '''
        if url is None:
            return None
        srs_url = url
        for proj, proj_url in projects_dict.iteritems():
            # Only check for keywords which correspond to mirror or SRS format
            if re.search(r'.x(|m)$', proj):
                if re.search(proj_url, url):
                # checking given url against urls in the projects_dict
                    shared_project = re.sub(r'm$', r'', proj)
                    if shared_project in projects_dict:
                        mirror_url = proj_url
                        shared_url = projects_dict[shared_project]
                        srs_url = re.sub(mirror_url, shared_url, url, count=1)
                        break
                elif re.search("fcm:" + proj + r'[^m]', url):
                # Looking for an fcm: shorthand notation based on keyword.
                    shared_project = re.sub(r'm$', r'', proj)
                    if shared_project in projects_dict:
                        # if the fcm keyword ends in '_tr' it's on the trunk
                        if re.match(r'fcm:' + proj + r'_tr', url):
                            srs_url = re.sub(r'fcm:' + proj + r'_tr',
                                     projects_dict[shared_project] + r'/trunk',
                                     url, count=1)
                        # if the fcm keyword ends in '_br' it's from branches
                        elif re.match(r'fcm:' + proj + r'_br', url):
                            srs_url = re.sub(r'fcm:' + proj + r'_br',
                                  projects_dict[shared_project] + r'/branches',
                                  url, count=1)
                        # maintain keyword style, but convert to srs.
                        else:
                            srs_url = re.sub(proj, shared_project, url,
                                             count=1)
                        break
        return srs_url

    @staticmethod
    def convert_to_keyword(url, projects):
        '''Takes url and project dictionary.
           Convert a the URL to a keyword based version, if a keyword exists
           in the project dictionary provied.
           Returns None if no keyword is defined.
        '''
        if url is None:
            return None
        keyword_url = None
        for proj, proj_url in projects.iteritems():
            if proj_url in url:
                old_proj_url = projects[proj]
                new_proj_url = "fcm:{0:s}".format(proj)
                keyword_url = re.sub(old_proj_url, new_proj_url, url,
                                     count=1)
                keyword_url = re.sub(r'/trunk', '_tr', keyword_url, count=1)
                keyword_url = re.sub(r'/branches', '_br', keyword_url, count=1)
                break
            if "fcm:" + proj in url:
                keyword_url = url
                break
        return keyword_url

    def generate_link(self, url):
        '''Given a URL, see if it can be made into a shared repository link.
           Returns a link as a str or None'''
        link = None
        if url is not None:
            # Look for a matching part of the URL in the list of projects
            for _, svn in self.projects.iteritems():
                if re.search(svn, url):
                    link = _url_to_trac_link(url)
                    break
        return link

    @staticmethod
    def link_from_loc_layout(url, mirror_url, fcm_exec):
        '''Attempt to generate a link to a url using a bunch of assumptions
           we know mostly hold due to working practices at the Met Office.
           Takes url, mirror url and fcm exec path as strings.
           Returns Link as string or None'''
        link = None
        if url is None or mirror_url is None or re.search(r'^file:/', url):
            return None
        _, stdout, _ = _run_command([fcm_exec, "loc-layout", mirror_url])
        path = None
        root = None
        lproject = None
        revision = None
        find_path = re.compile(r'^path:\s*')
        find_root = re.compile(r'^root:\s*')
        find_project = re.compile(r'^project:\s*')
        find_peg_rev = re.compile(r'^peg_rev:\s*')
        for line in stdout:
            if find_path.match(line):
                path = find_path.sub(r'', line)
                continue
            if find_root.match(line):
                root = find_root.sub(r'', line)
                continue
            if find_project.match(line):
                lproject = find_project.sub(r'', line)
                continue
            if find_peg_rev.match(line):
                revision = find_peg_rev.sub(r'', line)
        if (root is not None and
            lproject is not None
            and path is not None):

            # Convert to a trac url.
            if re.search(r'/svn/', url) :
                url = re.sub(r'svn', r'trac', url)
                elements = url.split('/')
                elements.insert(elements.index('trac')+2, 'browser')
                url = '/'.join(elements)
                if revision is not None:
                    link = url + '?rev={0:s}'.format(revision)
                else:
                    link = url
        return link

    @staticmethod
    def revision_from_loc_layout(mirror_url, fcm_exec):
        '''Attempt to recover a revision number using a url to the mirror
           repository. Also used to translate vn4.3 into 1234'''
        if mirror_url is None:
            return None
        _, stdout, _ = _run_command([fcm_exec, "loc-layout", mirror_url])
        revision = None
        find_peg_rev = re.compile(r'^peg_rev:\s*')
        for line in stdout:
            if find_peg_rev.match(line):
                revision = find_peg_rev.sub(r'', line)
                break
        return revision

    def generate_project_table(self):
        '''Returns a trac-formatted table containing the project source
           trees used in this suite.
           Method of SuiteReport object.
           Returns list of trac formatted table rows.'''

        lines = [" || '''Project''' || '''Tested Source Tree''' || " +
                 "'''Repository Location''' || '''Branch Parent''' || " +
                 "'''Ticket number''' || '''Uncommitted Changes''' ||"]

        def gen_table_element(text_list, link, bold=False):
            '''Takes list of items (strings or Nones) in preference order.
               Calls _select_preferred to get the first non None entry in list.
               Optional Bool "bold" turns on Trac formatting of bold text.
               Formats text as a Trac link if link is not None.
               Returns 'highlighted' text/link and "||" to form a single
               element of a trac table.'''
            text = _select_preferred(text_list)
            highlight = "'''" if bold else ""
            if text is not None and link is not None:
                element = " {2:s}[{0:s} {1:s}]{2:s} || ".format(link, text,
                                                                highlight)
            elif text is not None:
                element = " {1:s}{0:s}{1:s} || ".format(_escape_svn(text),
                                                        highlight)
            else:
                element = " || "
            return element

        for project in sorted(self.job_sources):
            proj_dict = self.job_sources[project]

            line = " || {0:s} || ".format(project)
            line += gen_table_element(
                             [proj_dict["tested source"]],
                             None)
            line += gen_table_element(
                             [proj_dict["human repo loc"],
                              proj_dict["repo loc"]],
                              proj_dict["repo link"])
            line += gen_table_element(
                             [proj_dict["human parent"],
                              proj_dict["parent loc"]],
                              proj_dict["parent link"])
            if "ticket no" in proj_dict:
                if proj_dict["ticket no"] is not None:
                    project_ticket_link = ["{0:s}:{1:s}".format(project,
                                proj_dict["ticket no"])]
                else:
                    project_ticket_link = [None]
            else:
                project_ticket_link = [None]
            line += gen_table_element(project_ticket_link, None)
            wc_link = None
            wc_text = None
            if "working copy changes" in proj_dict:
                if proj_dict["working copy changes"]:
                    wc_text = "YES"
                    wc_link = (r"{0:s}/{1:s}/{2:s}/{3:s}?path=log/{4:s}"
                              .format(ROSE_BUSH_URL[self.site], 'view',
                                      self.suite_owner, self.suitename,
                                      proj_dict["version file"]))
            line += gen_table_element([wc_text], wc_link, bold=True)
            lines.append(line)
        return lines

    @staticmethod
    def get_wallclock_and_memory(filename):
        '''Given a UM output filename read and parse for the wallclock
           and memory.'''
        wallclock = "Unavailable"
        memory = "Unavailable"
        find_wallclock = re.compile(
                           r'PE\s*0\s*Elapsed Wallclock Time:\s*(\d+(\.\d+|))')
        find_total_mem = re.compile(r'Total Mem\s*(\d+)')
        find_um_atmos_exe = re.compile(r'um-atmos.exe')
        check_for_percentage = re.compile('[0-9]+[%]')
        find_mem_n_units = re.compile(
                                  r'(?P<num>[0-9]*\.[0-9]*)(?P<unit>[A-Za-z])')
        try:
            for line in _read_file(filename):
                result = find_wallclock.search(line)
                if result:
                    wallclock = int(round(float(result.group(1))))
                result = find_total_mem.search(line)
                if result:
                    memory = int(result.group(1))
                if find_um_atmos_exe.match(line):
                    split_line = line.split()
                    if check_for_percentage.search(split_line[6]):
                        mem = find_mem_n_units.search(split_line[5])
                        memory = float(mem.group('num'))
                        if mem.group('unit') == 'G':
                            memory *= 1000000
                        elif mem.group('unit') == 'M':
                            memory *= 1000
                    else:
                        memory = int(line.split()[6])
        except Exception as err:
            wallclock = "Failure processing EOJ"
            memory = "Failure processing EOJ"
            stacktr = traceback.format_exc()
            print("[ERROR] Processing wallclock and memory use :\n{0:s}".
                  format(stacktr))
            print("Error type : {0:s}".format(type(err)))
            print(err)
        return wallclock, memory

    @staticmethod
    def query_database(suite_db_file):
        '''Query the database and return a dictionary of states.'''
        database = sqlite3.connect(suite_db_file)
        cursor = database.cursor()
        cursor.execute('select name, status from task_states;')
        data = {}
        for row in cursor:
            data[row[0]] = row[1]
        database.close()
        return data

    @staticmethod
    def ascertain_ticket_number(mirror_url, fcm_exec):
        '''Try and work out the ticket number from the Trac log.
           Takes URL on local (mirror) repository and fcm_exec path.
           Uses 'fcm log'
           Relies on commit line starting with '#[0-9]+' - meto working
           practices for commit says "start with ticket number"
           Returns ticket number as string or None.'''
        ticket_number = None
        if (re.search('/trunk[/@$]', mirror_url) or
            re.search('[fs][cv][mn]:\w+(.xm|.x|)_tr[/@$]', mirror_url)):
            return ticket_number
        _, stdout, _ = _run_command([fcm_exec, "log", "-l", "1", mirror_url])
        for line in stdout:
            result = re.search(r'^\s*(#\d+)', line)
            if result:
                ticket_number = result.group(1)
        return ticket_number

    @staticmethod
    def generate_groups(grouplist):
        '''Convert the list of groups run into a trac-formatted string.'''
        output = ''
        for group in grouplist[:-1]:
            output += "{0:s} [[br]]".format(_remove_quotes(group))
        output += "{0:s}".format(_remove_quotes(grouplist[-1]))
        return output

    def print_report(self):
        ''''Prints a Trac formatted report of the suite_report object'''
        try:
            trac_log = []
            ticket_nos = ""

            # Check to see if any of the soucres have associated tickets and
            # put links to them in the header if so.
            for project, url_dict in self.job_sources.iteritems():
                try:
                    if url_dict["ticket no"] is not None:
                        ticket_nos += ("{0:s}:{1:s} ".format(project,
                                                      url_dict["ticket no"]))
                except KeyError:
                    pass
            if ticket_nos != "":
                trac_log.append(" = Ticket {0:s} ".format(ticket_nos)
                               + "Testing Results - rose-stem output = ")
            else:
                trac_log.append(" = Testing Results - rose-stem output = ")
            trac_log.append("")

            trac_log.append(" || Suite Name: || {0:s} || ".format(
                                                               self.suitename))

            trac_log.append(" || Suite Owner: || {0:s} || ".format(
                                                             self.suite_owner))

            if self.trustzone:
                trac_log.append(" || Trustzone: || {0:s} || ".format(
                                                             self.trustzone))

            if self.fcm:
                trac_log.append(" || FCM version: || {0:s} || ".format(
                                                             self.fcm))

            if self.rose:
                trac_log.append(" || Rose version: || {0:s} || ".format(
                                                             self.rose))

            if self.cylc:
                trac_log.append(" || Cylc version: || {0:s} || ".format(
                                                             self.cylc))

            trac_log.append(" || Report Generated: || {0:s} || "
                            .format(self.creation_time))

            trac_log.append(" || Rose-Bush: || {0:s}/{1:s}/{2:s}/{3:s} || ".
               format(ROSE_BUSH_URL[self.site], 'taskjobs', self.suite_owner,
                      self.suitename))

            trac_log.append(" || Site: || {0:s} || ".format(self.site))
            trac_log.append(" || Groups Run: || {0:s} || ".format(
                            self.generate_groups(self.groups)))
            if self.rose_orig_host is not None:
                trac_log.append(" || ''ROSE_ORIG_HOST:'' || {0:s} || ".format(
                                self.rose_orig_host))
            trac_log.append("")


            if self.uncommitted_changes:
                trac_log.append("")
                trac_log.append("-----")
                word = "changes" if self.uncommitted_changes > 1 else "change"
                trac_log.append(" = WARNING !!! = ")
                trac_log.append("This rose-stem suite included " +
                        "{0:d} uncommitted".format(self.uncommitted_changes) +
                        " project {0:s} and is therefore ".format(word) +
                        "'''not valid''' for review")
                trac_log.append("-----")
                trac_log.append("")


            if len(self.multi_branches.keys()) > 0:
                trac_log.append("")
                trac_log.append("-----")
                trac_log.append(" = WARNING !!! = ")
                trac_log.append("This rose-stem suite included multiple "+ 
                        " branches in {0:d} projects:".format(len(
                                                 self.multi_branches.keys())))
                trac_log.append("")
                for project, branch_names in self.multi_branches.items():
                   trac_log.append("'''{0}'''".format(project))
                   for branch_name in "".join(branch_names).split():
                       trac_log.append(" * {0}".format(branch_name))
                trac_log.append("")
                trac_log.append("-----")
                trac_log.append("")

            trac_log.extend(self.generate_project_table())
            trac_log.append("")

            data = self.query_database(os.path.join(self.suite_path,
                                                    SUITE_DB_FILENAME))
            trac_log.extend(self.generate_task_table(data,
                                                     self.omit_housekeeping,
                                                     self.sort_by_status))
            trac_log.append("")
        except:
            try:
                suite_dir = self.suite_path
            except:
                suite_dir = "--cylc_suite_dir--"
            trac_log.extend(["There has been an exception in " +
                             "SuiteReport.print_report()",
                             "See output for more information",
                             "rose-stem suite output will be in the files :\n",
            "~/cylc-run/{0:s}/log/suite/out and ~/cylc-run/{0:s}/log/suite/err"
                                                 .format(suite_dir)])
        finally:
            # Pick up user specified log path if available,
            # otherwise default to cyclc suite dir.
            trac_log_path = "/No/Path/Provided"
            if self.log_path:
                trac_log_path = os.path.join(self.log_path, TRAC_LOG_FILE)
            else:
                trac_log_path = os.path.join(self.suite_path, TRAC_LOG_FILE)

            # Attempt to provide user with some output,
            # even in event of serious exceptions
            try:
                _write_file(trac_log_path, trac_log, newline=True)
            except IOError:
                print("[ERROR] Writing to {0:s} file : {1:s}".
                       format(TRAC_LOG_FILE, trac_log_path))
                print("{0:s} to this point ".format(TRAC_LOG_FILE) +
                      "would have read as follows :\n")
                print("----- Start of {0:s}.log -----".format(TRAC_LOG_FILE))
                for line in trac_log:
                    print(line)
                print("\n----- End of {0:s}.log -----\n\n"
                      .format(TRAC_LOG_FILE))
                raise
#==============================================================================
#    End of   "class.SuiteReport()"
#==============================================================================


def parse_options():
    '''Use OptionParser to parse options from command line.
       Would prefer ArgParse but currently unavailable in default site
       Pyhton installation at meto.'''
    expected_args = 0
    parser = OptionParser()
    parser.add_option("--include-housekeeping", "-H", dest="include_housekeep",
                      action="store_true", help="Include housekeeping tasks")
    parser.add_option("--status-sort", "-s", dest="sort_by_status",
                      action="store_true",
                      help="Sort task table by task status")
    parser.add_option("--suite-path", "-S", dest="suite_path",
                      action="store", help="Path to suite")
    parser.add_option("--log_path", "-L", dest="log_path", action="store",
                      default=None,
                      help="Output dir for {0:s}".format(TRAC_LOG_FILE))

    (opts, args) = parser.parse_args()

    # Find the suite database file
    if opts.suite_path is None:
        if 'CYLC_SUITE_RUN_DIR' in os.environ:
            opts.suite_path = os.environ['CYLC_SUITE_RUN_DIR']
            # cylc provides 3 args to the shutdown handler - so allow for them.
            expected_args = 3
        else:
            sys.exit("Path to suite not provided.")

    if len(args) != expected_args:
        parser.print_help()
        message = ("got {0:d} extra arguments, expected {1:d}\n"
                  .format(len(args), expected_args))
        message += "Extra arguments are :\n{0:}".format(args)
        print message
        sys.exit(message)
    return opts


def main():
    '''Main program.
       Sets up a SuiteReport object and calls it's print_report method.'''
    opts = parse_options()

    suite_report_obj = SuiteReport(suite_path=opts.suite_path,
                            log_path=opts.log_path,
                            omit_housekeeping=(not opts.include_housekeep),
                            sort_by_status=opts.sort_by_status)

    suite_report_obj.print_report()

if __name__ == '__main__':
    main()
