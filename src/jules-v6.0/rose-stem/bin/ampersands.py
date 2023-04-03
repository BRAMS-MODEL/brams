#!/usr/bin/env python
# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************
'''
This module contains various functions for putting continuation ampersands in
the same column throughout FORTRAN source code for JULES.

If there are comments after the continuation ampersand, then the whitespace is
adjusted as much as possible to reduce line to the required length.

Lines that are still too long will be identified, optionally with a message in
stdout.

Note that this code has made an effort to deal with the possibility of
amersands and exclamation marks within quoted strings or comments, but there
may be some cases which are missed. These lines will be left without applying
the ampersand shifting, and will be flagged, optionally with a message in
stdout.
'''
import sys
import re
from optparse import OptionParser

# Default line length for JULES.
DEFAULT_COL = 79


class ParsingError(Exception):
    '''
    Raised when an operation attempts to parse a line that doesn't look
    like it expects.
    '''
    def __init__(self):
        self.msg = 'Parsing Error.'
    pass


class CharError(ParsingError):
    '''
    Raised when there are an unexpected number of a certain char in a line.
    '''
    def __init__(self, char, number):
        self.number = number
        self.char = char
        self.msg = "There are {0:d} unquoted, uncommented " \
                   "\"{1:s}\" in this line.".format(number, char)
    pass


class QuoteError(ParsingError):
    '''
    Raised when there are an unexpected number of quote marks (single or
    double) in a line.
    '''
    def __init__(self, quote, number):
        self.number = number
        self.quote = quote
        self.quotemarks = {"'": "single quotes",
                           "\"": "double quotes"}

        self.msg = "There are an odd number ({0:d}) of {1:s} in " \
                   "this line.".format(number, self.quotemarks[quote])
    pass


def print_message(errtype, msg, iline=None, line=None, fname=None):
    '''
    Print a formatted message
    '''
    if fname is None:
        fnamestr = ""
    else:
        fnamestr = "{0:s}:".format(fname)

    if iline is None:
        if fnamestr is None:
            ilinestr = ""
        else:
            ilinestr = " "
    else:
        ilinestr = "L{0:05d}: ".format(iline)

    if line is None:
        linestr = ""
    else:
        linestr = ": {0:s}".format(line)

    print "{0:s}{1:s}{2:s} - {3:s}{4:s}".format(fnamestr, ilinestr, errtype,
                                                msg, linestr)


def find_quoted_char(line, char):
    '''
    Check if a particular string (char) is present inside double or single
    quotes within a line of text, and return locations of any occurrences
    within the line.

    e.g. for a line containing a quote,
    > line = "She said 'I like tea at Fortnum & Mason'. I prefer fish & chips."
                                            ^
    The function

    > locs = find_quoted_char(line, "&")

    will return the location of the first ampersand

    > print locs
    [32,]

    as it is within the single quoted speech. The location of the second
    ampersand is not returned, as it is not within a quote in the line.

    Note that this assumes that the line has been pre-processed to remove
    apostrophes. It will fail if there are an odd number of single or
    double quotes on the line.
    '''

    # First check there are any instances of char in this line. If not, just
    # leave without doing anything.
    if re.search(char, line) is None:
        return None

    # Then check if there are any quote marks. If not, just leave without doing
    # anything
    if re.search("['\"]", line) is None:
        return None

    # Now check if the search string contains single or double quotes. If so,
    # then can't search for it within those quotes, so just consider the
    # other type of quotes in the line. If the search string does not contain
    # a quote mark, then can search for it within single or double quotes.
    if re.search("'", char) is not None:
        # If char contains single quotes, then it can't exist inside single
        # quotes, so only check for double quoted text.
        dodouble = True
        dosingle = False
    elif re.search("\"", char) is not None:
        # If char contains double quotes, then it can't exist inside double
        # quotes, so only check for single quoted text.
        dodouble = False
        dosingle = True
    else:
        # If char contains no kind of quotes, then check if the line
        # contains any kind of quoted text.
        dodouble = True
        dosingle = True

    # Check for the existence of char within single quoted text.
    if dosingle:
        # If there are an odd number of single quotes, then something has gone
        # wrong with the parsing, or a corner case has been found.
        # Raise an error and stop processing this line.
        if len(re.findall("'", line)) % 2 != 0:
            raise QuoteError("'", len(re.findall("'", line)))

        # Find pairs of single quotes. This is assuming that if there are
        # two single quotes in a line, the text between is quoted.
        # This may give unexpected results if there are actually two
        # apostrophes.
        singlequotes = [singlequote
                        for singlequote in re.finditer("'.*?'", line)]

        # Find the location of any instances of the char in the quoted text.
        char_loc_sq = [loc.start() + singlequote.start()
                       for singlequote in singlequotes
                       for loc in re.finditer(char, singlequote.group())]
    else:
        # Not getting single quoted text, so return a zero-length list.
        char_loc_sq = []

    # Check for the existence of char within double quoted text.
    if dodouble:
        # If there are an odd number of double quotes, then something has gone
        # wrong with the parsing, or a corner case has been found.
        # Raise an error and stop processing this line.
        if len(re.findall("\"", line)) % 2 != 0:
            raise QuoteError("\"", len(re.findall("\"", line)))

        # Find pairs of double quotes. This is assuming that if there are
        # two double quotes in a line, the text between is quoted.
        # This may give unexpected results if this isn't really the case.
        doublequotes = [doublequote
                        for doublequote in re.finditer("\".*?\"", line)]

        # Find the location of any instances of the char in the quoted text.
        char_loc_dq = [loc.start() + doublequote.start()
                       for doublequote in doublequotes
                       for loc in re.finditer(char, doublequote.group())]
    else:
        # Not getting double quoted text, so return a zero-length list.
        char_loc_dq = []

    # Add both together and return.
    char_loc = char_loc_sq + char_loc_dq
    if len(char_loc) > 0:
        return char_loc
    else:
        return None


def find_commented_char(line, char):
    '''
    Check if a particular string (char) is present inside a Fortran comment
    within a line of text, and return locations of any occurrences within the
    line.

    e.g. for a line of Fortran code, say we want to find the opening
    parenthesis in the comment
    > line = "INTEGER, INTENT(IN) :: land_pts    ! No. land points (to run)"
                                                             ^
    The function

    > locs = find_commented_char(line, "\(")

    will return the location of the parenthesis in the comment

    > print locs
    [53,]

    as it is within the commented text. The location of the parenthesis
    around IN is not returned, as it before the quote starts.

    Note that this will give unexpected results if there are bangs that aren't
    comment markers. It assumes that the line has been pre-processed with
    find_quoted_char to remove any bangs that are within strings.
    '''

    # First check there are any instances of char, if not just leave
    # without processing the line.
    if re.search(char, line) is None:
        return None

    # Now check if there are comments in this line. Note that this assumes that
    # any bang present is a comment indicator. If there are any others (e.g.
    # within quoted text), then this may fail. It is recommended to pre-process
    # the line to make sure that the first bang on the line indicates the start
    # of a comment. In normal usage, run find_quoted_char before
    # find_commented_char to ensure it isn't a problem.
    if re.search("!", line) is None:
        # If there are no bangs, we can go.
        return None

    # Find comment location, but only if there is actually something after the
    # bang. If not, then there is no actual comment.
    # Note that this assumes that the first bang that is found is the start of
    # a comment (it assumes that any within quotes have already been dealt
    # with).
    comment = re.search("!.+", line)

    # If there is a comment, find any instances of the char in the comment.
    if comment is None:
        char_loc = []
    else:
        char_loc = [loc.start() + comment.start()
                    for loc in re.finditer(char, comment.group())]

    # Return the location of any commented characters.
    if len(char_loc) > 0:
        return char_loc
    else:
        return None


def replace_characters(line, locs, lens, replchar="X"):
    '''
    Replace characters in a line at particular locations.

    e.g. for the line
    > line = "She said 'I like tea at Fortnum & Mason'. I prefer fish & chips."

    And locations and lengths
    > locs = [32, ]
    > lens = [1, ]

    The result is
    > newline = replace_characters(line, locs, lens, replchar="+"):
    > print newline
    She said 'I like tea at Fortnum + Mason'. I prefer fish & chips.

    '''

    # This code requires that the replacement is a single string, which is
    # used as many times as necessary.
    if len(replchar) != 1:
        raise ValueError("Argument 'replchar' must be a single character")

    # Split the line into a list, as characters can't be changed in place in
    # a string.
    newline = list(line)

    # Loop through locations and lengths and replace each character with the
    # replacement character character.
    for loc, ln in zip(*[locs, lens]):
        for l in range(ln):
            newline[loc+l] = replchar

    # Return the newline, joined back together as a string.
    return "".join(newline)


def shift_ampersand(line, col=DEFAULT_COL):
    '''
    Check if the line contains an ampersand.
    If so then set location of ampersand to col so as to be consistent
    Sometimes there are comments after the ampersand, in this case keep
    the comment and ampersand where they are, unless the line is too long, in
    which case reduce any whitespace between the ampersand and the
    comment. If the line is still too long, then reduce whitespace between
    the end of the code and the ampersand until the comment fits within the
    required line length.
    '''

    # Pre-processor lines start with #. Ignore them completely.
    pre_proc = re.match("#", line.strip())

    # Lines that are completely commented start with a bang and are also
    # ignored completely.
    all_comment = re.match("!", line.strip())

    # Ignore empty lines, pre-processor directives and whole-line comments.
    if len(line.strip()) > 0 and pre_proc is None and all_comment is None:

        # First need to make sure that the continuation character and any
        # comments are found correctly. Can't just search for the
        # first ampersand or bang in the line as it is possible that strings
        # (text inside single or double quotes) and comments could contain
        # ampersands or bangs, and these must be ignored when processing the
        # line.
        # The solution is to first find any possibly misleading characters
        # and replace them with temporary substitutes. Then find the location
        # of any continuation character and/or comment in the line, before
        # replacing the substituted characters, and the process the line to
        # correctly align it.
        # The order of substitutions is
        # 1. Find any apostrophes that are within double quotes, and replace
        #    them with "W", so that any single quote marks that remain are
        #    paired quotes.
        #    (This will actually also find paired single quotes inside double
        #    quotes, but that is OK as any characters that need to be
        #    substituted inside the single quotes will also be found inside
        #    the outer double quotes later).
        #    e.g.
        #    CALL log_info("init_veg", "Doell's method", 'something else')
        #    becomes
        #    CALL log_info("init_veg", "DoellWs method", 'something else')
        # 2. Find any ampersands that are within pairs of double or single
        #    quotes and replace them with "X".
        #    e.g.
        #    CALL log_info("init_veg", "Doell & Siebert",   &
        #    becomes
        #    CALL log_info("init_veg", "Doell X Siebert",   &
        # 3. Find any bangs that are within pairs of double or single quotes
        #    and replace them with "Y".
        #    e.g.
        #    CALL log_info("init_veg", "Doell method!")  ! Describe method
        #    becomes
        #    CALL log_info("init_veg", "Doell methodY")  ! Describe method
        # 4. If 2. or 3. fail because of non-matching single quotes, then it
        #    could be that there is an apostrophe in a comment. Try to find
        #    it by finding the location of the comment and then replacing any
        #    commented apostrophes with "w".
        # 5. Find any ampersands that are within comments and replace them
        #    with "Z".
        #    e.g.
        #    CALL log_info("init_veg", "Doell-Siebert",   & ! Doell & Siebert
        #    becomes
        #    CALL log_info("init_veg", "Doell-Siebert",   & ! Doell Z Siebert
        # Once these replacements have been done, the line should only contain
        # at most one ampersand (as continuation character), and if it contains
        # any bangs, the first one is the start of the comment. Then find
        # the location of these to start any required processing.

        # Find where there are apostrophes inside double quotes and replace
        # them with "W" temporarily.
        try:
            quoted_apos_locs = find_quoted_char(line, "'")
        except:
            raise
        if quoted_apos_locs is not None:
            lens = len(quoted_apos_locs)*[1, ]
            line = replace_characters(line, quoted_apos_locs, lens,
                                      replchar="W")

        # Find where there are ampersands or bangs within single or double
        # quotes. If these searches fail, check in case there are apostrophes
        # after the quote mark.
        try:
            quoted_amp_locs = find_quoted_char(line, "&")
            quoted_bang_locs = find_quoted_char(line, "!")
            commented_apos_locs = None
        except QuoteError as e:
            # If either of these have failed because there are an odd number of
            # single quotes then it is possible that it is because of an
            # apostrophe in a comment (at this stage any apostrophes within
            # double quotes will have been found and dealt with). If they
            # failed because of an odd number of double quotes, there is
            # something wrong so just re-raise the error.
            if e.quote == "\"":
                raise
            else:
                # Find paired quote marks.
                singlequotes = [singlequote
                                for singlequote in re.finditer("'.*?'", line)]
                doublequotes = [doublequote
                                for doublequote in re.finditer("\".*?\"",
                                                               line)]
                # Loop through any bangs in the line and check if they are
                # between paired quotes.
                quoted_bangs = []
                for bang in re.finditer("!", line):
                    quoted = False
                    for singlequote in singlequotes:
                        if bang.start() >= singlequote.start() and \
                                bang.start() < singlequote.end():
                                    quoted = True
                                    break
                    for doublequote in doublequotes:
                        if bang.start() >= doublequote.start() and \
                                bang.start() < doublequote.end():
                                    quoted = True
                                    break
                    quoted_bangs.append(quoted)

                # Assume first unquoted bang is the comment.
                try:
                    comment_bang_loc = quoted_bangs.index(False)

                # If we can't find an unquoted bang, assume there is still
                # something wrong, so re-raise the quote error.
                except:
                    raise QuoteError(e.quote, e.number)

                # Otherwise, find the commented apostrophes, and replace
                # with "w" temporarily.
                commented_apos_locs = [apostrophe.start() + comment_bang_loc
                                       for apostrophe in
                                       re.finditer("'",
                                                   line[comment_bang_loc:])]
                if commented_apos_locs is not None:
                    lens = len(commented_apos_locs) * [1, ]
                    line = replace_characters(line, commented_apos_locs, lens,
                                              replchar="w")

                # Re-try finding quoted characters (if it still fails, will
                # raise another exception to be passed up).
                try:
                    quoted_amp_locs = find_quoted_char(line, "&")
                except:
                    raise
                try:
                    quoted_bang_locs = find_quoted_char(line, "!")
                except:
                    raise

        # If any ampersands or bangs were found within quotes, replace them
        # with "X" and "Y" respectively.
        if quoted_amp_locs is not None:
            lens = len(quoted_amp_locs) * [1, ]
            line = replace_characters(line, quoted_amp_locs, lens,
                                      replchar="X")
        if quoted_bang_locs is not None:
            lens = len(quoted_bang_locs) * [1, ]
            line = replace_characters(line, quoted_bang_locs, lens,
                                      replchar="Y")

        # Find where there are ampersands within comments and replace them with
        # "Z" temporarily.
        commented_amp_locs = find_commented_char(line, "&")
        if commented_amp_locs is not None:
            lens = len(commented_amp_locs) * [1, ]
            line = replace_characters(line, commented_amp_locs, lens,
                                      replchar="Z")

        # Check if there is still more than one ampersand in this line and
        # warn if there is.
        if (len(re.findall("&", line)) > 1):
            raise CharError("&", len(re.findall("&", line)))

        # Find the first ampersand in the line (if there is one).
        amp_loc = line.find("&")

        # Find the first bang in the line to see if it is a comment.
        comment_loc = line.find("!")

        # Now the locations of the characters that are needed have been found,
        # replace ampersands and bangs where they were.
        if quoted_apos_locs is not None:
            lens = len(quoted_apos_locs) * [1, ]
            line = replace_characters(line, quoted_apos_locs, lens,
                                      replchar="'")
        if quoted_amp_locs is not None:
            lens = len(quoted_amp_locs) * [1, ]
            line = replace_characters(line, quoted_amp_locs, lens,
                                      replchar="&")
        if quoted_bang_locs is not None:
            lens = len(quoted_bang_locs) * [1, ]
            line = replace_characters(line, quoted_bang_locs, lens,
                                      replchar="!")
        if commented_amp_locs is not None:
            lens = len(commented_amp_locs) * [1, ]
            line = replace_characters(line, commented_amp_locs, lens,
                                      replchar="&")
        if commented_apos_locs is not None:
            lens = len(commented_apos_locs) * [1, ]
            line = replace_characters(line, commented_apos_locs, lens,
                                      replchar="'")

        # If the line contains an unquoted or uncommented ampersand, need
        # to check if it is in the right place.
        if amp_loc != -1:

            # If there is no comment, then shift the ampersand (and remove any
            # trailing whitespace.
            if comment_loc == -1:

                # Put the ampersand at the end of an empty line.
                templine = list(" ") * col
                templine[-1] = "&"

                # Keep the part of the input line before the ampersand (without
                # white space).
                line = line[:amp_loc].rstrip()
                temp_length = len(line)

                # Don't require a space between the end of the code and the
                # continuation character, so maximum length of the code is
                # col - 1.
                if temp_length > col - 1:
                    # If the line is too long, just add the ampersand at the
                    # end (after a space).
                    line = list(line)
                    line.append(" &")
                    line = "".join(line)
                else:
                    # Otherwise, add the text to the empty line with the
                    # ampersand at the end.
                    line = list(line)
                    templine[0:temp_length] = line[0:temp_length]
                    line = "".join(templine)

            # If there is a comment after the ampersand then only do anything
            # if the line is too long, in which case can try to take white
            # space to get it down to length.
            elif comment_loc > amp_loc:

                # If the line is too long, first make sure there is only one
                # space between ampersand and comment, and get rid of any
                # trailing whitespace.
                if len(line.rstrip()) > col:
                    comment = line[comment_loc:].rstrip()
                    line = line[:amp_loc+1]
                    line = " ".join([line, comment])

                # If the line is still too long, see if ampersand can be moved
                # left to reduce the length sufficiently.
                if len(line.rstrip()) > col:

                    # How much do we need to reduce the length to get it short
                    # enough?
                    nloop = len(line) - col

                    # Get the current location of the ampersand.
                    amp_location = amp_loc

                    # Cast the line to a list as characters can't be edited in
                    # place in a string.
                    line = list(line)

                    # Try cutting down the whitespace one step at a time until
                    # there is only one space left.
                    for i in range(nloop):

                        if line[amp_location-1] == " ":
                            # If there is still whitespace that can be removed,
                            # remove it and update the ampersand location.
                            del line[amp_location-1]
                            amp_location -= 1
                        else:
                            # Ampersand is now next to no-blank text so place
                            # a single space and exit the loop.
                            line.insert(amp_location, " ")
                            break

                    # Join the line back up into a string.
                    line = "".join(line)

    return line


def check_line_len(line, maxlinelen=DEFAULT_COL, fname=None,
                   debug=False):
    '''
    Check line to see if it violates length requirements. If debugging,
    write some information to stdout.
    '''

    line_too_long = len(line) > maxlinelen

    return line_too_long


def apply_ampersand_shift(lines, col=DEFAULT_COL, fname=None, debug=False):
    '''
    For a lot of lines make sure any continuation ampersands are in the
    same column and return the result
    '''

    output_lines = []
    not_parsed = []
    for iline, line in enumerate(lines):
        try:
            outline = shift_ampersand(line, col)
        except ParsingError as e:
            if debug:
                print_message("PARSING ERROR",
                              "{0:s} Ampersand shifting has not been "
                              "applied".format(e.msg), iline+1, line=line,
                              fname=fname)
            outline = line
            not_parsed.append(iline)
        output_lines.append(outline)

    return output_lines, not_parsed


def apply_check_line_len(lines, fname=None, maxlinelen=DEFAULT_COL,
                         debug=False):
    '''
    For a lot of lines check if any lines are longer than required
    '''

    any_too_long = False
    ilines_too_long = []
    for iline, line in enumerate(lines):
        iline_too_long = check_line_len(line, maxlinelen, fname, debug)
        if iline_too_long:
            any_too_long = True
            ilines_too_long.append(iline)
            if debug:
                print_message("VIOLATION",
                              "Line > {0:d} columns".format(maxlinelen),
                              iline+1, line=line, fname=fname)

    if any_too_long:
        return ilines_too_long
    else:
        return None


def main():
    '''
    Main toplevel function for testing
    '''
    parser = OptionParser(usage="""
    %prog [--column col] [--debug] file_1 [file_2 [file_3] ... ]

    This script will attempt to manipulate white space to make sure ampersands
    are consistently in the same column in each line.

    If the line is still too long, it will minimise its length.

    The optional --column tells which column should be used (default=79)
    """)
    parser.add_option("--column", dest="col", type="int", default=DEFAULT_COL,
                      help="Column in which ampersands should appear")
    parser.add_option("--debug", action="store_true",
                      help="Report useful information")

    (opts, args) = parser.parse_args()

    input_file = args[0]

    with open(input_file, "r+") as file_in:
        lines_in = file_in.read().split("\n")
        new_lines, not_parsed = apply_ampersand_shift(lines_in, opts.col,
                                                      fname=input_file,
                                                      debug=opts.debug)
        if opts.debug:
            if len(not_parsed) > 0:
                print_message("WARNING",
                              "Ampersand alignment failed for some lines "
                              "due to parsing errors. Please check lines and "
                              "make sure they are correct.",
                              fname=input_file)
            ilines_too_long = apply_check_line_len(new_lines, input_file,
                                                   maxlinelen=opts.col,
                                                   debug=True)
            if ilines_too_long is not None:
                print_message("WARNING",
                              "Some lines are longer than {0:d} characters. "
                              "Please check and make them "
                              "shorter.".format(opts.col), fname=input_file)

        file_in.seek(0)
        file_in.write("\n".join(new_lines))
        file_in.truncate()

if __name__ == "__main__":
    main()
