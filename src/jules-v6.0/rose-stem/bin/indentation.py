#!/usr/bin/env python
#
# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************
'''
Module containing various functions used to apply UMDP3 style indentation
to Fortran source code
'''

import re
import sys

# Number of spaces of indent to apply per indentation level
INDENT = 2

# Patterns which match any line which indicates the following line should
# be indented by the above amount
INDENTATION_START = [
    # DO statement - possible label prefix followed by "DO"
    "(^\s*|^\s*\w+\s*:\s*)DO(\s+|$)",
    # SELECT statement - possible label prefix followed by "SELECT CASE"
    "(^\s*|^\s*\w+\s*:\s*)SELECT\s+(CASE|TYPE)(\s*|\()",
    # CASE statement 
    "^\s*CASE(\s*|\()",
    # CLASS IS statement
    "^\s*CLASS\s*IS(\s*|\()",
    # IF statement - must end in "THEN" with possible label prefix
    "^\s*(|\w+\s*:|[0-9]*)\s*IF\s*\(.*?\)\s*THEN\s*$",
    # ELSE statement - as above but allow "ELSE IF"
    "^\s*ELSE\s*(IF|$)",
    # TYPE definition statement - must be followed by a word but can 
    # contain optional arguments
    "^\s*TYPE\s*((,\s*\w+(|\s*\(\s*\w+\s*\)\s*)\s*)*\s*::|)\s*\w+\s*$",
    # WHERE statement - only if missing assignment statement
    "^\s*WHERE\s*\([^\(]*?\)\s*$",
    # ELSEWHERE statement
    "^\s*ELSE\s*WHERE\s*(\(.*?\)[^\w]*$|$)",
]

# Patterns which match any line which signifies that it and subsequent
# lines should be dedented by the above amount
INDENTATION_END = [
    # END DO statement - possibly followed by label suffix
    "^\s*([0-9]*\s*)END\s*DO(\s+\w+\s*$|\s*$)",
    # END SELECT statement - possibly followed by label suffix
    "^\s*END\s*SELECT(\s+\w+\s*$|\s*$)",
    # CASE statement 
    "^\s*CASE(\s*|\()",
    # CLASS IS statement
    "^\s*CLASS\s*IS(\s*|\()",
    # ELSE statement - as above, counts as beggining and end
    "^\s*ELSE\s*(IF|$)",
    # END IF statement - possibly followed by label suffix
    "^\s*END\s*IF(\s+\w+\s*$|\s*$)",
    # END TYPE statement
    "^\s*END\s*TYPE(\s+\w+|)\s*$",
    # END WHERE statement
    "^\s*END\s*WHERE\s*$",
    # ELSEWHERE statement
    "^\s*ELSE\s*WHERE\s*(\(.*?\)[^\w]*$|$)",
]

def simplify_line(lines):
    "A pre-processor for the lines to make them easier to handle"

    # Note we will be passed a slice to include the lines after the
    # current line to the end of the file to allow handling of 
    # continutations
    iline = 0
    line = lines[iline]

    # Strings - try to handle the possibilty that a single-quote may appear 
    # as an apostrophe in a double-quoted string
    line = re.sub("\".*?\"", "", line)
    line = re.sub("'.*?'", "", line)

    # Blank out trailing comments
    if "!" in line:
        line = line[:line.index("!")]

    # Pull any continuation lines into this line
    while re.search(".*&\s*(!.*|)$", line, flags=re.IGNORECASE):
        # Blank out possible trailing comments
        if "!" in line:
            line = line[:line.index("!")]
        # If this results in an empty line we are done
        if line.strip() == "":
            return line
        # Skip following lines if they contain only comments or are empty
        iline += 1
        while (re.search("^\s*!.*$", lines[iline], flags=re.IGNORECASE)
               or re.search("^\s*$", lines[iline], flags=re.IGNORECASE)
               or re.search("^#\s*(if|else|end).*$", lines[iline], 
                            flags=re.IGNORECASE)
               or re.search("^#(\sif).*$", lines[iline],
                            flags=re.IGNORECASE)):
            iline += 1
        line = re.sub("&\s*", "", line) + lines[iline]

    # Repeat the substitution for strings (in case any were pulled 
    # into the string above)
    line = re.sub("\".*?\"", "", line)
    line = re.sub("'.*?'", "", line)

    # Comments - blank out any comments on the trailing line too
    if "!" in line:
        line = line[:line.index("!")]

    # Brackets - remove any nested brackets (i.e. leave only top level brackets)
    # this is to aid with pattern matching where brackets are included
    bracket_nest_level = 0
    new_line = ""
    for char in (line):
        if char == "(":
            bracket_nest_level += 1
            if bracket_nest_level > 1:
                new_line += " "
                continue
        if char == ")":
            if bracket_nest_level > 1:
                new_line += " "
                bracket_nest_level -= 1
                continue
            bracket_nest_level -= 1
        new_line += char

    line = new_line
    return line

def get_current_indent(line):
    "Returns the white-space at the start of a given line"
    return re.search("^(?P<space>\s*)([^\s]|$)", line).group("space")

def indent_line(line, indentation):
    "Returns the given line adjusted by the required amount"
    indentation_str = " "*abs(indentation)
    if indentation > 0:
        return indentation_str + line
    elif indentation < 0:
        if re.search("^"+indentation_str, line):
            return re.sub(indentation_str, "", line, count=1)
        else:
            return line.lstrip()
    else:
        return line

def apply_indentation(lines, debug=False):
    "Apply the indentation rules to a list of lines"
    indentation = 0
    relative_indent = 0
    continuation = False
    new_lines = []

    for iline, line in enumerate(lines):

        if debug:
            print "{0:d}: \"{1:s}\"".format(iline, line)

        # # Ignore line if an OMP directive, or an fcm DEPENDS ON comment
        if (re.search("^\s*!\$.*", line, flags=re.IGNORECASE)):
            if debug:
                print "    (OMP comment)"
            new_lines.append(re.sub(r"^\s*!", "!", line))

        # If the entire line is a comment with the comment character
        # in the first position on the line, indent just the comment
        # text to line up with the current indentation level
        elif re.search("^!.*$", line, flags=re.IGNORECASE) and indentation > 0:
            if debug:
                print "    (Comment only line (from zero indent))"

            # Get the current indentation level of the line
            current_indent = get_current_indent(line)

            # Calculate the relative indent (how far the line must be
            # shifted by to meet the desired indentation level)
            relative_indent = indentation - len(current_indent)

            new_lines.append(indent_line(line, relative_indent))

        # If the line contains a pre-processor directive, don't mess
        # with the indentation but preserve the previous indent etc.
        elif re.search("^#\w+", line, flags=re.IGNORECASE):
            if debug:
                print "    (Pre-processor line)"
            new_lines.append(line)

        # If the line is entirely blank, don't mess with the settings
        # for continuation or the indentation value
        elif re.search("^\s*$", line, flags=re.IGNORECASE):
            if debug:
                print "    (Blank line)"
            new_lines.append(line)

        # If the line is only a comment not starting in the first column
        # copy the indentation of the previous line, unless that would 
        # push the comment to a lower indentation level than the current
        # level (in which case re-indent it like a normal line of code)
        elif re.search("^\s*!.*$", line, flags=re.IGNORECASE):
            if debug:
                print "    (Comment only line)"

            # Get the current indentation level of the line
            current_indent = get_current_indent(line)

            # Check to see if the current relative indent will push this
            # line below the current level (and adjust the relative 
            # indent if necessary)
            if relative_indent + len(current_indent) < indentation:
                relative_indent = indentation - len(current_indent)

            new_lines.append(indent_line(line, relative_indent))

        # If this line is continuing a previous line, instead of
        # trying to calculate its indentation just apply the same
        # indentation as the previous line and move on
        elif continuation:
            new_lines.append(indent_line(line, relative_indent))
            # If the next line is not a continuation reset the flag
            if not re.search(".*&\s*(!.*|)$", line, flags=re.IGNORECASE):
                if debug:
                    print "    (End of continuation)"
                continuation = False

        else:
            # Generate a simplified version of the line for use in
            # pattern matches
            simple_line = simplify_line(lines[iline:])

            if debug:
                print "   \"{0:s}\"".format(simple_line)

            # Check for ending statements first - since the indentation
            # shift for the end of a block must also be applied to the 
            # line containing the block ending statement 
            for pattern in INDENTATION_END:
                if re.search(pattern, simple_line, flags=re.IGNORECASE):
                    if debug:
                        print "    (End, matches {0:s})".format(pattern)
                    indentation -= INDENT

            # Get the current indentation level of the line
            current_indent = get_current_indent(line)

            # Calculate the relative indent (how far the line must be
            # shifted by to meet the desired indentation level)
            relative_indent = indentation - len(current_indent)

            # Indent the line by the required amount and save it to
            # the output array
            new_lines.append(indent_line(line, relative_indent))

            # Now check for starting statements - these will affect
            # the indentation level of future lines (but not the
            # current line)
            for pattern in INDENTATION_START:
                if re.search(pattern, simple_line, flags=re.IGNORECASE):
                    if debug:
                        print "    (Start, matches {0:s})".format(pattern)
                    indentation += INDENT

            # Finally, detect if the following line is a continuation
            # of this one (and therefore requires no indentation)
            if "!" in line:
                line = line[:line.index("!")]
            if re.search(".*&\s*(!.*|)$", line, flags=re.IGNORECASE):
                if debug:
                    print "    (Next line continues)"
                continuation = True

    # Sanity check - indentation should be back to 0 by the end, 
    # if not then error
    if indentation != 0:
        print "Final indentation level non-zero"
        return None

    return new_lines

def main():
    '''Main toplevel function for testing'''
    input_file = sys.argv[1]
    
    with open(input_file, "r+") as file_in:
        lines_in = file_in.read().split("\n")
        new_lines = apply_indentation(lines_in, debug=len(sys.argv) > 2)
        file_in.seek(0)
        file_in.write("\n".join(new_lines))
        file_in.truncate()

if __name__ == "__main__":
    main()
