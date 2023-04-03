#! /usr/bin/env python
# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************
'''
This module contains various functions for applying UMDP3 styling to
Fortran source code for JULES
'''
import re
import sys
from indentation import simplify_line

CODE_REPLACEMENTS = [
    # Replace Fortran 77 style conditional keywords
    (r'\.eq\.', ' == '), 
    (r'\.ne\.', ' /= '), 
    (r'\.gt\.', ' >  '),
    (r'\.lt\.', ' <  '), 
    (r'\.ge\.', ' >= '), 
    (r'\.le\.', ' <= '),
    # Ensure remaining comparitive logicals have spaces either side
    (r'([^\s])(?<!\()\.not\.', r'\g<1> .not.'),
    (r'\.not\.([^\s])', r'.not. \g<1>'),
    (r'([^\s])\.and\.', r'\g<1> .and.'),
    (r'\.and\.([^\s])', r'.and. \g<1>'),
    (r'([^\s])\.or\.', r'\g<1> .or.'),
    (r'\.or\.([^\s])', r'.or. \g<1>'),
    (r'([^\s])\.eqv\.', r'\g<1> .eqv.'),
    (r'\.eqv\.([^\s])', r'.eqv. \g<1>'),
    (r'([^\s])\.neqv\.', r'\g<1> .neqv.'),
    (r'\.neqv\.([^\s])', r'.neqv. \g<1>'),
    # Ensure hard-coded real numbers have a zero after the decimal point
    (r'([0-9])\.([^0-9]|$)', r'\g<1>.0\g<2>'),
    # Remove start of line ampersands, without changing the spacing of
    # the line they appear on
    (r'^(\s*)&', r'\g<1> '),
    # Make constructs which include brackets have exactly one space
    # between the construct and the bracket character
    (r'(^\s*)(\w+\s*:\s*|[0-9]+\s*|)((else|)\s*if)(|\s\s+)\(', 
      '\g<1>\g<2>\g<3> ('),
    (r'(^\s*)(\w+\s*:\s*|[0-9]+\s*|)where(|\s\s+)\(', '\g<1>\g<2>where ('),
    (r'(^\s*)case(|\s\s+)\(', '\g<1>case ('),
    (r'\)(|\s\s+)then(\W|$)', ') then\g<1>'),
    # Make block-ending keywords have exactly one space between them
    (r'(^\s*[0-9]*\s*)else(|\s\s+)if', '\g<1>else if'),
    (r'(^\s*[0-9]*\s*)end(|\s\s+)do', '\g<1>end do'),
    (r'(^\s*[0-9]*\s*)end(|\s\s+)forall', '\g<1>end forall'),
    (r'(^\s*[0-9]*\s*)end(|\s\s+)function', '\g<1>end function'),
    (r'(^\s*[0-9]*\s*)end(|\s\s+)interface', '\g<1>end interface'),
    (r'(^\s*[0-9]*\s*)end(|\s\s+)module', '\g<1>end module'),
    (r'(^\s*[0-9]*\s*)end(|\s\s+)program', '\g<1>end program'),
    (r'(^\s*[0-9]*\s*)end(|\s\s+)select', '\g<1>end select'),
    (r'(^\s*[0-9]*\s*)end(|\s\s+)subroutine', '\g<1>end subroutine'),
    (r'(^\s*[0-9]*\s*)end(|\s\s+)type', '\g<1>end type'),
    (r'(^\s*[0-9]*\s*)end(|\s\s+)where', '\g<1>end where'),
    (r'(^\s*[0-9]*\s*)end(|\s\s+)if', '\g<1>end if'),
    (r'(^\s*)(\w+\s*:\s*|[0-9]+\s*|)select(|\s\s+)case', 
      '\g<1>\g<2>select case'),
    (r'(^\s*)(\w+\s*:\s*|[0-9]+\s*|)select(|\s\s+)type', 
      '\g<1>\g<2>select type'),
    (r'(^|\W)go(|\s\s+)to(\W|$)', '\g<1>go to\g<3>'),
    # Make intent statements contain no extra spacing inside the brackets
    (r'(.*intent\s*)\(\s*in\s*\)(.*)', r'\g<1>(in)\g<2>'),
    (r'(.*intent\s*)\(\s*out\s*\)(.*)', r'\g<1>(out)\g<2>'),
    (r'(.*intent\s*)\(\s*inout\s*\)(.*)', r'\g<1>(inout)\g<2>'),
    # Make module USE, ONLY statments have exactly no space between the ONLY
    # and the colon character after it
    (r'^(\s*)use(\s*,\s*\w+\s*::|)(\s+\w+\s*,\s*)only\s*:(.*)$',
     r'\g<1>use\g<2>\g<3>only:\g<4>'),
    # Replace tabs with 1 space no matter what the tabs represented, indentation
    # to be run after tabs replaced
    (r'\t', ' '),
    # Spaces around math operators
    (r'(\w|\))(-|\+)(\w|\()', r'\g<1> \g<2> \g<3>'),
    (r'(\de)\s(-|\+)\s(\d)', r'\g<1>\g<2>\g<3>'),
    (r'(\))(-|\+)(\s*)', r'\g<1> \g<2> '),

    (r'(\S)={2}(\S)', r'\g<1> == \g<2>'),

    (r'(\w)([^\*])\*(\w|\D)', r'\g<1>\g<2> * \g<3>'),
    (r'(\*\s\*)', r'**'),

    (r'(\w)([^=])=(\w)', r'\g<1>\g<2> = \g<3>'),

    (r'(=)(\s)(\*)(\s)', r'\g<1>\g<3>'),

    (r'(\w)(>|<)\s=\s', r'\g<1> \g<2>= '),

    (r'(\w|\))([^/])/(\w|\()', r'\g<1>\g<2> / \g<3>'),
    (r'(\))/(\()', r'\g<1> / \g<2>'),

    (r'(\d)(\/)(\w)', r'\g<1> / \g<3>'),
    (r'(\w)(\/)(\()', r'\g<1> / \g<3>'),

    (r'(\d)\*(\(|\w)', r'\g<1> * \g<2>'),
    (r'(\))\*(\d)', r'\g<1> * \g<2>'),

    (r'(,\S*)(\s)(-|\+)(\s)(\d)', r'\g<1>\g<3>\g<5>'),

    (r'(\w|\))\*(\w|\()', r'\g<1> * \g<2>'),

    (r'(\(\S*)(\s)\+(\s)(\d)', r'\g<1>+\g<4>'),
    (r'(\(\S*)(\s)\-(\s)(\d)', r'\g<1>-\g<4>'),
    (r'(\S*)(\s)\+(\s)(\d)(\))', r'\g<1>+\g<4>\g<5>'),
    (r'(\S*)(\s)\-(\s)(\d)(\))', r'\g<1>-\g<4>\g<5>'),

    (r'(\))(\-|\+)(\d)', r'\g<1> \g<2> \g<3>'), 

    (r'(\s\s)\*(\s\s)', r' * '),

    (r'(\)|\d)(\s)(\*\*)', r'\g<1>\g<3>'),
    (r'(\s)(\*\*)', r'\g<2>'),
    (r'(\w)(\s)(\+)(\s)(\d)(,)', r'\g<1>+\g<5>,'),

    (r'(\w)(>|<)(\w)', r'\g<1> \g<2> \g<3>'),
    (r'(\s)<(\s)(commentmarker)' ,r'<\g<3>'),

    (r'(,)(\s)(\*)(\s)', r'\g<1>\g<3>'),

    (r'(\w)=(\.)', r'\g<1> = \g<2>'),

    (r'(CHARACTER)(\s)(\*)(\s)', r'\g<1>\g<3>'),
    (r'([DO\s\w*])=(\w)', r'\g<1> = \g<2>'),

    (r'(\s+)FILE(\s*)= ', r'\g<1>FILE='),
    (r'(KIND)(\s)=(\s)', r'\g<1>='),
    (r'(\()(LEN)(\s)=(\s)', r'\g<1>\g<2>='),

    (r'(\))(\/)', r'\g<1> \g<2>'),
    (r'(\()(\.)', r'\g<1> \g<2>'),
    (r'(CHARACTER)(\s)(\()(LEN)', r'\g<1>\g<3>\g<4>'),

    (r'(DATA)(\s)(\S*)(/)(\S*)(/)', r'\g<1> \g<3> \g<4> \g<5> \g<6>'),
    (r'(\d)(/)', r'\g<1> \g<2>'),

    (r'(\W\s+)\*\*', r'\g<1> **'),
    (r'(\)\s+)\*\*',r'\g<1>**'),
    (r'(\w)(\s+)(\s)\*(\s)(\s+)(&)', r'\g<1>\g<2>  *  \g<5>\g<6>'),
]

COMMENT_REPLACEMENTS = [
    # DEPENDS ON fcm constructions
    (r'^(\s*!)\s*depends\s*on\s*:\s*',r'\g<1> DEPENDS ON: '),
]

KEYWORDS = set([
    "abs",
    "abstract",
    "access",
    "achar",
    "acos",
    "action",
    "adjustl",
    "adjustr",
    "advance",
    "aimag",
    "aint",
    "algama",
    "all",
    "allocatable",
    "allocate",
    "allocated",
    "alog",
    "alog10",
    "amax0",
    "amax1",
    "amin0",
    "amin1",
    "amod",
    "and",
    "anint",
    "any",
    "asin",
    "associated",
    "atan",
    "atan2",
    "backspace",
    "bind",
    "bit_size",
    "blank",
    "btest",
    "c_alert",
    "c_associated",
    "c_backspace",
    "c_bool",
    "c_carriage_return",
    "c_char",
    "c_double",
    "c_double_complex",
    "c_f_pointer",
    "c_f_procpointer",
    "c_float",
    "c_float128",
    "c_float128_complex",
    "c_float_complex",
    "c_form_feed",
    "c_funloc",
    "c_funptr",
    "c_horizontal_tab",
    "c_int",
    "c_int128_t",
    "c_int16_t",
    "c_int32_t",
    "c_int64_t",
    "c_int8_t",
    "c_int_fast128_t",
    "c_int_fast16_t",
    "c_int_fast32_t",
    "c_int_fast64_t",
    "c_int_fast8_t",
    "c_int_least128_t",
    "c_int_least16_t",
    "c_int_least32_t",
    "c_int_least64_t",
    "c_int_least8_t",
    "c_intmax_t",
    "c_intptr_t",
    "c_loc",
    "c_long",
    "c_long_double",
    "c_long_double_complex",
    "c_long_long",
    "c_new_line",
    "c_null_char",
    "c_null_funptr",
    "c_null_ptr",
    "c_ptr",
    "c_ptrdiff_t",
    "c_short",
    "c_signed_char",
    "c_size_t",
    "c_sizeof",
    "c_vertical_tab",
    "cabs",
    "call",
    "case",
    "ccos",
    "cdabs",
    "cdcos",
    "cdexp",
    "cdlog",
    "cdsin",
    "cdsqrt",
    "ceiling",
    "cexp",
    "char",
    "character",
    "class",
    "clog",
    "close",
    "cmplx",
    "common",
    "complex",
    "conjg",
    "contains",
    "continue",
    "cos",
    "cosh",
    "count",
    "cpu_time",
    "cqabs",
    "cqcos",
    "cqexp",
    "cqlog",
    "cqsin",
    "cqsqrt",
    "cshift",
    "csin",
    "csqrt",
    "cycle",
    "dabs",
    "dacos",
    "dasin",
    "data",
    "datan",
    "datan2",
    "date_and_time",
    "dble",
    "dcmplx",
    "dconjg",
    "dcos",
    "dcosh",
    "ddim",
    "deallocate",
    "default",
    "deferred",
    "delim",
    "derf",
    "derfc",
    "dexp",
    "dfloat",
    "dgamma",
    "digits",
    "dimag",
    "dimension",
    "dint",
    "direct",
    "dlgama",
    "dlog",
    "dlog10",
    "dmax1",
    "dmin1",
    "dmod",
    "dnint",
    "do",
    "dot_product",
    "dprod",
    "dsign",
    "dsin",
    "dsinh",
    "dsqrt",
    "dtan",
    "dtanh",
    "elemental",
    "else",
    "elsewhere",
    "end",
    "endfile",
    "eor",
    "eoshift",
    "epsilon",
    "equivalence",
    "eqv",
    "erf",
    "erfc",
    "exist",
    "exit",
    "exp",
    "exponent",
    "extends",
    "external",
    "false",
    "file",
    "final",
    "float",
    "floor",
    "fmt",
    "forall",
    "form",
    "format",
    "formatted",
    "fraction",
    "function",
    "gamma",
    "go",
    "huge",
    "iabs",
    "iachar",
    "iand",
    "ibclr",
    "ibits",
    "ibset",
    "ichar",
    "idim",
    "idint",
    "idnint",
    "ieor",
    "if",
    "ifix",
    "implicit",
    "import",
    "in",
    "include",
    "index",
    "inout",
    "inquire",
    "int",
    "integer",
    "intent",
    "interface",
    "intrinsic",
    "iomsg",
    "ior",
    "iostat",
    "iqint",
    "is",
    "ishft",
    "ishftc",
    "isign",
    "iso_c_binding",
    "kind",
    "lbound",
    "len",
    "len_trim",
    "log",
    "log10",
    "logical",
    "matmul",
    "max",
    "max0",
    "max1",
    "maxexponent",
    "maxloc",
    "maxval",
    "merge",
    "min",
    "min0",
    "min1",
    "minexponent",
    "minloc",
    "minval",
    "mod",
    "module",
    "modulo",
    "move_alloc",
    "mvbits",
    "name",
    "named",
    "namelist",
    "nearest",
    "neqv",
    "new_line",
    "nextrec",
    "nint",
    "nml",
    "non_intrinsic",
    "none",
    "nopass",
    "not",
    "null",
    "nullify",
    "number",
    "only",
    "open",
    "opened",
    "optional",
    "or",
    "out",
    "pack",
    "pad",
    "parameter",
    "pass",
    "pointer",
    "position",
    "precision",
    "present",
    "print",
    "private",
    "procedure",
    "protected",
    "product",
    "program",
    "public",
    "pure",
    "qabs",
    "qacos",
    "qasin",
    "qatan",
    "qatan2",
    "qcmplx",
    "qconjg",
    "qcos",
    "qcosh",
    "qdim",
    "qerf",
    "qerfc",
    "qexp",
    "qgamma",
    "qimag",
    "qlgama",
    "qlog",
    "qlog10",
    "qmax1",
    "qmin1",
    "qmod",
    "qnint",
    "qsign",
    "qsin",
    "qsinh",
    "qsqrt",
    "qtan",
    "qtanh",
    "radix",
    "random_number",
    "random_seed",
    "range",
    "read",
    "readwrite",
    "real",
    "rec",
    "recl",
    "recursive",
    "repeat",
    "reshape",
    "result",
    "return",
    "rewind",
    "rrspacing",
    "save",
    "scale",
    "scan",
    "select",
    "selected_int_kind",
    "selected_real_kind",
    "sequence",
    "sequential",
    "set_exponent",
    "shape",
    "sign",
    "sin",
    "sinh",
    "size",
    "sngl",
    "spacing",
    "spread",
    "sqrt",
    "status",
    "stop",
    "subroutine",
    "sum",
    "system_clock",
    "tan",
    "tanh",
    "target",
    "then",
    "tiny",
    "to",
    "transfer",
    "transpose",
    "trim",
    "true",
    "type",
    "ubound",
    "unformatted",
    "unit",
    "unpack",
    "use",
    "value",
    "verif",
    "where",
    "while",
    "write",
    ])


def pre_process_line(line):
    ''' Returns a line of code with any strings or comments replaced by markers
    so that they don't confuse the processing.

    Any single quoted strings will be replaced with <single##marker>,
    where ## indicates the index in the array singlequotes returned.

    Any double quoted strinigss will be replaced with <double##marker>,
    where ## indicates the index in the array doublequotes returned.

    Comments will be replaced with <commentmarker>, The comment is returned
    in comment

    Once returned line is processed the Quoted strings and comments can be
    replaced by calling ReplaceStringComments'''

    # 1) For single quotes 
    singlequotes = None
    singlequotes = re.findall("'.*?'", line)
    if len(singlequotes) > 0:
        for i in range(len(singlequotes)):
            line = line.replace(
                singlequotes[i], '<single{0:d}marker>'.format(i),1)


    # 2) For double quotes 
    doublequotes = None
    doublequotes = re.findall("\".*?\"" , line)
    if len(doublequotes) > 0:
        for i in range(len(doublequotes)):
            line = line.replace(
                doublequotes[i],'<double{0:d}marker>'.format(i),1)


    # 3) See if any comments are left 
    comment = None
    found_comments = re.search("!.*", line)
    if found_comments != None:
        comment = found_comments.group(0)

        # Replace the Commented Text with a marker
        line = line.replace(comment,'<commentmarker>')


    return line, singlequotes, doublequotes, comment


def post_process_line(
    line, singlequotes, doublequotes, comment):
    ''' Returns a line of code with single/double quote markers where 
    the markers have been replaced by given strings. Used to put back strings 
    in code after processing. Intended to be used with pre_process_line.'''

    # 1) Now put text back in reverse order
    if comment != None:
        # Replace the commented text marker with comment text
        line = line.replace('<commentmarker>', comment)


    # 2) For double quotes 
    if doublequotes != None:
        for i in range(len(doublequotes)):
            line = line.replace(
                '<double{0:d}marker>'.format(i),doublequotes[i],1)


    # 3) For single quotes 
    if singlequotes != None:
        for i in range(len(singlequotes)):
            line = line.replace(
                '<single{0:d}marker>'.format(i),singlequotes[i],1)

    return line


def replace_patterns(line):
    '''Replace various patterns according to the styling guidelines on
    the provided line, returning the result'''

    stripline = line.strip()

    if len(stripline) == 0 or stripline[0] == "#":
        return line
    
    line, singlequotes, doublequotes, comment = pre_process_line(line)

    # Replace patterns found in code
    for pattern, replacement in CODE_REPLACEMENTS:
        if re.search(pattern, line, flags=re.IGNORECASE):
            line = re.sub(r"(?i)"+pattern, replacement, line)

    # Replace patterns found in comments
    if comment and line.strip() == "<commentmarker>" :
        for pattern, replacement in COMMENT_REPLACEMENTS:
            if re.search(pattern, comment, flags=re.IGNORECASE):
                comment = re.sub(r"(?i)"+pattern, replacement, comment)

    line = post_process_line(line, singlequotes, doublequotes, comment)

    return line


def upcase_keywords(line):
    '''Upper-case any Fortran keywords on the given line, and down-case any
    all capital words which aren't keywords, returning the result'''

    stripline = line.strip()

    if len(stripline) == 0 or stripline[0] == "!" or stripline[0] == "#":
        return line

    line, singlequotes, doublequotes, comment = pre_process_line(line)

    # Split the line of code into a set of words
    line_words = set(re.findall(r"[\w]+", line))

    for word in line_words:
        # Exclude ifdefs and special "__FILE__" or "__LINE__" directives
        if (word.isupper() and 
            not re.search("defined\(\s*{0:s}\s*\)".format(word), line) and
            not re.match("__\w+__", word)):
            for _ in range(line.count(word)):
                line = re.sub(r'(^|\W){0:s}(\W|$)'.format(word), 
                              r'\g<1>{0:s}\g<2>'.format(word.lower()), 
                              line)

    line_words = set([word.lower() for word in line_words])
    words_to_upcase = list(line_words.intersection(KEYWORDS))
    for keyword in words_to_upcase:
        for _ in range(line.lower().count(keyword)):
            line = re.sub(r'(?i)(^|\W){0:s}(\W|$)'.format(keyword), 
                          r'\g<1>{0:s}\g<2>'.format(keyword.upper()), 
                            line)

    # Now add back any comments/strings
    line = post_process_line(line, singlequotes, doublequotes, comment)

    return line

def declaration_double_colon(iline, lines):
    '''Attempt to add the double colon to definition lines which do not already
    have it'''

    line = lines[iline]

    for declare_type in ["REAL", "LOGICAL", "CHARACTER", "INTEGER"]:
        if re.search("^\s*{0:s}\W".format(declare_type), 
                     line, flags=re.IGNORECASE):
            # Pre-process the line to pull in any continuation lines
            simple_line = simplify_line(lines[iline:])
            if not re.search("\s+FUNCTION(,|\s|\()", 
                             simple_line, flags=re.IGNORECASE):
                # The presence of declaration attributes (ALLOCATABLE,
                # PUBLIC, POINTER, etc) are only valid when used with
                # the double-colon.  Therefore after allowing for the 
                # presence of either a (KIND/LEN=...) statement or an
                # older "*INT" declaration the first character should
                # not be a comma
                search = re.search(
                    "^(\s*{0:s}\s*(\(.*?\)|\*\s*[0-9]+|))\s+\w+".format(
                        declare_type), simple_line, flags=re.IGNORECASE)
                if search:
                    # Group 1 contains everything up to the start of the 
                    # variable definition
                    statement = search.group(1).strip()
                    # Attempt to fit the double-colon into an existing space to
                    # preserve indentation, otherwise just add it to the line
                    line = re.sub(r"{0:s}\s\s\s\s".format(re.escape(statement)),
                                  r"{0:s} :: ".format(statement), line, count=1)
                    line = re.sub(r"{0:s}\s*((?<!\s\s\s\s)\w)".format(
                            re.escape(statement)), r"{0:s} :: \g<1>".format(
                            statement), line, count=1)

    return line


def apply_styling(lines):
    '''For a lot of lines apply UMDP3 styling to each line and return 
    the result'''

    output_lines = []
    for iline, line in enumerate(lines):
        line = declaration_double_colon(iline, lines)
        line = replace_patterns(line)
        line = upcase_keywords(line)
        output_lines.append(line)

    return output_lines

def main():
    '''Main toplevel function for testing'''
    input_file = sys.argv[-1]
    with open(input_file, "r+") as file_in:
        lines_in = file_in.read().split("\n")
        new_lines = apply_styling(lines_in)
        file_in.seek(0)
        file_in.write("\n".join(new_lines))
        file_in.truncate()
    

if __name__ == '__main__':
    main()
