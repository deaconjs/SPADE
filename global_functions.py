import string
import os
def translate_filename_to_DOS_8_3_format(filename):
    """ Assumes working on Windows. clustalw doesn't seem to like opening files in normal
    format when a directory in the path has a space in it. Change the format to 8.3 format,
    where each filename has at most 8 chars in the filename (using ~1 format)."""
    filename = os.path.abspath(filename)
    components = string.split(filename, '\\')
    out = ""
    for i in range(len(components)):
        stuff = string.split(components[i], ' ')
        components[i] = string.join(stuff, '')
    for i in range(len(components)-1):
        if len(components[i]) > 8:
            components[i] = components[i][:6] + '~1'
        out = out + components[i] + '\\'
    x = string.rfind(components[-1],'.')
    if x == -1 or x<8:
        new_string = components[-1]
    else:       # if there's a dot and the base is > 8 characters
        new_string = components[-1][:6] + '~1' + components[-1][x:]
    out = out + new_string
    return out

# given command line arguments, returns a dict with keys ie. '-i' or '--help' and values
def getopts(argv):
    opts = {}
    while argv:
        if argv[0][0] == '-':
            opts[argv[0]] = argv[1]
            argv = argv[2:]
        else:
            argv = argv[1:]
    return opts
