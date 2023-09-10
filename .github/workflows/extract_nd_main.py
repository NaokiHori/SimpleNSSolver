import os
import sys
import glob
import enum


def get_filenames(root):
    results = glob.glob(f"{root}/**", recursive=True)
    retvals = list()
    for result in results:
        if result.endswith(".c") or result.endswith(".h"):
            retvals.append(result)
    return retvals


class NdimsType(enum.Enum):
    IN_2D = enum.auto()
    IN_3D = enum.auto()
    OTHER = enum.auto()


def extract_given_dim(ndims, lines):
    state = NdimsType.OTHER
    if_level = 0
    if_level_ndims = 0
    newlines = list()
    for line in lines:
        is_on_ndims_macro = False
        if "#if" in line:
            # found "if", increase nest counter
            if_level += 1
            if " NDIMS" in line:
                if "NDIMS==2" in line.replace(" ", ""):
                    is_on_ndims_macro = True
                    # now in 2D condition
                    state = NdimsType.IN_2D
                    if_level_ndims = if_level
                if "NDIMS==3" in line.replace(" ", ""):
                    is_on_ndims_macro = True
                    # now in 3D condition
                    state = NdimsType.IN_3D
                    if_level_ndims = if_level
        elif "#else" in line:
            # check this "else" is for ndims
            if if_level == if_level_ndims:
                is_on_ndims_macro = True
                # if it is, swap state (3d if now 2d, vice versa)
                if state == NdimsType.IN_2D:
                    state = NdimsType.IN_3D
                elif state == NdimsType.IN_3D:
                    state = NdimsType.IN_2D
                else:
                    print("else found but if not found beforehand")
                    sys.exit()
        elif "#endif" in line:
            if if_level == if_level_ndims:
                is_on_ndims_macro = True
                state = NdimsType.OTHER
            # found "endif", reduce nest counter
            if_level -= 1
        if not is_on_ndims_macro:
            # we do not include macro about ndims
            if ndims == 2 and state != NdimsType.IN_3D:
                newlines.append(line)
            if ndims == 3 and state != NdimsType.IN_2D:
                newlines.append(line)
    return newlines


def modify_comments(lines):
    """
        there are weird comments which are used by Sphinx, which look like
          // <comment> | <number of lines><cr>
        I use this function to modify this kind of stuffs as
          // <comment><cr>
    """
    delim = " | "
    newlines = list()
    for line in lines:
        if "//" in line and delim in line:
            line = line.split(delim)[0] + "\n"
        newlines.append(line)
    return newlines


def adjust_blank_lines(lines):
    """
        this function merges two (and more) successive blank lines
        into one blank
    """
    nitems = len(lines)
    flags = [True for _ in range(nitems)]
    for n in range(1, nitems):
        # check two neighbouring lines
        l0 = lines[n - 1]
        l1 = lines[n]
        if "\n" == l0 and "\n" == l1:
            flags[n] = False
    newlines = list()
    for line, flag in zip(lines, flags):
        if flag:
            newlines.append(line)
    return newlines


def main():
    argv = sys.argv
    # sanitise input
    assert 2 == len(argv)
    ndims = int(argv[1])
    # input source files
    fnames = list()
    fnames += get_filenames("src")
    fnames += get_filenames("include")
    for fname in fnames:
        with open(fname, "r") as f:
            lines = f.readlines()
        lines = extract_given_dim(ndims, lines)
        lines = modify_comments(lines)
        lines = adjust_blank_lines(lines)
        if 0 == len(lines):
            # nothing remains, delete file
            os.system(f"rm {fname}")
            continue
        # dump
        lines = "".join(lines)
        with open(fname, "w") as f:
            f.write(lines)


if __name__ == "__main__":
    main()
