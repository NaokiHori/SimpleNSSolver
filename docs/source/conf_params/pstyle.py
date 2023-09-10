# ref: alabaster/support.py

from pygments.style import Style
from pygments.token import (
    Keyword,
    Name,
    Comment,
    String,
    Error,
    Number,
    Operator,
    Generic,
    Whitespace,
    Punctuation,
    Other,
    Literal,
)

# color fliped

class MyAlabaster(Style):
    background_color = "#070707"  # doesn't seem to override CSS 'pre' styling?
    default_style = ""
    styles = {
        # No corresponding class for the following:
        # Text:                     "", # class:  ''
        # Whitespace: "underline #070707",  # class: 'w'
        Error: "#5BFFFF border:#10D6D6",  # class: 'err'
        Other: "#FFFFFF",  # class 'x'
        Comment: "italic #70A6FD",  # class: 'c'
        Comment.Preproc: "noitalic",  # class: 'cp'
        Keyword: "bold #FFBB9E",  # class: 'k'
        Keyword.Constant: "bold #FFBB9E",  # class: 'kc'
        Keyword.Declaration: "bold #FFBB9E",  # class: 'kd'
        Keyword.Namespace: "bold #FFBB9E",  # class: 'kn'
        Keyword.Pseudo: "bold #FFBB9E",  # class: 'kp'
        Keyword.Reserved: "bold #FFBB9E",  # class: 'kr'
        Keyword.Type: "bold #FFBB9E",  # class: 'kt'
        Operator: "#A7D7FF",  # class: 'o'
        Operator.Word: "bold #FFBB9E",  # class: 'ow' - like keywords
        Punctuation: "bold #FFFFFF",  # class: 'p'
        # because special names such as Name.Class, Name.Function, etc.
        # are not recognized as such later in the parsing, we choose them
        # to look the same as ordinary variables.
        Name: "#FFFFFF",  # class: 'n'
        Name.Attribute: "#3B5FFF",  # class: 'na' - to be revised
        Name.Builtin: "#FFBB9E",  # class: 'nb'
        Name.Builtin.Pseudo: "#CB9A5B",  # class: 'bp'
        Name.Class: "#FFFFFF",  # class: 'nc' - to be revised
        Name.Constant: "#FFFFFF",  # class: 'no' - to be revised
        Name.Decorator: "#888",  # class: 'nd' - to be revised
        Name.Entity: "#31A3FF",  # class: 'ni'
        Name.Exception: "bold #33FFFF",  # class: 'ne'
        Name.Function: "#FFFFFF",  # class: 'nf'
        Name.Property: "#FFFFFF",  # class: 'py'
        Name.Label: "#0A86FF",  # class: 'nl'
        Name.Namespace: "#FFFFFF",  # class: 'nn' - to be revised
        Name.Other: "#FFFFFF",  # class: 'nx'
        Name.Tag: "bold #FFBB9E",  # class: 'nt' - like a keyword
        Name.Variable: "#FFFFFF",  # class: 'nv' - to be revised
        Name.Variable.Class: "#FFFFFF",  # class: 'vc' - to be revised
        Name.Variable.Global: "#FFFFFF",  # class: 'vg' - to be revised
        Name.Variable.Instance: "#FFFFFF",  # class: 'vi' - to be revised
        Number: "#66FFFF",  # class: 'm'
        Literal: "#FFFFFF",  # class: 'l'
        Literal.Date: "#FFFFFF",  # class: 'ld'
        String: "#B165F9",  # class: 's'
        String.Backtick: "#B165F9",  # class: 'sb'
        String.Char: "#B165F9",  # class: 'sc'
        String.Doc: "italic #70A6FD",  # class: 'sd' - like a comment
        String.Double: "#B165F9",  # class: 's2'
        String.Escape: "#B165F9",  # class: 'se'
        String.Heredoc: "#B165F9",  # class: 'sh'
        String.Interpol: "#B165F9",  # class: 'si'
        String.Other: "#B165F9",  # class: 'sx'
        String.Regex: "#B165F9",  # class: 'sr'
        String.Single: "#B165F9",  # class: 's1'
        String.Symbol: "#B165F9",  # class: 'ss'
        Generic: "#FFFFFF",  # class: 'g'
        Generic.Deleted: "#5BFFFF",  # class: 'gd'
        Generic.Emph: "italic #FFFFFF",  # class: 'ge'
        Generic.Error: "#10D6D6",  # class: 'gr'
        Generic.Heading: "bold #FFFF7F",  # class: 'gh'
        Generic.Inserted: "#FF5FFF",  # class: 'gi'
        Generic.Output: "#888",  # class: 'go'
        Generic.Prompt: "#8BACCB",  # class: 'gp'
        Generic.Strong: "bold #FFFFFF",  # class: 'gs'
        Generic.Subheading: "bold #7FFF7F",  # class: 'gu'
        Generic.Traceback: "bold #5BFFFF",  # class: 'gt'
    }

