import os
from typing import List
from docutils import nodes
from docutils.nodes import Node
from docutils.parsers.rst import directives
from sphinx.directives import optional_int
from sphinx.util.docutils import SphinxDirective
from sphinx.util.typing import OptionSpec
from sphinx.directives.code import LiteralIncludeReader

def get_lines(filename, tag):
    # find "/**" and tag (declaration of function)
    # l0: line number of the upper-closest "/**"
    # l1: line number of the end of function
    l0 = -1
    l1 = -1
    with open(filename, "r") as f:
        in_def = False
        for cnt, line in enumerate(f):
            if "/**" in line:
                l0 = cnt
            if f"{tag}(" in line:
                # function name is found
                in_def = True
            if in_def and "){" in line:
                l1 = cnt
                break
    # index in python starts from 0, while line number is from 1
    l0 += 1
    l1 += 1
    if l0 <= 0:
        msg = f"non-positive l0: {l0}"
        raise RuntimeError(msg)
    if l1 <= 0:
        msg = f"non-positive l1: {l1}"
        raise RuntimeError(msg)
    if l1 <= l0:
        msg = f"l1 ({l1}) < l0 {l0}"
        raise RuntimeError(msg)
    return "{}-{}".format(l0, l1)

def reshape_text(text):
    # replace last "){" with ");"
    return text.replace("){", ");")

class MyDeclare(SphinxDirective):
    """
    Customised LiteralInclude class
    """
    has_content = False
    required_arguments = 1
    optional_arguments = 0
    final_argument_whitespace = True
    option_spec: OptionSpec = {
        "language": directives.unchanged_required,
        "tag": directives.unchanged_required,
    }
    def run(self) -> List[Node]:
        document = self.state.document
        if not document.settings.file_insertion_enabled:
            return [document.reporter.warning("File insertion disabled",
                                              line=self.lineno)]
        try:
            rel_filename, filename = self.env.relfn2path(self.arguments[0])
            self.env.note_dependency(rel_filename)
            self.options["lines"] = get_lines(filename, self.options["tag"])
            reader = LiteralIncludeReader(filename, self.options, self.config)
            text, lines = reader.read()
            text = reshape_text(text)
            retnode: Element = nodes.literal_block(text, text, source=filename)
            retnode["force"] = "force" in self.options
            self.set_source_info(retnode)
            retnode["language"] = self.options["language"]
            retnode["linenos"] = False
            self.add_name(retnode)
            return [retnode]
        except Exception as exc:
            return [document.reporter.warning(exc, line=self.lineno)]

def setup(app):
    app.add_directive("mydeclare", MyDeclare)
    return {
        "version": "0.1",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }

