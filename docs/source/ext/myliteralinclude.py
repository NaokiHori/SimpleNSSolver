import os
import re
from typing import List
from docutils import nodes
from docutils.nodes import Node
from docutils.parsers.rst import directives
from docutils.statemachine import StringList
from sphinx.directives import optional_int
from sphinx.util.docutils import SphinxDirective
from sphinx.util.typing import OptionSpec
from sphinx.directives.code import LiteralIncludeReader


def get_lines(filename, tag):
    # load all file contents
    with open(filename, "r") as f:
        lines = f.readlines()
    # check number of tags in the document
    n_tags = 0
    for line in lines:
        if tag in line:
            n_tags += 1
    # there should be at least 1
    if n_tags == 0:
        return
    # check pairs of "start lineno" and "end lineno"
    ss = list()
    es = list()
    for cnt, line in enumerate(lines):
        if f" {tag} " in line:
            if not "|" in line:
                msg = "delimiter | is not found"
                raise RuntimeError(msg)
            if not "//" in line and not ("/*" in line and "*/" in line):
                msg = "// nor /* */ are found"
                raise RuntimeError(msg)
            num = line.split("|")[1]
            num = re.sub(r"[^0-9]", "", num)
            num = int(num)
            s = cnt + 1
            e = s + num - 1
            s_string = lines[s].strip()
            e_string = lines[e].strip()
            if s_string == "{":
                if e_string == "}":
                    s += 1
                    e -= 1
                else:
                    # starting from "{", but does not end with "}"
                    # maybe making a mistake in "num" setting
                    assert(0 == 1)
            # python index starts from 0, while line number starts from 1
            s += 1
            e += 1
            ss.append(s)
            es.append(e)
    # check number of "{" and "}" are same
    for s, e in zip(ss, es):
        part_of_lines = lines[s-1:e]
        num_curly_s = 0
        num_curly_e = 0
        for l in part_of_lines:
            num_curly_s += l.count("{")
            num_curly_e += l.count("}")
        assert(num_curly_s == num_curly_e)
    # pack result and return
    retval = list()
    for s, e in zip(ss, es):
        retval.append((f"{s}-{e}", s))
    return retval

def remove_head_spaces(lines):
    lines = lines.split("\n")[:-1]
    nspaces = 0
    for lcnt, line in enumerate(lines):
        for ccnt, char in enumerate(line):
            if char != " ":
                if lcnt == 0:
                    nspaces = ccnt
                else:
                    nspaces = min(nspaces, ccnt)
                break
    newlines = list()
    for line in lines:
        newlines.append(line[nspaces:])
    newlines = "\n".join(newlines)
    return newlines

def container_wrapper(directive: SphinxDirective, literal_node: Node, caption: str) -> nodes.container:  # NOQA
    container_node = nodes.container('', literal_block=True,
                                     classes=['literal-block-wrapper'])
    parsed = nodes.Element()
    directive.state.nested_parse(StringList([caption], source=''),
                                 directive.content_offset, parsed)
    if isinstance(parsed[0], nodes.system_message):
        msg = __('Invalid caption: %s' % parsed[0].astext())
        raise ValueError(msg)
    elif isinstance(parsed[0], nodes.Element):
        caption_node = nodes.caption(parsed[0].rawsource, '',
                                     *parsed[0].children)
        caption_node.source = literal_node.source
        caption_node.line = literal_node.line
        container_node += caption_node
        container_node += literal_node
        return container_node
    else:
        raise RuntimeError  # never reached

class MyLiteralInclude(SphinxDirective):
    """
    Customised LiteralInclude class
    """
    has_content = False
    required_arguments = 1
    optional_arguments = 0
    final_argument_whitespace = True
    option_spec: OptionSpec = {
        'language': directives.unchanged_required,
        'tag': directives.unchanged_required,
        'caption': directives.unchanged,
    }
    def run(self) -> List[Node]:
        document = self.state.document
        if not document.settings.file_insertion_enabled:
            return [document.reporter.warning('File insertion disabled',
                                              line=self.lineno)]
        try:
            location = self.state_machine.get_source_and_line(self.lineno)
            rel_filename, filename = self.env.relfn2path(self.arguments[0])
            self.env.note_dependency(rel_filename)
            pairs = get_lines(filename, self.options['tag'])
            retnodes = list()
            for pair in pairs:
                self.options['lines'], self.options['lineno-start'] = pair[0], pair[1]
                reader = LiteralIncludeReader(filename, self.options, self.config)
                text, lines = reader.read(location=location)
                text = remove_head_spaces(text)
                retnode: Element = nodes.literal_block(text, text, source=filename)
                retnode['force'] = 'force' in self.options
                self.set_source_info(retnode)
                retnode['language'] = self.options['language']
                retnode['linenos'] = True
                extra_args = retnode['highlight_args'] = {}
                extra_args['linenostart'] = reader.lineno_start
                ## caption, always show filename
                # /../../src/<filename>.c -> src/<filename>.c
                if "src/" in self.arguments[0]:
                    delimiter = "src/"
                elif "include/" in self.arguments[0]:
                    delimiter = "include/"
                caption = delimiter + self.arguments[0].split(delimiter)[1]
                retnode = container_wrapper(self, retnode, caption)
                ##
                self.add_name(retnode)
                retnodes.append(retnode)
            return retnodes
        except Exception as exc:
            return [document.reporter.warning(exc, line=self.lineno)]

def setup(app):
    app.add_directive("myliteralinclude", MyLiteralInclude)

    return {
        'version': '0.1',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
