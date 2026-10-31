"""Cyclomatic complexity per function, Python standard library only.

    python3 tools/complexity.py <file_or_dir> [...] [--min 10]

Complexity = 1 + decision points: if/elif, for/while (and their else),
except handlers, boolean operators (each extra operand), conditional
expressions, comprehension ifs, match cases, assert. This is McCabe's count
as radon computes it, close enough to rank the functions to split.
"""
import ast
import os
import sys


class _Counter(ast.NodeVisitor):
    def __init__(self):
        self.n = 1

    def generic_visit(self, node):
        if isinstance(node, (ast.If, ast.For, ast.AsyncFor, ast.While, ast.IfExp,
                             ast.ExceptHandler, ast.Assert)):
            self.n += 1
            if isinstance(node, (ast.For, ast.AsyncFor, ast.While)) and node.orelse:
                self.n += 1
        elif isinstance(node, ast.BoolOp):
            self.n += len(node.values) - 1
        elif isinstance(node, ast.comprehension):
            self.n += len(node.ifs)
        elif hasattr(ast, "match_case") and isinstance(node, ast.match_case):
            self.n += 1
        # Nested function definitions are counted on their own.
        for child in ast.iter_child_nodes(node):
            if not isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef, ast.Lambda)):
                self.visit(child)


def functions(path):
    tree = ast.parse(open(path).read(), filename=path)
    out = []

    def walk(node, prefix):
        for child in ast.iter_child_nodes(node):
            if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)):
                c = _Counter()
                for stmt in child.body:
                    c.visit(stmt)
                name = prefix + child.name
                out.append((c.n, path, child.lineno, name,
                            (child.end_lineno or child.lineno) - child.lineno + 1))
                walk(child, name + ".")
            elif isinstance(child, ast.ClassDef):
                walk(child, prefix + child.name + ".")
    walk(tree, "")
    return out


def main(argv):
    mn = 10
    if "--min" in argv:
        i = argv.index("--min"); mn = int(argv[i + 1]); argv = argv[:i] + argv[i + 2:]
    files = []
    for a in argv:
        if os.path.isdir(a):
            for r, _, fs in os.walk(a):
                files += [os.path.join(r, f) for f in fs if f.endswith(".py")]
        else:
            files.append(a)
    rows = sorted((x for f in files for x in functions(f)), reverse=True)
    for n, path, line, name, length in rows:
        if n >= mn:
            print("%4d  %4d lines  %s:%d  %s" % (n, length, os.path.relpath(path), line, name))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
