"""Small mathematical expression interpreter for user-defined component laws.

No Python eval/exec, attributes, indexing, imports, comprehensions, or objects
are available. Names and ordinary scalar arithmetic are explicitly checked.
"""
from __future__ import annotations

import ast
import math
import operator
from numbers import Real


FUNCTIONS = {
    "sqrt": math.sqrt, "exp": math.exp, "log": math.log,
    "log10": math.log10, "abs": abs, "min": min, "max": max,
    "sin": math.sin, "cos": math.cos, "tanh": math.tanh,
}
CONSTANTS = {"pi": math.pi, "e": math.e}
OPERATORS = {ast.Add: operator.add, ast.Sub: operator.sub,
             ast.Mult: operator.mul, ast.Div: operator.truediv,
             ast.Pow: operator.pow}


class ExpressionError(ValueError):
    """An invalid or numerically undefined user expression."""


def _finite(value):
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ExpressionError("Expression must produce a real number.")
    try:
        value = float(value)
    except (ValueError, OverflowError) as error:
        raise ExpressionError("Expression is outside the numerical range.") from error
    if not math.isfinite(value):
        raise ExpressionError("Expression must produce a finite number.")
    return value


class ScalarExpression:
    """Compile a bounded AST once, then interpret with supplied scalar names."""

    def __init__(self, expression: str, names=()):
        if not isinstance(expression, str) or not expression.strip():
            raise ExpressionError("Enter a mathematical expression.")
        if len(expression) > 512:
            raise ExpressionError("Expression exceeds 512 characters.")
        self.expression = expression
        self.names = frozenset(names)
        try:
            self.tree = ast.parse(expression, mode="eval").body
        except (SyntaxError, RecursionError) as error:
            raise ExpressionError("Invalid mathematical expression syntax.") from error
        if sum(1 for _ in ast.walk(self.tree)) > 128:
            raise ExpressionError("Expression is too complex.")
        self._validate(self.tree)

    def _validate(self, node, depth=0):
        if depth > 24:
            raise ExpressionError("Expression nesting is too deep.")
        if isinstance(node, ast.Constant):
            _finite(node.value)
        elif isinstance(node, ast.Name):
            if node.id not in self.names and node.id not in CONSTANTS:
                raise ExpressionError(f"Unknown variable: {node.id}.")
        elif isinstance(node, ast.BinOp) and type(node.op) in OPERATORS:
            self._validate(node.left, depth + 1)
            self._validate(node.right, depth + 1)
        elif isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
            self._validate(node.operand, depth + 1)
        elif isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
            if node.func.id not in FUNCTIONS or node.keywords or not 1 <= len(node.args) <= 8:
                raise ExpressionError("Only listed mathematical functions are allowed.")
            for argument in node.args:
                self._validate(argument, depth + 1)
        else:
            raise ExpressionError("Only scalar arithmetic and mathematical functions are allowed.")

    def __call__(self, **variables):
        values = {name: _finite(value) for name, value in variables.items()}
        missing = self.names - values.keys()
        # Unused available variable names need not be passed.
        required = {n.id for n in ast.walk(self.tree) if isinstance(n, ast.Name)}
        if missing & required:
            raise ExpressionError("Missing variables: " + ", ".join(sorted(missing & required)))

        def visit(node):
            if isinstance(node, ast.Constant):
                return float(node.value)
            if isinstance(node, ast.Name):
                return values[node.id] if node.id in values else CONSTANTS[node.id]
            if isinstance(node, ast.UnaryOp):
                value = visit(node.operand)
                return -value if isinstance(node.op, ast.USub) else value
            if isinstance(node, ast.BinOp):
                left, right = visit(node.left), visit(node.right)
                if isinstance(node.op, ast.Pow) and abs(right) > 64:
                    raise ExpressionError("Exponent magnitude may not exceed 64.")
                return _finite(OPERATORS[type(node.op)](left, right))
            return _finite(FUNCTIONS[node.func.id](*(visit(arg) for arg in node.args)))

        try:
            return _finite(visit(self.tree))
        except ExpressionError:
            raise
        except (ArithmeticError, ValueError, TypeError, KeyError) as error:
            raise ExpressionError(f"Expression is undefined: {error}") from error


def evaluate_expression(expression: str, variables: dict) -> float:
    return ScalarExpression(expression, variables.keys())(**variables)
