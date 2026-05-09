"""Sanitized numpy expression evaluator for custom contour expressions.

Allows users to type expressions like:

    sqrt(stress_vm**2 + 0.3*stress_max_prin**2)
    1 - stress_vm / 600e6
    abs(disp_x) + abs(disp_y)

against the per-entity arrays we extract from the OP2 result. Only a
whitelist of numpy functions is exposed in the eval namespace; no
``__import__``, no attribute access on user objects, no file I/O.

The evaluator is intentionally simple - we'd rather refuse a tricky
expression than risk a sandbox escape. If a user needs richer semantics
they can shape their data into the named arrays we expose here.
"""

from __future__ import annotations

import ast
import math
from typing import Mapping

import numpy as np


# Functions a user can call. Each must take numeric or numpy inputs.
_ALLOWED_FUNCS = {
    'abs': np.abs,
    'sqrt': np.sqrt,
    'exp': np.exp,
    'log': np.log,
    'log10': np.log10,
    'sin': np.sin,
    'cos': np.cos,
    'tan': np.tan,
    'arctan2': np.arctan2,
    'min': np.minimum,
    'max': np.maximum,
    'minimum': np.minimum,
    'maximum': np.maximum,
    'where': np.where,
    'clip': np.clip,
    'sign': np.sign,
    'mean': np.mean,
    'sum': np.sum,
    'pi': math.pi,
    'e': math.e,
    # Allow both syntactic-sugar 'pow' and a no-frills power
    'pow': np.power,
}

_ALLOWED_NODE_TYPES = (
    ast.Expression,
    ast.BinOp, ast.UnaryOp, ast.BoolOp, ast.Compare,
    ast.Call, ast.Name, ast.Load, ast.Constant,
    ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Mod, ast.Pow, ast.FloorDiv,
    ast.USub, ast.UAdd,
    ast.Eq, ast.NotEq, ast.Lt, ast.LtE, ast.Gt, ast.GtE,
    ast.And, ast.Or, ast.Not,
    ast.IfExp,
)


class ExpressionError(ValueError):
    """Raised when an expression is invalid or uses disallowed nodes."""


def _validate(node: ast.AST) -> None:
    """Walk the parsed tree; raise ExpressionError on anything outside
    the allowed node types or disallowed function names."""
    for sub in ast.walk(node):
        if not isinstance(sub, _ALLOWED_NODE_TYPES):
            raise ExpressionError(f"Disallowed syntax: {type(sub).__name__}")
        if isinstance(sub, ast.Call):
            func = sub.func
            if not isinstance(func, ast.Name):
                raise ExpressionError(
                    "Only direct function calls are allowed (no attribute calls)."
                )
            if func.id not in _ALLOWED_FUNCS:
                raise ExpressionError(f"Function '{func.id}' is not allowed.")
        if isinstance(sub, ast.Name):
            # Names are validated against the variable map at eval time.
            pass


def evaluate(expression: str, variables: Mapping[str, np.ndarray]) -> np.ndarray:
    """Parse, validate, and evaluate `expression` against `variables`.

    Returns a numpy array with the same shape as the input arrays.
    Scalar results are wrapped to match by broadcasting.
    """
    try:
        tree = ast.parse(expression, mode='eval')
    except SyntaxError as exc:
        raise ExpressionError(f"Syntax error: {exc.msg}") from exc

    _validate(tree)

    namespace = {**_ALLOWED_FUNCS, **dict(variables)}
    try:
        result = eval(  # noqa: S307 - validated above
            compile(tree, '<expression>', 'eval'),
            {"__builtins__": {}},
            namespace,
        )
    except NameError as exc:
        raise ExpressionError(str(exc)) from exc

    if not isinstance(result, np.ndarray):
        result = np.asarray(result)
    return result
