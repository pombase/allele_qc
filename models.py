"""
Pydantic models to define syntax rules. This is useful for type-hinting. An example grammar can be found in grammar.py.

We define a grammar as a list of SyntaxRule objects, see the readme.
"""
from typing import Callable
from pydantic import BaseModel


class SyntaxRule(BaseModel):
    type: str
    rule_name: str
    regex: str
    apply_syntax: Callable[[list[str]], str] = lambda g: ''
    check_invalid: Callable[[list[str]], str] = lambda g: ''
    check_sequence: Callable[[list[str], dict], str] = lambda g, gg: ''
    coordinate_indexes: tuple[int, ...] = ()


def find_rule(grammar: list[SyntaxRule], rule_type, rule_name) -> SyntaxRule:
    for rule in grammar:
        if rule.type == rule_type and rule.rule_name == rule_name:
            return rule
    return None
