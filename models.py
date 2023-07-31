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
    check_sequence: Callable[[list[str], dict], str] = lambda g, gg: ''
    further_check: Callable[[list[str]], bool] = lambda g: True


class AllowedTypes(BaseModel):
    allowed_types: dict[frozenset, str]
    composed_types: dict[str, list[str]]

    def split_composed_types(self, types: list[str]) -> list[str]:
        """
        Split a list of types into a list of types that are not composed types, and a list of composed types.
        """
        splitted_types = []
        for type in types:
            if type in self.composed_types:
                splitted_types += self.composed_types[type]
            else:
                splitted_types.append(type)
        return frozenset(splitted_types)

    def __getitem__(self, types: frozenset) -> str:
        splitted_types = self.split_composed_types(types)
        return self.allowed_types[splitted_types]



def find_rule(grammar: list[SyntaxRule], rule_type, rule_name) -> SyntaxRule:
    for rule in grammar:
        if rule.type == rule_type and rule.rule_name == rule_name:
            return rule
    return None
