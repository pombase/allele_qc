"""
Pydantic models to define syntax rules. This is useful for type-hinting. An example grammar can be found in grammar.py.

We define a grammar as a list of SyntaxRule objects, see the readme.
"""
from typing import Callable
from pydantic import BaseModel
import re


class SyntaxRule(BaseModel):
    type: str
    rule_name: str
    regex: str
    apply_syntax: Callable[[list[str]], str] = lambda g: ''
    check_sequence: Callable[[list[str], dict], str] = lambda g, gg: ''
    further_check: Callable[[list[str], dict], bool] = lambda g, gg: True
    format_for_transvar: Callable[[list[str], dict], list[str]] = lambda g, gg: []

    def get_groups(self, allele_sub_string: str, gene: dict) -> list[str]:
        """
        Match an allele description with the regex of this syntax rule (should match entire string), and further_checks.
        Returns the match.groups().
        """
        match = re.match(self.regex, allele_sub_string)
        if match is None:
            raise ValueError(f'allele_substring {allele_sub_string} does not match regex of {self.type}:{self.rule_name}')

        groups = match.groups()
        if self.further_check(groups, gene):
            return groups
        raise ValueError(f'allele_substring {allele_sub_string} does not match further_check rule {self.type}:{self.rule_name}')


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
