from typing import Callable
from pydantic import BaseModel


class SyntaxRule(BaseModel):
    type: str
    rule_name: str
    regex: str
    apply_syntax: Callable[[list[str]], str] = lambda g: ''
    check_invalid: Callable[[list[str]], str] = lambda g: ''
    check_sequence: Callable[[list[str, dict]], str] = lambda g, gg: ''
