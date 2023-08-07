from transvar_main_script import parser_add_annotation, parser_add_mutation, parser_add_general
from transvar.anno import main_anno
import argparse
from functools import partial
import io
from contextlib import redirect_stdout
from pydantic import BaseModel


class TransvarAnnotation(BaseModel):
    input: str
    transcript: str
    gene: str
    strand: str
    coordinates: str
    region: str
    info: str

    def from_list(t: list[str]) -> 'TransvarAnnotation':
        return TransvarAnnotation(input=t[0], transcript=t[1], gene=t[2], strand=t[3], coordinates=t[4], region=t[5], info=t[6])


def parse_transvar_string(transvar_str: str) -> list[TransvarAnnotation]:

    # We ommit the first and last lines of the output, which are the header and the empty line
    transvar_list = transvar_str.split('\n')[1:-1]
    if len(transvar_list) == 0:
        raise ValueError("Invalid variant description")

    return [TransvarAnnotation.from_list(t.split('\t')) for t in transvar_list]


class TransvarCustomString(str):
    """Hacky class to circunvent https://github.com/zwdzwd/transvar/issues/59
    """
    def upper(self):
        return self

    def strip(self, __chars=None):
        return TransvarCustomString(str.strip(self, __chars))

    def split(self, __sep=None, __maxsplit=-1):
        return [TransvarCustomString(x) for x in str.split(self, __sep, __maxsplit)]


def get_transvar_str_annotation(variant_type: str, variant_description: str) -> str:

    if variant_type not in ['ganno', 'canno', 'panno']:
        raise ValueError("variant_type must be one of 'ganno', 'canno', 'panno'")

    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers()
    p = subparsers.add_parser('ganno', help='annotate gDNA element')
    parser_add_annotation(p)
    parser_add_mutation(p)
    parser_add_general(p)
    p.set_defaults(func=partial(main_anno, at='g'))

    p = subparsers.add_parser("canno", help='annotate cDNA elements')
    parser_add_annotation(p)
    parser_add_mutation(p)
    parser_add_general(p)
    p.set_defaults(func=partial(main_anno, at='c'))

    p = subparsers.add_parser("panno", help='annotate protein element')
    parser_add_annotation(p)
    parser_add_mutation(p)
    parser_add_general(p)
    p.set_defaults(func=partial(main_anno, at='p'))

    args = parser.parse_args([variant_type, '-i', variant_description, '--ensembl', 'data/pombe_genome.gtf.transvardb', '--reference', 'data/pombe_genome.fa'])
    args.i = TransvarCustomString(args.i)

    output_stream = io.StringIO()
    with redirect_stdout(output_stream):
        args.func(args)
    return output_stream.getvalue()

