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

    result = [TransvarAnnotation.from_list(t.split('\t')) for t in transvar_list]
    # Scan for errors
    # for t in result:
    #     if 'Error=' in t.info:
    #         raise ValueError(t.info)

    return result


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

    # We set the -v argument to 2 (verbose), to raise errors
    args = parser.parse_args([variant_type, '-i', variant_description, '--ensembl', 'data/pombe_genome.gtf.transvardb', '--reference', 'data/pombe_genome.fa', '-v', '2'])

    # We include this to pause on errors
    args.suspend = True
    args.i = TransvarCustomString(args.i)

    output_stream = io.StringIO()
    with redirect_stdout(output_stream):
        args.func(args)

    output_str = output_stream.getvalue()
    # Extra error handling
    if variant_type == 'panno':
        # In the case where the indicated positions don't match any transcript of the gene, transvar returns coordinates(gDNA/cDNA/protein) = `././.`
        # and info = 'no_valid_transcript_found'. Maybe there is a special case where the info is different?

        transvar_fields_first_row = output_str.split('\n')[1].split('\t')
        if transvar_fields_first_row[-3] == '././.':
            if transvar_fields_first_row[-1] == 'no_valid_transcript_found':
                raise ValueError('no_valid_transcript_found')
            else:
                raise ValueError('Unknown error: ', transvar_fields_first_row[-1])

    return output_str

