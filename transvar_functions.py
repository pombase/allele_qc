from transvar_main_script import parser_add_annotation, parser_add_mutation, parser_add_general
from transvar.anno import read_config, main_one, AnnoDB, print_header
import argparse
from functools import partial
import io
from contextlib import redirect_stdout, redirect_stderr
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


def get_anno_db() -> AnnoDB:
    annotation_parser = argparse.ArgumentParser(description=__doc__)
    parser_add_annotation(annotation_parser)
    annotation_args = annotation_parser.parse_args(['--ensembl', 'data/pombe_genome.gtf.transvardb', '--reference', 'data/pombe_genome.fa'])
    config = read_config()
    return AnnoDB(annotation_args, config)


def get_transvar_str_annotation(variant_type: str, variant_description: str, db: AnnoDB = None) -> str:

    if db is None:
        db = get_anno_db()

    if variant_type not in ['ganno', 'canno', 'panno']:
        raise ValueError("variant_type must be one of 'ganno', 'canno', 'panno'")

    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers()
    p = subparsers.add_parser('ganno', help='annotate gDNA element')
    parser_add_annotation(p)
    parser_add_mutation(p)
    parser_add_general(p)
    p.set_defaults(func=partial(main_one, db=db, at='g'))

    p = subparsers.add_parser("canno", help='annotate cDNA elements')
    parser_add_annotation(p)
    parser_add_mutation(p)
    parser_add_general(p)
    p.set_defaults(func=partial(main_one, db=db, at='c'))

    p = subparsers.add_parser("panno", help='annotate protein element')
    parser_add_annotation(p)
    parser_add_mutation(p)
    parser_add_general(p)
    p.set_defaults(func=partial(main_one, db=db, at='p'))

    # We set the -v argument to 2 (verbose), to raise errors
    args = parser.parse_args([variant_type, '-i', variant_description, '--ensembl', 'data/pombe_genome.gtf.transvardb', '--reference', 'data/pombe_genome.fa', '-v', '2'])

    # We include this to pause on errors
    args.suspend = True
    args.i = TransvarCustomString(args.i)

    output_stream = io.StringIO()
    error_stream = io.StringIO()
    with redirect_stderr(error_stream):
        with redirect_stdout(output_stream):
            if (not args.vcf) and (not args.noheader):
                print(print_header(args))
            args.func(args)

    output_str = output_stream.getvalue()
    output_stream.close()
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

    # Some cases are printed to stderr but not raised... (err_warn with I:g.1878352_1878357delCCAAAAinsTTTGGG)
    error_str = error_stream.getvalue()
    error_stream.close()
    if error_str != '':
        raise ValueError(error_str)

    return output_str

