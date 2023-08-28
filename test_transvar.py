from transvar_functions import get_transvar_str_annotation, parse_transvar_string, get_anno_db
import unittest
from transvar.err import SequenceRetrievalError, InvalidInputError

db = get_anno_db('data/pombe_genome.gtf.transvardb', 'data/pombe_genome.fa')


class TransvarTest(unittest.TestCase):

    def test_ganno(self):

        # A correct annotation that should return three transcripts (SPBC1198.04c.1, SPBC1198.04c.2 and SPNCRNA.1321.1)
        variant_list = parse_transvar_string(get_transvar_str_annotation('ganno', 'II:g.178497T>A', db))
        self.assertEqual(len(variant_list), 3)

        # A wrong chromosome / sequence returns a SequenceRetrievalError
        try:
            get_transvar_str_annotation('ganno', 'blah:g.1000A>T', db)
        except SequenceRetrievalError as e:
            self.assertEqual(e.args[0], 'chromosome blah not found in reference')
        else:
            self.fail("SequenceRetrievalError not raised")

        # A wrong position returns a generic exception
        try:
            get_transvar_str_annotation('ganno', 'I:g.1000A>T', db)
        except Exception as e:
            self.assertEqual(e.args[0], "invalid_reference_base_A;expect_G")
        else:
            self.fail("Exception not raised")

        # A wrongly formatted ganno returns an InvalidInputError
        try:
            get_transvar_str_annotation('ganno', 'dummy', db)
        except InvalidInputError as e:
            self.assertEqual(e.args[0], "invalid_mutation_string: dummy (type:g)")
        else:
            self.fail("InvalidInputError not raised")

        try:
            get_transvar_str_annotation('ganno', 'I:blah', db)
        except InvalidInputError as e:
            self.assertEqual(e.args[0], "invalid_mutation_string_blah")
        else:
            self.fail("InvalidInputError not raised")

        # An empty insertion returns an InvalidInputError
        try:
            get_transvar_str_annotation('ganno', 'I:g.1000_1001ins', db)
        except InvalidInputError as e:
            self.assertEqual(e.args[0], "insertion_without_inserted_sequence_g.1000_1001ins")
        else:
            self.fail("InvalidInputError not raised")

    def test_panno(self):

        # A wrong chromosome / sequence returns a SequenceRetrievalError
        try:
            get_transvar_str_annotation('panno', 'blah:p.T566S', db)
        except Exception as e:
            self.assertEqual(e.args[0], 'invalid_gene_blah')
        else:
            self.fail("Exception not raised")

        # A wrong position returns a ValueError (this is not the case by default, I added it,
        # see transvar_functions.py)
        try:
            get_transvar_str_annotation('panno', 'SPAPB1A10.09:p.A2F', db)
        except ValueError as e:
            self.assertEqual(e.args[0], "no_valid_transcript_found")
        else:
            self.fail("ValueError not raised")

        # A wrongly formatted panno returns a InvalidInputError
        try:
            get_transvar_str_annotation('panno', 'SPAPB1A10.09:blah', db)
        except InvalidInputError as e:
            self.assertEqual(e.args[0], "invalid_mutation_string_blah")
        else:
            self.fail("InvalidInputError not raised")

        # An empty insertion does not return an error in panno
        get_transvar_str_annotation('panno', 'SPAPB1A10.09:p.10_11ins', db)

        # A correct annotation that should return two transcripts
        variant_list = parse_transvar_string(get_transvar_str_annotation('panno', 'SPBC1198.04c:p.N3A', db))
        self.assertEqual(len(variant_list), 2)


