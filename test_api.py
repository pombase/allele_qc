from api import app, CheckAlleleDescriptionResponse, CheckModificationResponse, AlleleFix
from fastapi.testclient import TestClient
import unittest

client = TestClient(app)

# NOTE: these tests may fail if the genome data is updated and should be updated accordingly.


# test /check_modification endpoint from api.py
class GetGenomeRegionTest(unittest.TestCase):

    def test_all_formats(self):
        for format in ['fasta', 'genbank', 'embl']:
            response = client.get("/genome_region", params={'systematic_id': 'SPAPB1A10.09', 'upstream': 0, 'downstream': 0, 'format': format})
            self.assertEqual(response.status_code, 200)

    def test_multiple_transcripts(self):
        # Should return the longest transcript, in this case SPAC22A12.08c.1
        response1 = client.get("/genome_region", params={'systematic_id': 'SPAC22A12.08c', 'upstream': 0, 'downstream': 0, 'format': 'fasta'})
        self.assertEqual(response1.status_code, 200)

        response2 = client.get("/genome_region", params={'systematic_id': 'SPAC22A12.08c.1', 'upstream': 0, 'downstream': 0, 'format': 'fasta'})
        self.assertEqual(response2.status_code, 200)
        self.assertEqual(response1.text, response2.text)

    def test_dummy_id(self):
        # Should return 404
        response = client.get("/genome_region", params={'systematic_id': 'dummy', 'upstream': 0, 'downstream': 0, 'format': 'fasta'})
        self.assertEqual(response.status_code, 404)


# test /check_allele endpoint from api.py
class CheckAlleleTest(unittest.TestCase):

    def test_dummy_id(self):
        # Should return 404
        response = client.get("/check_allele", params={'systematic_id': 'dummy', 'allele_description': 'V123A,PLR-140-AAA,150-600', 'allele_type': 'partial_amino_acid_deletion'})
        self.assertEqual(response.status_code, 404)

    def test_valid_id(self):
        # Should contain a valid response
        response = client.get("/check_allele", params={'systematic_id': 'SPBC359.03c', 'allele_description': 'V123A,PLR-140-AAA,150-600', 'allele_type': 'partial_amino_acid_deletion'})
        self.assertEqual(response.status_code, 200)

        try:
            resp = CheckAlleleDescriptionResponse.parse_obj(response.json())
        except Exception as e:
            self.fail(e)

        self.assertEqual(resp.invalid_error, '')
        self.assertEqual(resp.rules_applied, 'amino_acid_mutation:single_aa|amino_acid_mutation:multiple_aa|partial_amino_acid_deletion:multiple_aa')
        self.assertEqual(resp.change_type_to, 'amino_acid_deletion_and_mutation')
        self.assertEqual(resp.change_description_to, '')
        self.assertIsInstance(resp.user_friendly_fields, CheckAlleleDescriptionResponse)


# test /check_modification endpoint from api.py
class CheckModificationTest(unittest.TestCase):

    def test_dummy_id(self):
        # Should return 404
        response = client.get("/check_modification", params={'systematic_id': 'dummy', 'sequence_position': 'V123; V124,V125'})
        self.assertEqual(response.status_code, 404)

    def test_valid_id(self):
        # Should contain a valid response
        response = client.get("/check_modification", params={'systematic_id': 'SPBC359.03c', 'sequence_position': 'V123; V124,V125'})
        self.assertEqual(response.status_code, 200)

        try:
            resp = CheckModificationResponse.parse_obj(response.json())
        except Exception as e:
            self.fail(e)

        self.assertEqual(resp.change_sequence_position_to, 'V123,V124,V125')
        self.assertEqual(resp.sequence_error, 'V123|V124|')
        self.assertIsInstance(resp.user_friendly_fields, CheckModificationResponse)


# test fixing endpoints from api.py
class FixTest(unittest.TestCase):

    entrypoints = ['/multi_shift_fix', '/old_coords_fix', '/histone_fix']

    # Common part of the tests
    def test_dummy_id(self):
        # Should return 404
        for ep in self.entrypoints:
            response = client.get(ep, params={'systematic_id': 'dummy', 'targets': 'S123,A124,N125'})
        self.assertEqual(response.status_code, 404)

    def test_targets_format_error(self):
        # Should return 422
        for ep in self.entrypoints:
            response = client.get(ep, params={'systematic_id': 'SPAPB1A10.09', 'targets': 'S123,A124,N125, 123'})
            self.assertEqual(response.status_code, 422)

    def test_few_positions(self):
        # Should return empty list if 2 or less positions are provided
        for ep in self.entrypoints:
            response = client.get(ep, params={'systematic_id': 'SPAPB1A10.09', 'targets': 'S123,A124'})
            self.assertEqual(response.status_code, 200)
            self.assertEqual(response.json(), [])

        # Also with the alternative syntax
        for ep in self.entrypoints:
            response = client.get(ep, params={'systematic_id': 'SPAPB1A10.09', 'targets': 'SA-123-AV'})
            self.assertEqual(response.status_code, 200)
            self.assertEqual(response.json(), [])

    def test_valid_id(self):
        # Should contain a valid response
        for ep in self.entrypoints:
            response = client.get(ep, params={'systematic_id': 'SPAPB1A10.09', 'targets': 'S123,A124,N125'})
            self.assertEqual(response.status_code, 200)

            try:
                resp = [AlleleFix.parse_obj(ele) for ele in response.json()]
            except Exception as e:
                self.fail(e)

    def test_multi_shift(self):
        # Should contain a valid response
        pass
    # # There should be two possible fixes for this example
    # self.assertEqual(len(resp), 2)

    # # This syntax should also work
    # response = client.get("/multi_shift_fix", params={'systematic_id': 'SPAPB1A10.09', 'targets': 'SAN-123-AAA'})

    # try:
    #     resp = [AlleleFix.parse_obj(ele) for ele in response.json()]
    # except Exception as e:
    #     self.fail(e)
    # self.assertEqual(len(resp), 2)
