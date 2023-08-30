from api import app, CheckAlleleDescriptionResponse, CheckModificationResponse, AlleleFix, OldCoordsFix
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
        response = client.get("/check_allele", params={'systematic_id': 'dummy', 'allele_description': 'V123A,PLR140AAA,150-600', 'allele_type': 'partial_amino_acid_deletion', 'allele_name': ''})
        self.assertEqual(response.status_code, 404)

    def test_valid_id(self):
        # Should contain a valid response
        response = client.get("/check_allele", params={'systematic_id': 'SPBC359.03c', 'allele_description': 'V123A,PLR140AAA,150-600', 'allele_type': 'partial_amino_acid_deletion', 'allele_name': ''})
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
        # The mod code is not really used in the test
        response = client.get("/check_modification", params={'systematic_id': 'dummy', 'sequence_position': 'V123; V124,V125', 'mod_code': 'MOD:00046'})
        self.assertEqual(response.status_code, 404)

    def test_valid_id(self):
        # Should contain a valid response
        response = client.get("/check_modification", params={'systematic_id': 'SPBC359.03c', 'sequence_position': 'V123; V124,V125', 'mod_code': 'MOD:00046'})
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

    # Common part of the tests ================================================

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

    def test_valid_id(self):
        # Should contain a valid response and support genes with multiple transcripts
        for systematic_id in ['SPBC359.03c', 'SPAC22A12.08c', 'SPAC22A12.08c.2']:
            for ep in self.entrypoints:
                response = client.get(ep, params={'systematic_id': systematic_id, 'targets': 'S123,A124,N125'})
                self.assertEqual(response.status_code, 200)
                try:
                    resp = [AlleleFix.parse_obj(ele) for ele in response.json()]
                except Exception as e:
                    self.fail(e)

    # Endpoint-specific part of the tests ================================================

    def test_multi_shift(self):
        response = client.get("/multi_shift_fix", params={'systematic_id': 'SPAPB1A10.09', 'targets': 'S123,A124,N125'})

        try:
            resp = [AlleleFix.parse_obj(ele) for ele in response.json()]
        except Exception as e:
            self.fail(e)
        # There should be two possible fixes for this example
        self.assertEqual(len(resp), 2)

        # This syntax should also work
        response = client.get("/multi_shift_fix", params={'systematic_id': 'SPAPB1A10.09', 'targets': 'SAN-123-AAA'})

        try:
            resp = [AlleleFix.parse_obj(ele) for ele in response.json()]
        except Exception as e:
            self.fail(e)
        self.assertEqual(len(resp), 2)

        # Should return empty list if 2 or less positions are provided
        response = client.get("/multi_shift_fix", params={'systematic_id': 'SPAPB1A10.09', 'targets': 'S123,A124'})
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json(), [])

        # Also with the alternative syntax
        response = client.get("/multi_shift_fix", params={'systematic_id': 'SPAPB1A10.09', 'targets': 'SA-123-AV'})
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json(), [])

        # With the DNA case
        response = client.get("/multi_shift_fix", params={'systematic_id': 'SPAPB1A10.09', 'targets': 'GTATGTAATCGT-424-AAAAAAAAAAAA', 'dna_or_protein': 'dna'})
        self.assertTrue(response.json()[0]['values'].startswith('G425A'))
        self.assertEqual(response.status_code, 200)

    def test_old_coords(self):
        response = client.get("/old_coords_fix", params={'systematic_id': 'SPBC1706.01', 'targets': 'P170A,V223A,F225A,AEY-171-LLL'})

        try:
            resp = [OldCoordsFix.parse_obj(ele) for ele in response.json()]
        except Exception as e:
            self.fail(e)
        self.assertEqual(len(resp), 1)
        self.assertEqual(resp[0].values, 'P182A,V235A,F237A,A183L,E184L,Y185L')
        self.assertEqual(resp[0].revision, '20110324')
        self.assertEqual(resp[0].location, '588765..591194')

    def test_histone(self):
        response = client.get("/histone_fix", params={'systematic_id': 'SPAC1834.04', 'targets': 'ART-1-LLL,K9A,K14R,K14A'})

        try:
            resp = [AlleleFix.parse_obj(ele) for ele in response.json()]
        except Exception as e:
            self.fail(e)
        self.assertEqual(len(resp), 1)
        self.assertEqual(resp[0].values, 'A2L,R3L,T4L,K10A,K15R,K15A')


class ResidueAtPositionTest(unittest.TestCase):

    def test_valid_id(self):
        # Request protein
        response = client.get("/residue_at_position", params={'systematic_id': 'SPAPB1A10.09', 'position': 1, 'dna_or_protein': 'protein'})
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.text, 'M')

        # Request DNA on +1 strand, positive residues
        for i, res in enumerate('ATGCAAACAG'):
            response = client.get("/residue_at_position", params={'systematic_id': 'SPAPB1A10.09', 'position': i + 1, 'dna_or_protein': 'dna'})
            self.assertEqual(response.status_code, 200)
            self.assertEqual(response.text, res)

        # Negative residues
        for i, res in enumerate('TCAACCGGTTCA'[::-1]):
            response = client.get("/residue_at_position", params={'systematic_id': 'SPAPB1A10.09', 'position': -(i + 1), 'dna_or_protein': 'dna'})
            self.assertEqual(response.status_code, 200)
            self.assertEqual(response.text, res)

        # Request DNA on -1 strand, positive residues
        for i, res in enumerate('ATGTCGGCTCAG'):
            response = client.get("/residue_at_position", params={'systematic_id': 'SPAPB1A10.10c', 'position': i + 1, 'dna_or_protein': 'dna'})
            self.assertEqual(response.status_code, 200)
            self.assertEqual(response.text, res)

        # Negative residues GACAGAATATACTTCA
        for i, res in enumerate('GACAGAATATACTTCA'[::-1]):
            response = client.get("/residue_at_position", params={'systematic_id': 'SPAPB1A10.10c', 'position': -(i + 1), 'dna_or_protein': 'dna'})
            self.assertEqual(response.status_code, 200)
            self.assertEqual(response.text, res)

class TransvarEntryPointsTest(unittest.TestCase):

    def test_allele_entrypoint(self):

        response = client.get("/allele_transvar_coordinates", params={'systematic_id': 'SPBC359.03c', 'allele_description': 'A3V,SEA23PPP,150-200', 'allele_type': 'amino_acid_deletion_and_mutation', 'allele_name': 'blah'})
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.text, 'II:g.113041G>A/c.8C>T/p.A3V|II:g.112974_112982delinsAGGAGGAGG/c.67_75delinsCCTCCTCCT/p.S23_A25delinsPPP|II:g.112449_112601del153/c.452_604del153/p.D151_I201del51')

        response = client.get("/allele_transvar_coordinates", params={'systematic_id': 'SPAC688.08', 'allele_description': 'S1137A', 'allele_type': 'amino_acid_mutation', 'allele_name': 'blah'})
        self.assertEqual(response.status_code, 400)

    def test_modification_entrypoint(self):

        response = client.get('/protein_modification_transvar_coordinates', params={'systematic_id': 'SPBC359.03c', 'sequence_position': 'A3,K4'})
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.text, 'II:g.113040_113042/c.7_9/p.3A|II:g.113037_113039/c.10_12/p.4K')

        response = client.get('/protein_modification_transvar_coordinates', params={'systematic_id': 'SPAC688.08', 'sequence_position': 'S1137'})
        self.assertEqual(response.status_code, 400)

    def test_anno_entry_points(self):
        response = client.get('/ganno', params={'variant_description': 'II:g.178497T>A'})
        self.assertEqual(response.status_code, 200)
        self.assertEqual(len(response.json()), 3)

        response = client.get('/panno', params={'variant_description': 'SPBC1198.04c:p.N3A'})
        self.assertEqual(response.status_code, 200)
        self.assertEqual(len(response.json()), 2)

        response = client.get('/canno', params={'variant_description': 'SPAC3F10.09:c.5A>T'})
        self.assertEqual(response.status_code, 200)
        self.assertEqual(len(response.json()), 1)
