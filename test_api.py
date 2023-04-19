from api import app
from fastapi.testclient import TestClient
import unittest

client = TestClient(app)


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
