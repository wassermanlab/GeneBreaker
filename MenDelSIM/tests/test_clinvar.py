import unittest
from MenDelSIM.src.clinvar import ClinVar
from MenDelSIM.src.transcript import Transcript
from MenDelSIM.src.api_helper import *


class clinvarTests(unittest.TestCase):
    ISG15_uid = get_all_transcripts("ISG15", "hg38")[0]["qualifiers"]["uid"]
    transcript = Transcript(ISG15_uid, "hg38")

    # test 1
    def test_copy_number_error(self):
        clinvar = {
            "TYPE": "ClinVar",
            "REGION": "GENIC",
            "IMPACT": {
                "clinvar_id": 475283,
                "START": 1014042,
                "REF": "G",
                "ALT": "A",
            },
            "ZYGOSITY": "HETEROZYGOUS"}
        clinvar = ClinVar(clinvar, self.transcript)
        row = clinvar.get_vcf_row().split("\t")
        self.assertEqual(row[0], "chr1")
        self.assertEqual(row[1], "1014042")
        self.assertEqual(
            row[3], "G")
        self.assertEqual(
            row[4], "A")
