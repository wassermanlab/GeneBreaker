import unittest
from GeneBreaker.src.clinvar import ClinVar
from GeneBreaker.src.transcript import Transcript
from GeneBreaker.src.api_helper import *


class clinvarTests(unittest.TestCase):
    SOX9_uid = get_all_transcripts("SOX9", "hg38")[0]["qualifiers"]["uid"]
    transcript = Transcript(SOX9_uid, "hg38")

    # test 1
    def test_clinvar(self):
        clinvar = {
            "TYPE": "CLINVAR",
            "REGION": "GENIC",
            "IMPACT": {
                "CLINVAR_ID": 450265,
            },
            "ZYGOSITY": "HETEROZYGOUS"}
        clinvar = ClinVar(clinvar, self.transcript)
        row = clinvar.get_vcf_row()
        self.assertEqual(row["chrom"], "chr17")
        self.assertEqual(row["pos"], "72121409")
        self.assertEqual(row["ref"], "C")
        self.assertEqual(row["alt"], "T")
