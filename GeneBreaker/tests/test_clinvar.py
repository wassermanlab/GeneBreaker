import unittest
from GeneBreaker.src.clinvar import ClinVar
from GeneBreaker.src.transcript import Transcript
from GeneBreaker.src.api_helper import *


class clinvarTests(unittest.TestCase):
    SOX9_uid = get_all_transcripts("SOX9", "grch38")[0]["qualifiers"]["uid"]
    transcript = Transcript(SOX9_uid, "grch38")

    # test 1
    def test_clinvar(self):
        clinvar = {
            "TYPE": "CLINVAR",
            "REGION": "GENIC",
            "IMPACT": {"CLINVAR_ID": 507525},
            "ZYGOSITY": "HETEROZYGOUS"}
        clinvar = ClinVar(clinvar, self.transcript)
        row = clinvar.get_vcf_row()
        self.assertEqual(row["chrom"], "17")
        self.assertEqual(row["pos"], "72121409")
        self.assertEqual(row["ref"], "C")
        self.assertEqual(row["alt"], "T")
    
     # test 2
    def test_clinvar_wrong(self):
        clinvar = {
            "TYPE": "CLINVAR",
            "REGION": "GENIC",
            "IMPACT": {"CLINVAR_ID": 2000000000},
            "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            ClinVar(clinvar, self.transcript)

