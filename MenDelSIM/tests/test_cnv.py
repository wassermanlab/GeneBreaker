import unittest
from MenDelSIM.src.copy_number_variant import CopyNumberVariant
from MenDelSIM.src.transcript import Transcript
from MenDelSIM.src.api_helper import *

class IndelCreationTests(unittest.TestCase):
    XKR8_uid = get_all_transcripts("XKR8", "hg38")[0]["qualifiers"]["uid"]
    transcript = Transcript(XKR8_uid, "hg38")
    # test 1
    def test_copy_number_error(self):
        print(self.transcript.get_regionmap())
        # cnv = {
        #     "TYPE": "CNV",
        #     "REGION": "GENIC",
        #     "IMPACT": {
        #         "CHROM": "chr1",
        #         "START": 27959462,
        #         "END": 27959762,
        #         "COPY_CHANGE": 1,  
        #     },
        #     "ZYGOSITY": "HETEROZYGOUS"}
        # with self.assertRaises(ValueError):
        #     CopyNumberVariant(cnv, self.transcript)
        # cnv["IMPACT"]["COPY_CHANGE"] = 0
        # with self.assertRaises(ValueError):
        #     CopyNumberVariant(cnv, self.transcript)
        # cnv["IMPACT"]["COPY_CHANGE"] = -3
        # with self.assertRaises(ValueError):
        #     CopyNumberVariant(cnv, self.transcript)