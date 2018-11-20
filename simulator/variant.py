class Variant: 
    # assume var_template is of type dict already
    def __init__(self, var_template):
        try: 
            self.type = var_template["TYPE"]
            self.region = var_template["REGION"]
            self.impact = var_template["IMPACT"]
            self.location = var_template["LOCATION"]
            if self.location not in ["ANY", "NONE"]:
                self.location = long(var_template["LOCATION"])
            if self.type not in ['SNV', 'INDEL', 'CNV', 'SV', 'MEI', 'STR']:
                raise Exception('type not one of the recognized types: SNV, INDEL, CNV, SV, MEI, STR')
            if self.region not in ['CODING', 'UTR', 'INTRONIC', 'PROMOTER', 'ENHANCER']:
                raise Exception('region not one of the recognized regions: CODING, UTR, INTRONIC, PROMOTER, ENHANCER') 
            if type(self.location) is not long and self.location!="ANY" and self.location !="NONE":
                raise Exception("location is not valid")
        except:
            print("check that each variant follows the variant schema")
    

    def get_type(self):
        return self.type


    def get_region(self):
        return self.region
