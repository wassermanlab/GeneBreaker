import gene

class Variant: 
# "VAR1": {
#         "TYPE": "SNV",
#         "REGION": "INTRONIC",
#         "IMPACT": -8332
#     }
    # assume var_template is of type dict already
    def __init__(self, var_template):
        try: 
            self.type = var_template["TYPE"]
            self.region = var_template["REGION"]
            self.impact = var_template["IMPACT"]
        except:
            print("check that each variant follows the variant schema")