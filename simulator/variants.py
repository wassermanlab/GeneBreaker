class Variants: 

    def __init__(self, gene, inheretance, var1_template, var2_template):
        """ 
        creates variants object which has gene, inheretance, var1, var2 
        """
        self.gene = gene # may have to expand to phenotype 
        self.inheretance = inheretance
        self.var1_template = var1_template
        self.var2_template = var2_template


    def template_2_variant(): 
        """ turns variant template into variant, returning a dict with the 
        following keys: chrom, pos, id, ref, alt """
        return {"chrom": 1, 
                "pos": 1, 
                "ID": "ABC",
                "REF": "A"
                "ALT": "C"}


    def save_vcf_output(file_name, trio)
        """ saves a vcf output of the variants
        if  trio == false then it just saves the child
        if trio == true it saves a child and 2 parents"""
        return False