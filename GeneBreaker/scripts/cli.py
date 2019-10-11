#! /usr/bin/env python

import argparse
import datetime
import json
from GeneBreaker.src.variants import Variants

def parse_arguments():
    """Parses inputted arguments as described"""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-c', '--config', help="config file", type=str, default= "config.json")
    parser.add_argument(
        '-o', '--output', help='output vcf',
        type=str, default="output.vcf")
    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()
    config = args.config
    output = args.output
    # parse config so that it has the right output
    try: 
        variants_json = config
        variants = Variants(variants_json)
        variants.save_vcf_output(output)
    except: 
        print ("Check that your config is formatted the correct way")

if __name__ == "__main__":
    main()