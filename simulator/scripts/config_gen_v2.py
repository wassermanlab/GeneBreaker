import json
import re
import os
import sys
from GUD.ORM import Gene
from GUD.ORM import ShortTandemRepeat
from GUD.ORM import ClinVar
from GUD.ORM import Chrom
from GUD.ORM import CNV
from GUD.ORM.genomic_feature import GenomicFeature
from simulator.src.transcript import Transcript
from . import establish_GUD_session
from prompt_toolkit.shortcuts import message_dialog, input_dialog, yes_no_dialog, button_dialog

session = establish_GUD_session()

def get_main_dict():
    # ## file name
    # file_name = None
    # while file_name is None:
    #     file_name = input_dialog(
    #         title='Configuration File Creator',
    #         text='Enter the name of your config file (no spaces allowed):',
    #         ok_text = "Enter",
    #         cancel_text = "Rename")
        
    #     if file_name is not None: ## check that filename has no spaces
    #         if " " in file_name or file_name is "":
    #             message_dialog(
    #                 title='Configuration File Creator',
    #                 text='Empty file names and spaces in file names are not permitted, try again.')
    #             file_name = None      
    # ## sex
    # sex = button_dialog(
    #     title='Button dialog example',
    #     text='What sex would you like your proband to be?',
    #     buttons=[
    #         ('XY Male', "XY-MALE"),
    #         ('XX Female', "XX-FEMALE")
    #     ])
    # gene 
    gene_symbol = None
    while gene_symbol is None:
        gene_symbol = input_dialog(
            title='Configuration File Creator',
            text='Enter the name of your gene:',
            ok_text = "Enter",
            cancel_text = "Re-Enter")
        if gene_symbol is not None and gene_symbol is not "": ## check that gene is not empty
            gene = Gene()
            genes = gene.select_by_name(session, gene_symbol)
            if genes:
                pass
            else:
                message_dialog(
                    title='Configuration File Creator',
                    text='Gene name does not exist, try again.')
                gene_symbol = None
               

def main():
    satisfied = False
    while satisfied is False:
        config_dict = get_main_dict()
        satisfied = yes_no_dialog(
            title='Configuration File Creator',
            text='Are you satisfied with your config?') ## todo add full file 


if __name__ == "__main__":
    main()