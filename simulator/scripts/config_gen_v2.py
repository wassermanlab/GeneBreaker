import json, re, os, sys, cmd
from GUD.ORM import Gene
from GUD.ORM import ShortTandemRepeat
from GUD.ORM import ClinVar
from GUD.ORM import Chrom
from GUD.ORM import CNV
from GUD.ORM.genomic_feature import GenomicFeature
from simulator.src.transcript import Transcript
from . import establish_GUD_session
# from prompt_toolkit.shortcuts import message_dialog, input_dialog, yes_no_dialog, button_dialog
 
# def get_main_dict(cmd.Cmd):
#     intro   = "=====Configuration File Generator====="
#     prompt  = "(config_gen)"

#     # ## file name
#     # file_name = None
#     # while file_name is None:
#     #     file_name = input_dialog(
#     #         title='Configuration File Creator',
#     #         text='Enter the name of your config file (no spaces allowed):',
#     #         ok_text = "Enter",
#     #         cancel_text = "Rename")
        
#     #     if file_name is not None: ## check that filename has no spaces
#     #         if " " in file_name or file_name is "":
#     #             message_dialog(
#     #                 title='Configuration File Creator',
#     #                 text='Empty file names and spaces in file names are not permitted, try again.')
#     #             file_name = None      
#     # ## sex
#     # sex = button_dialog(
#     #     title='Button dialog example',
#     #     text='What sex would you like your proband to be?',
#     #     buttons=[
#     #         ('XY Male', 'XY-MALE'),
#     #         ('XX Female', 'XX-FEMALE')
#     #     ])
#     # gene 
#     # gene_symbol = None
#     # while gene_symbol is None:
#     #     gene_symbol = input_dialog(
#     #         title='Configuration File Creator',
#     #         text='Enter the name of your gene:',
#     #         ok_text = "Enter",
#     #         cancel_text = "Re-Enter")
#     #     if gene_symbol is not None and gene_symbol is not "": ## check that gene is not empty
#     #         gene = Gene()
#     #         genes = gene.select_by_name(session, gene_symbol, True)
#     #         if genes:
#     #             gene_selections = []
#     #             gene_texts = ""
#     #             for g in genes:
#     #                 button_text = "Start: {}\tEnd: {}\tAccession: {}\tUID: {}\n".format(g.start, g.end,\
#     #                      g.qualifiers["name"], g.qualifiers["uid"])
#     #                 gene_selections.append((str(g.qualifiers["uid"]), str(g.qualifiers["uid"])))
#     #                 gene_texts = gene_texts + button_text
#     #             gene_symbol = button_dialog(
#     #                 title='Button dialog example',
#     #                 text='Select the UID of the transcript of {} you would like?\n'.format(gene_symbol) + gene_texts,
#     #                 buttons=gene_selections)
#     #         else:
#     #             message_dialog(
#     #                 title='Configuration File Creator',
#     #                 text='Gene name does not exist, try again.')
#     #             gene_symbol = None
               

# def main(cmd.Cmd):
#     config_dict = get_main_dict()


# if __name__ == "__main__":
#     main()