from flask import render_template, flash, redirect, url_for, jsonify, request
from . import app
from GUD.ORM import Gene
from . import establish_GUD_session
session = establish_GUD_session()
# from variant_gen import get_genes

@app.route('/')
@app.route('/home')
def home():
    return render_template('home.html', title='Home')

@app.route('/variant_gen', methods=['GET'])
def variant_gen():
    return render_template('variant_gen.html', title='Variant Generator')

@app.route('/_get_transcript', methods=['GET'])
def get_transcript():
    gene_name = request.args.get('gene_name', "", type=str)
    gene = Gene()
    genes = gene.select_by_name(session, gene_name, True)
    response = []
    if len(genes) > 0:
        for g in genes:
            response.append({'uid': int(g.qualifiers["uid"]),
                            'accession': str(g.qualifiers["name"]),
                            'name': str(g.qualifiers["name2"]),
                            'chrom': str(g.chrom)})
    return jsonify(response)

