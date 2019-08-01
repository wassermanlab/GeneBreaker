from flask import render_template, flash, redirect, url_for, jsonify, request
from . import app

@app.route('/')
@app.route('/home')
def home():
    return render_template('home.html')

@app.route('/variants')
def variants():
    return render_template('variants.html')

@app.route('/family')
def family():
    var1 = request.args.get('var1', default=None) 
    var2 = request.args.get('var2', default=None) 
    sex = request.args.get('sex', default=None) 
    return render_template('family.html', sex = sex, var1=var1, var2=var2)





















# @app.route('/variant_gen2', methods=['GET'])
# def variant_gen2():
#     return render_template('variant_gen.html', title='Variant Generator')

# @app.route('/variant_gen', methods=['GET'])
# def variant_gen():
#     return render_template('variants.html', title='Variant Generator')

# @app.route('/variant_gen/var1', methods=['GET'])
# def var1():
#     gene_uid = request.args.get('gene_uid', type=str)
#     chrom = request.args.get('chrom', type=str)
#     sex = request.args.get('sex', type=str)
#     return render_template('variant_specifics.html',  
#     title='Variant Generator', gene_uid=gene_uid, chrom=chrom, sex=sex)

# @app.route('/_get_transcript', methods=['GET'])
# def get_transcript():
#     session = establish_GUD_session()
#     gene_name = request.args.get('gene_name', "", type=str)
#     gene = Gene()
#     genes = gene.select_by_name(session, gene_name, True)
#     session.close()
#     response = []
#     if len(genes) > 0:
#         for g in genes:
#             response.append({'uid': int(g.qualifiers["uid"]),
#                             'accession': str(g.qualifiers["name"]),
#                             'name': str(g.qualifiers["name2"]),
#                             'chrom': str(g.chrom)})
#     return jsonify(response)

