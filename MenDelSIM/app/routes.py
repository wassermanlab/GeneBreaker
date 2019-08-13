from MenDelSIM.app import app
from flask import request, jsonify, render_template, url_for, redirect, send_from_directory
from werkzeug.exceptions import HTTPException, NotFound, BadRequest
from werkzeug.utils import secure_filename
import json
import sys, os
from MenDelSIM.src.variants import Variants

@app.route('/')
@app.route('/home')
def home():
    return render_template('home.html')

@app.route('/contact')
def contact():
    return "CONTACT"

@app.route('/design_variants', methods=['GET', 'POST'])
def design_variants():
    if request.method == 'POST':
        variant_data = request.json
        print(variant_data)
        return 'HELLO'

    # try: 
    #     variants_json = config
    #     variants = Variants(variants_json)
    #     variants.save_vcf_output(output)
    # except: 
    #     print ("Check that your config is formatted the correct way")


    # return "design_variants"








# def allowed_file(filename):
#     return '.' in filename and \
#            filename.rsplit('.', 1)[1].lower() in {'vcf'}

# @app.route('/design_family', methods=['GET', 'POST'])
# def design_family():
#     if request.method == 'GET':
#         return render_template('family.html')
#     elif request.method == 'POST':
#         print("POST")
#         # check if the post request has the file part
#         if 'file' not in request.files:
#             return "NO_FILE"
#         file = request.files['file']
#         # if user does not select file, browser also
#         if file and allowed_file(file.filename):
#             filename = secure_filename(file.filename)
#             file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
#             return redirect(url_for('design_family_file', filename=filename))
#         else: 
#             return "invalid file"

# @app.route('/design_family/<filename>')
# def design_family_file(filename):
#     return send_from_directory(app.config['UPLOAD_FOLDER'],filename)

# @app.route('/get_family')
# def get_family():
#     return "get_family"