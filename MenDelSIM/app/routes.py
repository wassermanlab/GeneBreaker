from MenDelSIM.app import app
from flask import request, jsonify, render_template, url_for, redirect, send_from_directory, send_file
from werkzeug.exceptions import HTTPException, NotFound, BadRequest
from werkzeug.utils import secure_filename
import json
import sys, os
from MenDelSIM.src.variants import Variants
from datetime import datetime

@app.route('/')
@app.route('/home')
def home():
    return render_template('home.html')

@app.route('/contact')
def contact():
    return render_template('contact.html')

@app.route('/design_variants', methods=["GET", "POST"])
def design_variants():
    # if POST then return the variant file sent with the JSON file 
    if request.method == 'POST': 
        if not os.path.exists('/tmp/simulator'):
            os.makedirs('/tmp/simulator')
        try: 
            now = datetime.now()
            dt_string = now.strftime("%d-%m-%Y_%H:%M:%S")
            Variants(request.json).save_vcf_output("/tmp/simulator/"+ dt_string + ".vcf")
            return send_file("/tmp/simulator/"+ dt_string + ".vcf", attachment_filename=dt_string + ".vcf")
        except Exception as e:
            return BadRequest("cannot produce requested variant: " + str(e))
    else: # render the template to make variants  
        return render_template('variants.html')

@app.route('/design_family', methods=["GET", "POST"])
def design_family():
    return "family"