from MenDelSIM.app import app
from flask import request, jsonify, render_template, url_for, redirect, send_from_directory, send_file
from werkzeug.exceptions import HTTPException, NotFound, BadRequest
from werkzeug.utils import secure_filename
import json
import sys, os
from MenDelSIM.src.transcript import Transcript
from MenDelSIM.src.variants import Variants
from MenDelSIM.src.api_helper import *
from datetime import datetime

@app.route('/get_transcripts/<genome>/<name>')
def get_transcripts_api(genome,name):
    # TODO: remove cross origin
    response = jsonify(get_all_transcripts(name, genome))
    response.headers.add('Access-Control-Allow-Origin', '*')
    return response

@app.route('/get_str/<genome>/<transcript_uid>/<region>')
@app.route('/get_clinvar/<genome>/<transcript_uid>/<region>')
@app.route('/get_clingen/<genome>/<transcript_uid>/<region>')
def get_clingen_clinvar_str_api(genome, transcript_uid, region):
    transcript = Transcript(transcript_uid, genome)
    res = []
    if region in ["UTR", "INTRONIC", "GENIC", "CODING"]:
        region = transcript.get_requested_region(region)
    else:
        start = int(region.split(":")[1].split("-")[0])-1
        end = int(region.split(":")[1].split("-")[1])
        region = [(start, end)]
    location = "overlapping"
    if region in ["UTR", "INTRONIC"]:
        location = "within"
    for r in region:
        if request.path.startswith("/get_str"):
            res = res + get_strs(r[0]+1, r[1], transcript.get_chr(), genome, location) #add 1 to start to make 1 based for api call
        elif request.path.startswith("/get_clinvar"):
            res = res + get_clinvars(r[0]+1, r[1], transcript.get_chr(), genome, location) #add 1 to start to make 1 based for api call
        else:
            res = res + get_cnvs(r[0]+1, r[1], transcript.get_chr(), genome, location) #add 1 to start to make 1 based for api call
    response = jsonify(res)
    response.headers.add('Access-Control-Allow-Origin', '*')
    return response

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
    else: # render the template to make variants \
        return render_template('variants.html')

@app.route('/design_family', methods=["GET", "POST"])
def design_family():
    return "family"