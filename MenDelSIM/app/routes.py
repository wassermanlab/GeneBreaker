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
def get_strs_api(genome, transcript_uid, region):
    transcript = Transcript(transcript_uid, genome)
    region = transcript.get_requested_region(region)
    strs = []
    for r in region:
        strs = strs + get_strs(r[0]+1, r[1], transcript.get_chr(), genome, "within") #add 1 to start to make 1 based for api call
    response = jsonify(strs)
    response.headers.add('Access-Control-Allow-Origin', '*')
    return response

@app.route('/get_clinvar/<genome>/<transcript_uid>/<region>')
def get_clinvar_api(genome, transcript_uid, region):
    transcript = Transcript(transcript_uid, genome)
    region = transcript.get_requested_region(region)
    clinvar = []
    for r in region:
        clinvar = clinvar + get_clinvars(r[0]+1, r[1], transcript.get_chr(), genome, "within") #add 1 to start to make 1 based for api call
    response = jsonify(clinvar)
    response.headers.add('Access-Control-Allow-Origin', '*')
    return response

@app.route('/get_clingen/<genome>/<transcript_uid>/<region>')
def get_clingen_api(genome, transcript_uid, region):
    transcript = Transcript(transcript_uid, genome)
    region = transcript.get_requested_region(region)
    clingen = []
    for r in region:
        clingen = clingen + get_cnvs(r[0]+1, r[1], transcript.get_chr(), genome, "within") #add 1 to start to make 1 based for api call
    response = jsonify(clingen)
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