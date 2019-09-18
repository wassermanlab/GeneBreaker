from MenDelSIM.app import app
from flask import request, jsonify, render_template, url_for, redirect, send_from_directory, send_file, make_response
from werkzeug.exceptions import HTTPException, NotFound, BadRequest
from werkzeug.utils import secure_filename
import json
import sys, os
from MenDelSIM.src.transcript import Transcript
from MenDelSIM.src.variants import Variants
from MenDelSIM.src.api_helper import *
from datetime import datetime
from MenDelSIM.src.file_outputs import make_ped, make_vcf

@app.route('/get_transcripts/<genome>/')
@app.route('/get_transcripts/<genome>/<name>')
def get_transcripts_api(genome,name = None):
    # TODO: remove cross origin
    if (name is None):
        return jsonify([])
    response = get_all_transcripts(name, genome)
    if response is False: 
        return jsonify([])
    response = jsonify(response)
    return response

@app.route('/get_str/<genome>/<transcript_uid>/')
@app.route('/get_clinvar/<genome>/<transcript_uid>/')
@app.route('/get_clingen/<genome>/<transcript_uid>/')
@app.route('/get_str/<genome>/<transcript_uid>/<region>')
@app.route('/get_clinvar/<genome>/<transcript_uid>/<region>')
@app.route('/get_clingen/<genome>/<transcript_uid>/<region>')
def get_clingen_clinvar_str_api(genome, transcript_uid, region=None):
    if (region is None ):
        response = jsonify([])
        return response
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
    return response

@app.route('/design_variants', methods=["POST"])
def design_variants():
    # if POST then return the variant file sent with the JSON file 
    if request.method == 'POST': 
        try: 
            response = Variants(request.json).save_vcf_output()
            response = jsonify(response)
            return response

        except Exception as e:
            response = jsonify({"error": str(e)})
            return response

@app.route('/get_file', methods=["POST"])
def get_file():
    # takes in a json file like this {var1: {}, var2: {}, family: {}}
    if request.method == 'POST': 
        if not os.path.exists('/tmp/simulator'):
            os.makedirs('/tmp/simulator')
        try: 
            filetype = request.args.get('filetype')
            request_payload = request.json 
            print(request_payload)
            dt_string = datetime.now().strftime("%d-%m-%Y_%H:%M:%S")
            if (filetype == "ped"):
                filename = dt_string + ".ped"
                make_ped(request_payload, filename)
            elif (filetype == "vcf"): 
                filename = dt_string + ".vcf"
                make_vcf(request_payload, filename)
            else: 
                Exception("incorrect type")
            return send_file("/tmp/simulator/"+ filename)
        except Exception as e:
            response = jsonify(e)
            return response



