from GeneBreaker.app import app
from flask import request, jsonify, render_template, url_for, redirect, send_from_directory, send_file, make_response
from werkzeug.exceptions import HTTPException, NotFound, BadRequest
from werkzeug.utils import secure_filename
import json
import sys, os, time
from GeneBreaker.src.transcript import Transcript
from GeneBreaker.src.variants import Variants
from GeneBreaker.src.api_helper import *
from datetime import datetime
from GeneBreaker.src.file_outputs import make_ped, make_vcf

@app.route('/')
def test():
    return "TEST"

@app.route('/get_transcripts/<genome>/')
@app.route('/get_transcripts/<genome>/<name>')
def get_transcripts_api(genome, name = None):
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
        region_tuple = transcript.get_requested_region("GENIC")
    else:
        start = int(region.split(":")[1].split("-")[0])-1
        end = int(region.split(":")[1].split("-")[1])
        region_tuple = [(start, end)]
    chrom = transcript.get_chr()
    location = "overlapping"
    
    # get full genic or custom region
    if request.path.startswith("/get_str"):
        res = get_strs(region_tuple[0][0]+1, region_tuple[0][1], chrom, genome, location) #add 1 to start to make 1 based for api call
    elif request.path.startswith("/get_clinvar"):
        res = get_clinvars(region_tuple[0][0]+1, region_tuple[0][1], chrom, genome, location) #add 1 to start to make 1 based for api call
    else:
        res = get_cnvs(region_tuple[0][0]+1, region_tuple[0][1], chrom, genome, location) #add 1 to start to make 1 based for api call
  
    # cut if UTR, INTRONIC, or CODING 
    if region in ["UTR", "INTRONIC", "CODING"]:
        subregions = transcript.get_requested_region(region)
        res = extract_subregions(res, subregions, region)
    
    response = jsonify(res)
    return response

@app.route('/design_variants', methods=["POST"])
def design_variants():
    # if POST then return the variant file sent with the JSON file 
    if request.method == 'POST': 
        try: 
            response = Variants(request.json).get_variant_rows()
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
            dt_string = datetime.now().strftime("%d-%m-%Y_%H:%M:%S")
            if (filetype == "ped"):
                filename = "/tmp/simulator/"+ dt_string + ".ped"
                make_ped(request_payload, filename)
            elif (filetype == "vcf"): 
                filename = "/tmp/simulator/"+ dt_string + ".vcf"
                make_vcf(request_payload, filename)
            else: 
                Exception("incorrect type")
            return send_file(filename)
        except Exception as e:
            response = jsonify(e)
            return response



