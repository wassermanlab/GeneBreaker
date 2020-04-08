import requests
import json
import configparser
import time
config = configparser.ConfigParser()
config.read('GeneBreaker/src/config.ini')
host = config.get('DEFAULT', 'GUD_HOST')

def get_all_results(request_url):
    time.sleep(1)
    results = []
    next = True
    url = request_url
    # in this while loop wait 1 sec
    counter = 0 
    while(next):
        if (counter == 4):
            counter = 0 
            time.sleep(1)
        else: 
            counter = counter + 1
        req = requests.get(url)
        if req.status_code == 200:
            res = json.loads(req.content.decode('utf-8'))
            results = results + res['results']
            if ('next' in res):
                url = res['next']
            else: 
                next = False
        else:
            return False
    return results

def get_all_transcripts(gene_name, genome):
    url = host + '/api/v1/' + genome + '/genes?names='+gene_name
    print(url)
    return get_all_results(url)

def get_str(genome, uid):
    #takes in 1 based
    url = host + '/api/v1/' + genome + '/short_tandem_repeats?uids=' + str(uid)
    res = get_all_results(url)
    if res != False: 
        return res[0]
    return []

def get_strs(start, end, chrom, genome, location):
    #takes in 1 based
    url = host + '/api/v1/' + genome + '/short_tandem_repeats?start='+ str(start) + '&end='+ str(end) + '&chrom=' + chrom + '&location=' + location
    res = get_all_results(url)
    if res != False:
        return res
    return []

def get_clinvar(genome, uid):
    #takes in 1 based
    url = host + '/api/v1/' + genome + '/clinvar?clinvar_ids=' + str(uid)
    res = get_all_results(url)
    if res != False: 
        return res[0]
    return []

def get_clinvars(start, end, chrom, genome, location):
    #takes in 1 based
    url = host + '/api/v1/' + genome + '/clinvar?start='+ str(start) + '&end='+ str(end) + '&chrom=' + chrom + '&location=' +location
    res = get_all_results(url)
    if res != False: 
        return res
    return []

def get_cnvs(start, end, chrom, genome, location):
    #takes in 1 based
    url = host + '/api/v1/' + genome + '/copy_number_variants?start='+ str(start) + '&end='+ str(end) + '&chrom=' + chrom + '&location=' +location
    res = get_all_results(url)
    if res != False: 
        return res
    return []

def get_transcript(uid, genome):
    url = host + '/api/v1/' + genome + '/genes?uids='+ str(uid)
    res = get_all_results(url)
    print (res)
    if res != False: 
        return res[0]
    return []

def get_chromosome(chr, genome): 
    url = host + '/api/v1/' + genome + '/chroms'
    res = get_all_results(url)
    for i in res:
        if i["chrom"] == chr:
            return i
    return False
