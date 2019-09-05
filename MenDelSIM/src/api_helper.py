import requests
import json
host = 'http://127.0.0.1:5000'

def get_all_results(request_url):
    # http://127.0.0.1:5000/api/v1/hg38/chroms
    results = []
    next = True
    url = request_url
    while(next):
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

def get_all_genenames(genome):
    url = host + '/api/v1/' + genome + '/genes/symbols'
    return get_all_results(url)

def get_all_transcripts(gene_name, genome):
    url = host + '/api/v1/' + genome + '/genes?names='+gene_name
    return get_all_results(url)

def get_str(start, end, chrom, genome):
    #takes in 1 based
    url = host + '/api/v1/' + genome + '/short_tandem_repeats?start='+ str(start) + '&end='+ str(end) + '&chrom=' + chrom + '&location=exact'
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

def get_clinvars(start, end, chrom, genome, location):
    #takes in 1 based
    url = host + '/api/v1/' + genome + '/short_tandem_repeats?start='+ str(start) + '&end='+ str(end) + '&chrom=' + chrom + '&location=' +location
    res = get_all_results(url)
    if res != False: 
        return res
    return []

def get_cnvs(start, end, chrom, genome, location):
    #takes in 1 based
    url = host + '/api/v1/' + genome + '/copy_number_variants?start='+ str(start) + '&end='+ str(end) + '&chrom=' + chrom + '&location=' +location
    print(url)
    res = get_all_results(url)
    if res != False: 
        return res
    return []

def get_transcript(uid, genome):
    url = host + '/api/v1/' + genome + '/genes?uids='+ str(uid)
    res = get_all_results(url)
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

# print(get_chromosome("chr3", "hg19"))
# t = get_transcript("69540","hg38")
# print(t)        
# print(len(get_all_genenames("hg38")))
# print(len(get_all_results("http://127.0.0.1:5000/api/v1/hg38/chroms")))
