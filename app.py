from flask import abort, Flask, flash, render_template, request, redirect, send_from_directory, make_response, session
import subprocess, identifiers, hgvs_conversion, sqlite3, re
from werkzeug.contrib.cache import FileSystemCache
from hashlib import sha256

# create the application object
APP = Flask(__name__, static_folder='static')
#APP.secret_key = "vrtoolsuite"
APP.config['UPLOAD_FOLDER'] = 'static/uploads'
cache = FileSystemCache('/tmp/vr-suite')

# GENERAL ROUTING
@APP.route('/')
def index():
    return render_template('overview.html', tabs=compileTabs("overview"))

@APP.route('/upload', methods=['POST'])
def upload():
    r = make_response(redirect(request.referrer))
    #check if the post request has the file part
    if 'file' not in request.files:
        return r
    file = request.files['file']
    # TODO Run VT tool to normalize variants in file, FASTA files already in referenceFASTA folder
    # subprocess.call(['vt', 'normalize', file.filename, '-r', seq.fa, '-o', 'normalized' + file.filename])
    read = file.read()
    vcf = read.decode()
    fileKey = sha256(read).hexdigest()
    #Set cookies to track the filename and the file key (hash of the file contents)
    r.set_cookie("fileKey", fileKey)
    r.set_cookie("filename", file.filename)
    #Store the uploaded VCF in the cache with the file key
    cache.set(fileKey, vcf, timeout=604800000)
    return r

# TOOL 1 ROUTING
@APP.route('/vcf', methods=['GET', 'POST'])
def tool1():
    fileKey = get_fileKey()
    #Check for the uploaded VCF in the cache
    if not fileKey or not cache.has(fileKey):
        return render_template('vcf_tool.html', tabs=compileTabs("vcf"), status="No file detected")
    vcf = cache.get(fileKey)
    if request.method == 'POST':
        """Checks to see if the transformed VCF already exists in the downloads folder, generates it if not.
            Also sends the file path for downloading the transformed VCF file."""
        #Check if transformed VCF exists in the cache
        if not cache.has("vr." + fileKey) and not cache.has(fileKey + ".json"):
            #Get the identifiers for the VCF
            build = request.form.get('build')
            sqs = getSQs(build);
            vr_vcf,vr_models = identifiers.getIdentifiers(vcf,sqs)
            #Inserts the identifiers back into the file for download
            print("Saving file to cache...")
            cache.set("vr." + fileKey, vr_vcf, timeout=604800000)
            cache.set(fileKey + '.json', vr_models, timeout=604800000)
            print('Done!')
        else:
            vr_vcf = cache.get("vr." + fileKey)
            vr_models = cache.get(fileKey + '.json')
        return render_template('vcf_tool.html', tabs=compileTabs("vcf"))
    return render_template('vcf_tool.html', tabs=compileTabs("vcf"))

@APP.route('/vcf-download')
def vcf_download():
    filename = "vr." + get_filename()
    fileKey = "vr." + get_fileKey()
    if not cache.has(fileKey):
        return render_template('vcf_tool.html', tabs=compileTabs("vcf"), status="Error finding file")
    response = make_response(cache.get(fileKey))
    response.headers.set("Content-Disposition", "attachment; filename=" + filename )
    return response

@APP.route('/json-download', methods=['GET', 'POST'])
def json_download():
    filename = get_filename()[0:-4] + ".json"
    fileKey = get_fileKey() + ".json"
    if not cache.has(fileKey):
        abort(400)
    response = make_response(cache.get(fileKey))
    response.headers.set("Content-Disposition", "attachment; filename=" + filename )
    return response

#TOOL 2 ROUTING
@APP.route('/hgvs', methods=['GET', 'POST'])
def tool2():
    if request.method == 'POST':
        """ Uses conversions.py to convert a HGVS string to a VMC JSON bundle and displays it along with the example JSON schema."""
        hgvs = request.form['hgvs'].strip()
        acc = hgvs.split(":")[0]
        sequence_id = "Sequence identifier not found for " + acc
        found = False
        with sqlite3.connect("VRSQ_DB.db") as db:
            sql = "SELECT VRSQ FROM accession WHERE ACCESSION=\'" + acc + "\'"
            cursor = db.cursor()
            cursor.execute(sql)
            rows = cursor.fetchall()
            for row in rows:
                found = True
                sequence_id = row[0]
        if found:
            vr_model = hgvs_conversion.from_hgvs(sequence_id,hgvs)
            return render_template('hgvs_tool.html', tabs=compileTabs("hgvs"), hgvs=hgvs, vr_model=vr_model)
        else:
            response = "No sequence identifier found for " + acc
            return render_template('hgvs_tool.html', tabs=compileTabs("hgvs"), hgvs=hgvs, vr_model=response)
    return render_template('hgvs_tool.html',tabs=compileTabs("hgvs"))

#TOOL 3 ROUTING
@APP.route('/digest', methods=['GET', 'POST'])
def tool3():
    if request.method == 'POST':
        build = request.form['build'].strip()
        chr = request.form['chr'].strip()
        start = request.form['start'].strip()
        end = request.form['end'].strip()
        state = request.form['state'].strip()
        with sqlite3.connect("VRSQ_DB.db") as db:
            sql = "SELECT VRSQ,ACCESSION FROM accession WHERE BUILD=\'" + build + "\' AND CHROM=\'" + chr + "\'"
            cursor = db.cursor()
            cursor.execute(sql)
            rows = cursor.fetchall()
            sequence_id = rows[0][0]
            acc = rows[0][1]
        sq,vsl,va,vr_model = identifiers.digestIdentifiers(sequence_id,acc,start,end,state)
        return render_template('digest_tool.html', tabs=compileTabs("digest"),start=start,end=end,state=state,sq=sq,vsl=vsl,va=va,vr_model=vr_model)
    else:
        return render_template('digest_tool.html',tabs=compileTabs("digest"))

def get_filename():
    """Returns the name of the uploaded VCF file"""
    return request.cookies.get("filename")

def get_fileKey():
    """Returns the key for the uploaded VCF file"""
    return request.cookies.get("fileKey")

def getSQs(build):
    #{1:(vrsq,accession)...}
    sqs = {}
    sql = 'select * from accession where BUILD=\"' + build + '\"'
    with sqlite3.connect("VRSQ_DB.db") as db:
        cursor = db.cursor()
        cursor.execute(sql)
        data = cursor.fetchall()
        for entry in data:
            sqs[entry[3]] = (entry[0],entry[1])
    return sqs

def compileTabs(option):
    if option == "vcf":
        return [(False,"Overview","/"), (True,"VCF","/vcf"), (False,"HGVS","/hgvs"), (False,"Digest","/digest")]
    elif option == "hgvs":
        return [(False,"Overview","/"), (False,"VCF","/vcf"), (True,"HGVS","/hgvs"), (False,"Digest","/digest")]
    elif option == "digest":
        return [(False,"Overview","/"), (False,"VCF","/vcf"), (False,"HGVS","/hgvs"), (True,"Digest","/digest")]
    else:
        return [(True,"Overview","/"), (False,"VCF","/vcf"), (False,"HGVS","/hgvs"), (False,"Digest","/digest")]

# start the server with the 'run()' method
if __name__ == '__main__':
    APP.run(debug=True, host="0.0.0.0",port=8080)
