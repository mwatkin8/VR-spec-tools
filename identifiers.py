import base64,hashlib,json
from multiprocessing import Pool

sqs = []
info = None;

def processLine(line):
    final = ''
    line = line.strip('\n')
    try:
        if line != "":
            if line[1] == '#':
                return [line,{}]
            elif line[0] == '#':
                out = """##INFO=<ID=VRSQ,Number=1,Type=String,Description="VR Sequence identifier">""" + '\n' + \
                        """##INFO=<ID=VRVSL,Number=1,Type=String,Description="VR Sequence Location identifier">""" + '\n' + \
                        """##INFO=<ID=VRVA,Number=1,Type=String,Description="VR Allele identifier">\n""" + line
                return [out,{}]
            else:
                ll = line.split('\t')
                chr = ll[0]
                if chr == 'MT':
                    return ['',{}]
                start = int(ll[1]) - 1
                ref = ll[3]
                alt = ll[4]
                end = start + len(ll[3])
                state = alt
                #indel
                if len(alt) > 1 and len(ref) > 1:
                    end = start + len(ref) - 1
                    state = alt
                else:
                    #ins
                    if len(alt) > 1:
                        state = alt[1:]
                        end = start + 1
                    #del
                    if len(ref) > 1:
                        start = start + 1
                        end = start + len(ref[1:])
                        state = ''
                sq = sqs[chr][0]
                vsl_model = getVSL(sq,start,end)
                vsl = 'ga4gh:VSL.' + digest(bytes(json.dumps(vsl_model), 'utf-8'))
                va_model = getVA(vsl,state)
                va = 'ga4gh:VA.' + digest(bytes(json.dumps(va_model), 'utf-8'))
                id_line = ";VRSQ=" + sq + ";VRVSL=" + vsl + ";VA=" + va
                model = assembleJSON(sqs[chr],vsl_model,va_model,va)
                l = '\t'.join(ll[0:info]) + id_line  + '\t'.join(ll[info:])
                return [l,model]
    except Exception as e:
        print(e)
        return ''

def getIdentifiers(vcf,sq_list):
    global sqs
    global info
    #Save the list of SQ identifiers as a global (accessible for the Pool)
    sqs = sq_list
    #Save the position of the INFO field as a global (accessible for the Pool)
    h = vcf.split('#CHROM')
    headers = h[1].split('\n')
    info_pos = 1
    for h in headers[0].split('\t'):
        if h == 'INFO':
            break
        else:
            info_pos += 1
    info = info_pos
    print("Getting identifiers...")
    #Send the file to the Pool for processing
    lines = vcf.split('\n')
    combined = Pool().map(processLine,lines)
    combined.remove(None)
    print('Done!')
    ids = []
    m = []
    print('Assembling outputs...')
    for entry in combined:
        ids.append(entry[0])
        m.append(entry[1])
    ids.remove('')
    m = list(filter(None,m))
    models = {"alleles":m}
    print('Done!')
    return '\n'.join(ids),json.dumps(models, ensure_ascii=False, indent=4)

def digest(blob, n=24):
    d = hashlib.sha512(blob).digest()
    result = base64.urlsafe_b64encode(d[:n]).decode("ascii")
    return result

def assembleJSON(sqs,vsl_model,va_model,va):
    seq_id = "refseq:" + sqs[1]
    t = {
      "_id":va,
      "location": {
        "interval": {
          "end": vsl_model["interval"]["end"],
          "start": vsl_model["interval"]["start"],
          "type": "SimpleInterval"
        },
        "sequence_id": seq_id,
        "type": "SequenceLocation"
      },
      "state": {
        "sequence": va_model["state"]["sequence"],
        "type": "SequenceState"
      },
      "type": "Allele"
    }
    return t

def getVSL(sq,start,end):
    seq = sq.split('.')[1]
    t = {
        "interval": {
            "end": start,
            "start": end,
            "type": "SimpleInterval"
        },
        "sequence_id": sq,
        "type": "SequenceLocation"
    }
    return t

def getVA(vsl,state):
    vsl = vsl.split('.')[1]
    t = {
        "location":vsl,
        "state":{
            "sequence":state,
            "type":"SequenceState"
        },
        "type":"Allele"
    }
    return t

def digestIdentifiers(sq,acc,start,end,state):
    vsl_model = getVSL(sq,start,end)
    vsl = 'ga4gh:VSL.' + digest(bytes(json.dumps(vsl_model), 'utf-8'))
    va_model = getVA(vsl,state)
    va = 'ga4gh:VA.' + digest(bytes(json.dumps(va_model), 'utf-8'))
    model = assembleJSON([sq,acc],vsl_model,va_model,va)
    return sq,vsl,va,json.dumps(model, ensure_ascii=False, indent=4)
