import base64,hashlib,sqlite3,json
from datetime import datetime
from collections import OrderedDict

def getIdentifiers(vcf,sqs):
    ids = []
    models = {"alleles":[]}
    print("Getting identifiers...")
    for line in vcf.split('\n'):
        try:
            if line != "":
                if line[0] == '#':
                    continue
                else:
                    ll = line.split('\t')
                    chr = ll[0]
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
                    ids.append((sq,vsl,va))
                    model = assembleJSON(sqs[chr],vsl_model,va_model,va)
                    models["alleles"].append(model)
        except:
            continue

    return ids,json.dumps(models, ensure_ascii=False, indent=4)

def digestIdentifiers(sq,acc,start,end,state):
    vsl_model = getVSL(sq,start,end)
    vsl = 'ga4gh:VSL.' + digest(bytes(json.dumps(vsl_model), 'utf-8'))
    va_model = getVA(vsl,state)
    va = 'ga4gh:VA.' + digest(bytes(json.dumps(va_model), 'utf-8'))
    model = assembleJSON([sq,acc],vsl_model,va_model,va)
    return sq,vsl,va,json.dumps(model, ensure_ascii=False, indent=4)

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

def digest(blob, n=24):
    d = hashlib.sha512(blob).digest()
    result = base64.urlsafe_b64encode(d[:n]).decode("ascii")
    return result

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

def insertIdentifiers(ids, vcf):
    final = ""
    idx = 0
    info = 0
    for line in vcf.split('\n'):
        if line != "":
            if line[1] == '#':
                final += line + '\n'
            elif line[0] == '#':
                headers = line.split('\t')
                for i in range(0, len(headers)):
                    if headers[i] == 'INFO':
                        info = i + 1
                final = final + """##INFO=<ID=VRSQ,Number=1,Type=String,Description="VR Sequence identifier">""" + '\n' + \
                        """##INFO=<ID=VRVSL,Number=1,Type=String,Description="VR Sequence Location identifier">""" + '\n' + \
                        """##INFO=<ID=VRVA,Number=1,Type=String,Description="VR Allele identifier">\n""" + line
            else:
                line_list = line.split('\t')
                id_line = ";VRSQ=" + ids[idx][0] + ";VRVSL=" + ids[idx][1] + ";VA=" + ids[idx][2]
                final = final + '\n' + '\t'.join(line_list[0:info]) + id_line  + '\t'.join(line_list[info:])
                idx += 1
    return final
