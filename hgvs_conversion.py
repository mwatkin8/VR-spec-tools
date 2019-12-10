import hgvs.parser, identifiers, json

def from_hgvs(SQ,hgvs_expr):
    parser = hgvs.parser.Parser()
    sv = parser.parse_hgvs_variant(hgvs_expr)

    if isinstance(sv.posedit.pos, hgvs.location.BaseOffsetInterval):
        if sv.posedit.pos.start.is_intronic or sv.posedit.pos.end.is_intronic:
            raise ValueError("Intronic HGVS variants are not supported ({sv.posedit})")

    if sv.posedit.edit.type == 'ins':
        start=sv.posedit.pos.start.base
        end=sv.posedit.pos.start.base
        state = sv.posedit.edit.alt

    elif sv.posedit.edit.type in ('sub', 'del', 'delins', 'identity'):
        start=sv.posedit.pos.start.base - 1
        end=sv.posedit.pos.end.base
        state = sv.posedit.edit.alt or ''

    else:
        raise ValueError(f"HGVS variant type {sv.posedit.edit.type} is unsupported")

    location = identifiers.getVSL(SQ,start,end)
    vsl = 'ga4gh:VSL.' + identifiers.digest(bytes(json.dumps(location), 'utf-8'))
    allele = identifiers.getVA(vsl,state)
    va = 'ga4gh:VA.' + identifiers.digest(bytes(json.dumps(allele), 'utf-8'))
    model = identifiers.assembleJSON([SQ,sv.ac],location,allele,va)

    return json.dumps(model, ensure_ascii=False, indent=4)
