import json
import os.path
from padmet.classes.padmetRef import PadmetRef

METACYC_PADMET = os.path.join('Data/metacyc_26.0_prot70.padmet')


def metacyc_file(input_padmet=METACYC_PADMET, output_name='MetaCyc', version='26_0'):
    output = output_name + version + '_classes.json'
    pref = PadmetRef(input_padmet)
    rels = pref.getAllRelation()
    classes = dict()
    for r in rels:
        if r.type == 'is_a_class':
            if r.id_in not in classes:
                classes[r.id_in] = set()
            classes[r.id_in].add(r.id_out)
    for c, p in classes.items():
        classes[c] = list(p)
    with open(output, 'w') as o:
        json.dump(fp=o, obj=classes)


metacyc_file()
