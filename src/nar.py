import os
import json 

BASEDIR = os.path.dirname(os.path.abspath(__file__))
PATH    = BASEDIR + '/REPR.json'


def get_repr(key:'str'):
    with open('{}/{}/{}.json'.format(BASEDIR, 'nuclrepr', key), 'r') as file:
        d = json.loads(file.read())
    return d
