import unittest
from unittest.mock import patch
import io
import sys
from functools import wraps
from ontosunburst.ontology import *
from ontosunburst.ontosunburst import *

# ==================================================================================================
# GLOBAL
# ==================================================================================================

CH_URL = 'http://localhost:3030/chebi/'
CH_LST = ['38028', '28604', '85146']
REF_CH = ['38028', '28604', '85146',
          '23066', '27803', '37565',
          '58215', '79983', '42639']


# ==================================================================================================
# FUNCTIONS UTILS
# ==================================================================================================

def save_fig_json(fig, file):
    fig = fig.to_dict()
    with open(file, 'w') as f:
        json.dump(fig, f)


def sort_fig_dict_lst(fig):
    """
    {data: [{hovertext: [str],
             ids: [str],
             labels: [str],
             marker: {colors: [int]},
             parents: [str],
             values: [int]}]}
    """
    lst_keys = ['hovertext', 'ids', 'labels', 'parents', 'values']
    for k in lst_keys:
        fig['data'][0][k] = sorted(fig['data'][0][k])
    fig['data'][0]['marker']['colors'] = sorted(fig['data'][0]['marker']['colors'])
    return fig


def fig_to_lines(fig_dict):
    lst_keys = ['hovertext', 'ids', 'labels', 'parents', 'values']
    lines = set()
    for i in range(len(fig_dict['data'][0]['ids'])):
        line = tuple(fig_dict['data'][0][k][i] for k in lst_keys)
        line = line + (str(fig_dict['data'][0]['marker']['colors'][i]),)
        lines.add(line)
    return lines


def are_fig_dict_equals(fig1, fig2_file):
    fig1 = fig1.to_dict()
    fig1_l = fig_to_lines(fig1)
    fig1 = json.dumps(sort_fig_dict_lst(fig1), sort_keys=True)
    with open(fig2_file, 'r') as f:
        fig2 = json.load(f)
        fig2_l = fig_to_lines(fig2)
        fig2 = json.dumps(sort_fig_dict_lst(fig2), sort_keys=True)
    return (fig1 == fig2) and (fig1_l == fig2_l)


def dicts_with_sorted_lists_equal(dict1, dict2):
    if dict1.keys() != dict2.keys():
        return False
    for key in dict1:
        if sorted(dict1[key]) != sorted(dict2[key]):
            return False
    return True


def test_for(func):
    def decorator(test_func):
        @wraps(test_func)
        def wrapper(*args, **kwargs):
            return test_func(*args, **kwargs)
        wrapper._test_for = func
        return wrapper
    return decorator


class DualWriter(io.StringIO):
    def __init__(self, original_stdout):
        super().__init__()
        self.original_stdout = original_stdout

    def write(self, s):
        super().write(s)
        self.original_stdout.write(s)


# ==================================================================================================
# UNIT TESTS
# ==================================================================================================

# TEST ONTOLOGY : CLASSES EXTRACTION
# --------------------------------------------------------------------------------------------------
class TestChEBIClassesExtraction(unittest.TestCase):
    @test_for(extract_chebi_roles)
    @patch('sys.stdout', new_callable=lambda: DualWriter(sys.stdout))
    def test_extract_chebi_roles(self, mock_stdout):
        all_roles, d_roles_ontology, names = extract_chebi_roles(CH_LST, CH_URL)
        output = mock_stdout.getvalue().strip()
        wanted_ontology = {'35703': ['24432'], '24432': ['role'], '38028': ['35703'],
                           '75767': ['75763'], '75763': ['25212'], '25212': ['52206'],
                           '52206': ['24432'], '84735': ['75763'], '76507': ['25212'],
                           '76924': ['75763'], '28604': ['75767', '76507', '84735', '76924'],
                           '50503': ['23888'], '23888': ['52217'], '52217': ['33232'],
                           '33232': ['role'], '63047': ['63046', '64047'],
                           '64047': ['33232', '78295'], '78295': ['52211'], '52211': ['24432'],
                           '63046': ['51086'], '51086': ['role'], '85146': ['50503', '63047']}
        wanted_names = {'35703': 'xenobiotic', '24432': 'biological role', '38028': 'cyanuric acid',
                        'role': 'role', '75767': 'animal metabolite',
                        '75763': 'eukaryotic metabolite', '28604': 'isofucosterol',
                        '25212': 'metabolite', '52206': 'biochemical role',
                        '84735': 'algal metabolite', '76507': 'marine metabolite',
                        '76924': 'plant metabolite', '50503': 'laxative', '23888': 'drug',
                        '85146': 'carboxymethylcellulose', '52217': 'pharmaceutical',
                        '33232': 'application', '63047': 'food emulsifier',
                        '64047': 'food additive', '63046': 'emulsifier', '78295': 'food component',
                        '52211': 'physiological role', '51086': 'chemical role'}
        wanted_roles = {'38028': {'24432', 'role', '35703'},
                        '28604': {'52206', '84735', '75767', 'role', '25212', '76507', '75763',
                                  '24432', '76924'},
                        '85146': {'52217', '50503', '63047', '51086', 'role', '63046', '23888',
                                  '33232', '78295', '64047', '52211', '24432'}}
        self.assertEqual(output, '3/3 chebi id with roles associated.')
        self.assertDictEqual(all_roles, wanted_roles)
        self.assertDictEqual(names, wanted_names)
        self.assertTrue(dicts_with_sorted_lists_equal(d_roles_ontology, wanted_ontology))

    @test_for(extract_classes)
    def test_extract_classes_chebi(self):
        ch_classes, d_classes_ontology, names = extract_classes(CHEBI, CH_LST, ROOTS[CHEBI], None,
                                                                CH_URL)
        wanted_ontology = {'35703': ['24432'], '24432': ['role'], '38028': ['35703'],
                           '75767': ['75763'], '75763': ['25212'], '25212': ['52206'],
                           '52206': ['24432'], '84735': ['75763'], '76507': ['25212'],
                           '76924': ['75763'], '28604': ['75767', '76507', '84735', '76924'],
                           '50503': ['23888'], '23888': ['52217'], '52217': ['33232'],
                           '33232': ['role'], '63047': ['63046', '64047'],
                           '64047': ['33232', '78295'], '78295': ['52211'], '52211': ['24432'],
                           '63046': ['51086'], '51086': ['role'], '85146': ['50503', '63047']}
        wanted_names = {'35703': 'xenobiotic', '24432': 'biological role', '38028': 'cyanuric acid',
                        'role': 'role', '75767': 'animal metabolite',
                        '75763': 'eukaryotic metabolite', '28604': 'isofucosterol',
                        '25212': 'metabolite', '52206': 'biochemical role',
                        '84735': 'algal metabolite', '76507': 'marine metabolite',
                        '76924': 'plant metabolite', '50503': 'laxative', '23888': 'drug',
                        '85146': 'carboxymethylcellulose', '52217': 'pharmaceutical',
                        '33232': 'application', '63047': 'food emulsifier',
                        '64047': 'food additive', '63046': 'emulsifier', '78295': 'food component',
                        '52211': 'physiological role', '51086': 'chemical role'}
        wanted_roles = {'38028': {'24432', 'role', '35703'},
                        '28604': {'52206', '84735', '75767', 'role', '25212', '76507', '75763',
                                  '24432', '76924'},
                        '85146': {'52217', '50503', '63047', '51086', 'role', '63046', '23888',
                                  '33232', '78295', '64047', '52211', '24432'}}
        self.assertEqual(ch_classes, wanted_roles)
        self.assertDictEqual(names, wanted_names)
        self.assertTrue(dicts_with_sorted_lists_equal(d_classes_ontology, wanted_ontology))

    @test_for(ontosunburst)
    def test_ontosunburst_ch1(self):
        fig = ontosunburst(interest_set=CH_LST, ontology=CHEBI, root='00',
                           abundances=None, reference_set=REF_CH, ref_abundances=None,
                           analysis='topology', output='test_ch1', write_output=False,
                           class_ontology=None, labels=None, endpoint_url=None,
                           ref_base=True, show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test_ch1.json')
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))
