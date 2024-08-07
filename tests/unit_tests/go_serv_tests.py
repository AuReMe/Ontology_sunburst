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

GO_LST = ['go:0043227', 'go:0043229', 'go:0043231', 'go:0044422']
GO_URL = 'http://localhost:3030/go/'


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
class TestGOClassesExtraction(unittest.TestCase):
    @test_for(extract_go_classes)
    @patch('sys.stdout', new_callable=lambda: DualWriter(sys.stdout))
    def test_extract_go_classes(self, mock_stdout):
        all_classes, d_classes_ontology, names = extract_go_classes(GO_LST, GO_URL)
        output = mock_stdout.getvalue().strip()
        wanted_ontology = {'go:0043227': ['go:0043226'],
                           'go:0043226': ['go:0110165'],
                           'go:0005575': ['GO'],
                           'go:0110165': ['go:0005575'],
                           'go:0043229': ['go:0043226'],
                           'go:0043231': ['go:0043229', 'go:0043227']}
        wanted_names = {'go:0043227': 'membrane-bounded organelle',
                        'go:0043226': 'organelle',
                        'go:0110165': 'cellular anatomical entity',
                        'go:0005575': 'cellular_component',
                        'go:0043229': 'intracellular organelle',
                        'go:0043231': 'intracellular membrane-bounded organelle'}
        wanted_classes = {'go:0043227': {'go:0043226', 'go:0110165', 'GO', 'go:0005575'},
                          'go:0043229': {'go:0043226', 'go:0110165', 'GO', 'go:0005575'},
                          'go:0043231': {'GO', 'go:0043227', 'go:0043226', 'go:0110165',
                                         'go:0005575', 'go:0043229'}}

        self.assertEqual(output, 'No GO class found for : go:0044422')
        self.assertDictEqual(all_classes, wanted_classes)
        self.assertDictEqual(names, wanted_names)
        self.assertTrue(dicts_with_sorted_lists_equal(d_classes_ontology, wanted_ontology))

    @test_for(extract_classes)
    @patch('sys.stdout', new_callable=lambda: DualWriter(sys.stdout))
    def test_extract_classes_go(self, mock_stdout):
        go_classes, d_classes_ontology, names = extract_classes(GO, GO_LST, ROOTS[GO], None, GO_URL)
        output = mock_stdout.getvalue().strip()
        wanted_ontology = {'go:0043227': ['go:0043226'],
                           'go:0043226': ['go:0110165'],
                           'go:0005575': ['GO'],
                           'go:0110165': ['go:0005575'],
                           'go:0043229': ['go:0043226'],
                           'go:0043231': ['go:0043229', 'go:0043227']}
        wanted_names = {'go:0043227': 'membrane-bounded organelle',
                        'go:0043226': 'organelle',
                        'go:0110165': 'cellular anatomical entity',
                        'go:0005575': 'cellular_component',
                        'go:0043229': 'intracellular organelle',
                        'go:0043231': 'intracellular membrane-bounded organelle'}
        wanted_classes = {'go:0043227': {'go:0043226', 'go:0110165', 'GO', 'go:0005575'},
                          'go:0043229': {'go:0043226', 'go:0110165', 'GO', 'go:0005575'},
                          'go:0043231': {'GO', 'go:0043227', 'go:0043226', 'go:0110165',
                                         'go:0005575', 'go:0043229'}}

        self.assertEqual(output, 'No GO class found for : go:0044422')
        self.assertDictEqual(go_classes, wanted_classes)
        self.assertDictEqual(names, wanted_names)
        self.assertTrue(dicts_with_sorted_lists_equal(d_classes_ontology, wanted_ontology))

    @test_for(ontosunburst)
    def test_ontosunburst_go1(self):
        fig = ontosunburst(interest_set=GO_LST, ontology=GO, root='00',
                           abundances=None, reference_set=None, ref_abundances=None,
                           analysis='topology', output='test_go1', write_output=False,
                           class_ontology=None, labels=None, endpoint_url=None, root_cut='uncut',
                           ref_base=False, show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test_go1.json')
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))
