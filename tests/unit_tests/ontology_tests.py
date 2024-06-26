import unittest
from unittest.mock import patch
import io
import sys
from ontosunburst.ontology import *

"""
Tests manually good file creation.
No automatic tests integrated.
"""


# ==================================================================================================
# GLOBAL
# ==================================================================================================

# GENERAL DICT ONTO (METACYC, KEGG)
# --------------------------------------------------------------------------------------------------

MET_LST = ['a', 'b', 'c']
MET_REF = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
MET_LAB = [1, 2, 3]
MET_RAB = [1, 2, 3, 4, 5, 6, 7, 8]
MC_ONTO = {'a': ['ab'], 'b': ['ab'], 'c': ['cde', 'cf'], 'd': ['cde'], 'e': ['cde', 'eg'],
           'f': ['cf'], 'g': ['gh', 'eg'], 'h': ['gh'],
           'ab': [ROOTS[METACYC]], 'cde': ['cdecf', 'cdeeg'], 'cf': ['cdecf'],
           'eg': [ROOTS[METACYC], 'cdeeg'], 'gh': [ROOTS[METACYC]],
           'cdecf': [ROOTS[METACYC]], 'cdeeg': ['cdeeg+'], 'cdeeg+': [ROOTS[METACYC]]}

# EC
# --------------------------------------------------------------------------------------------------
EC_LST = ['1.4.5.6', '1.4.6.7', '2.1.2.3', '1.5.3', '1.6.9.-', '1.-.-.-', '1.4.-.-']
EC_ONTO = {'1.4.5.-': ['1.4.-.-'], '1.4.6.-': ['1.4.-.-'], '2.1.2.-': ['2.1.-.-'],
           '1.5.3.-': ['1.5.-.-'], '1.6.9.-': ['1.6.-.-'],
           '1.4.-.-': ['1.-.-.-'], '2.1.-.-': ['2.-.-.-'], '1.5.-.-': ['1.-.-.-'],
           '1.6.-.-': ['1.-.-.-'], '1.-.-.-': [ROOTS[EC]], '2.-.-.-': [ROOTS[EC]]}

# CHEBI
# --------------------------------------------------------------------------------------------------
CH_URL = 'http://localhost:3030/chebi/'
CH_LST = ['38028', '28604', '85146']
REF_CH = ['38028', '28604', '85146',
          '23066', '27803', '37565',
          '58215', '79983', '42639']

# GO
# --------------------------------------------------------------------------------------------------
GO_LST = ['GO:0043227', 'GO:0043229', 'GO:0043231', 'GO:0044422']
GO_URL = 'http://localhost:3030/go/'


# ==================================================================================================
# FUNCTIONS UTILS
# ==================================================================================================
def dicts_with_sorted_lists_equal(dict1, dict2):
    if dict1.keys() != dict2.keys():
        return False
    for key in dict1:
        if sorted(dict1[key]) != sorted(dict2[key]):
            return False
    return True


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

# TEST CLASSES EXTRACTION
# --------------------------------------------------------------------------------------------------
class TestClassesExtraction(unittest.TestCase):
    @patch('sys.stdout', new_callable=lambda: DualWriter(sys.stdout))
    def test_extract_met_classes_1(self, mock_stdout):
        d_obj = extract_met_classes(MET_LST, MC_ONTO)
        output = mock_stdout.getvalue().strip()
        self.assertEqual(output, '3 metabolic objects to classify\n'
                                 '3/3 metabolic objects classified')
        self.assertEqual(d_obj, {'a': ['ab'], 'b': ['ab'], 'c': ['cde', 'cf']})

    @patch('sys.stdout', new_callable=lambda: DualWriter(sys.stdout))
    def test_extract_met_classes_2(self, mock_stdout):
        d_obj = extract_met_classes(MET_LST + ['x'], MC_ONTO)
        output = mock_stdout.getvalue().strip()
        self.assertEqual(output, '4 metabolic objects to classify\n'
                                 'x not classified.\n'
                                 '3/4 metabolic objects classified')
        self.assertEqual(d_obj, {'a': ['ab'], 'b': ['ab'], 'c': ['cde', 'cf']})

    @patch('sys.stdout', new_callable=lambda: DualWriter(sys.stdout))
    def test_extract_ec_classes_1(self, mock_stdout):
        d_obj, d_onto = extract_ec_classes(EC_LST, EC_ONTO)
        output = mock_stdout.getvalue().strip()
        wanted_d_obj = {'1.4.5.6': ['1.4.5.-'], '1.4.6.7': ['1.4.6.-'], '2.1.2.3': ['2.1.2.-'],
                        '1.5.3': ['1.5.-.-'], '1.6.9.-': ['1.6.-.-'], '1.-.-.-': ['Enzyme'],
                        '1.4.-.-': ['1.-.-.-']}
        self.assertEqual(output, '7 EC numbers to classify\n'
                                 '7/7 EC numbers classified')
        self.assertEqual(d_obj, wanted_d_obj)
        self.assertEqual(d_onto, {**EC_ONTO, **wanted_d_obj})

    @patch('sys.stdout', new_callable=lambda: DualWriter(sys.stdout))
    def test_extract_ec_classes_2(self, mock_stdout):
        d_obj, d_onto = extract_ec_classes(EC_LST + ['3.5.6.9', 'ecID'], EC_ONTO)
        output = mock_stdout.getvalue().strip()
        wanted_d_obj = {'1.4.5.6': ['1.4.5.-'], '1.4.6.7': ['1.4.6.-'], '2.1.2.3': ['2.1.2.-'],
                        '1.5.3': ['1.5.-.-'], '1.6.9.-': ['1.6.-.-'], '1.-.-.-': ['Enzyme'],
                        '1.4.-.-': ['1.-.-.-']}
        self.assertEqual(output, '9 EC numbers to classify\n'
                                 '3.5.6.9 not classified\n'
                                 'ecID not classified\n'
                                 '7/9 EC numbers classified')
        self.assertEqual(d_obj, wanted_d_obj)
        self.assertEqual(d_onto, {**EC_ONTO, **wanted_d_obj})

    def test_get_all_classes(self):
        leaf_classes = {'a': ['ab'], 'b': ['ab'], 'c': ['cde', 'cf']}
        all_classes_met = get_all_classes(leaf_classes, MC_ONTO, ROOTS[METACYC])
        wanted_all_classes = {'a': {'FRAMES', 'ab'}, 'b': {'FRAMES', 'ab'},
                              'c': {'cdeeg+', 'cde', 'cdeeg', 'FRAMES', 'cdecf', 'cf'}}
        self.assertEqual(all_classes_met, wanted_all_classes)

    def test_extract_classes(self):
        pass


# TEST CHEBI + GO : NEED DEDICATED SPARQL SERVER TO EXECUTE
# --------------------------------------------------------------------------------------------------
# ChEBI
class TestChEBIClassesExtraction(unittest.TestCase):
    @patch('sys.stdout', new_callable=lambda: DualWriter(sys.stdout))
    def test_extract_chebi_roles(self, mock_stdout):
        all_roles, d_roles_ontology = extract_chebi_roles(CH_LST, CH_URL)
        output = mock_stdout.getvalue().strip()
        wanted_ontology = {'xenobiotic': ['biological role'], 'biological role': ['role'],
                           '38028': ['xenobiotic'], 'algal metabolite': ['eukaryotic metabolite'],
                           'eukaryotic metabolite': ['metabolite'],
                           'metabolite': ['biochemical role'],
                           'biochemical role': ['biological role'],
                           'marine metabolite': ['metabolite'],
                           'plant metabolite': ['eukaryotic metabolite'],
                           'animal metabolite': ['eukaryotic metabolite'],
                           '28604': ['animal metabolite', 'plant metabolite', 'algal metabolite',
                                     'marine metabolite'],
                           'food emulsifier': ['food additive', 'emulsifier'],
                           'food additive': ['food component', 'application'],
                           'application': ['role'], 'food component': ['physiological role'],
                           'physiological role': ['biological role'],
                           'emulsifier': ['chemical role'], 'chemical role': ['role'],
                           'laxative': ['drug'], 'drug': ['pharmaceutical'],
                           'pharmaceutical': ['application'],
                           '85146': ['food emulsifier', 'laxative']}
        wanted_roles = {'38028': {'role', 'biological role', 'xenobiotic'},
                        '28604': {'algal metabolite', 'biochemical role', 'biological role',
                                  'eukaryotic metabolite', 'marine metabolite', 'plant metabolite',
                                  'metabolite', 'animal metabolite', 'role'},
                        '85146': {'emulsifier', 'food additive', 'pharmaceutical', 'application',
                                  'biological role', 'food component', 'chemical role', 'role',
                                  'physiological role', 'food emulsifier', 'drug', 'laxative'}}

        self.assertEqual(output, '3/3 chebi id with roles associated.')
        self.assertDictEqual(all_roles, wanted_roles)
        self.assertTrue(dicts_with_sorted_lists_equal(d_roles_ontology, wanted_ontology))


# GO
class TestGOClassesExtraction(unittest.TestCase):
    @patch('sys.stdout', new_callable=lambda: DualWriter(sys.stdout))
    def test_extract_go_classes(self, mock_stdout):
        all_classes, d_classes_ontology = extract_go_classes(GO_LST, GO_URL)
        output = mock_stdout.getvalue().strip()
        wanted_ontology = {'membrane-bounded organelle': ['organelle'],
                           'organelle': ['cellular anatomical entity'],
                           'cellular_component': ['GO'],
                           'cellular anatomical entity': ['cellular_component'],
                           'go:0043227': ['membrane-bounded organelle'],
                           'intracellular organelle': ['organelle'],
                           'go:0043229': ['intracellular organelle'],
                           'intracellular membrane-bounded organelle':
                               ['membrane-bounded organelle', 'intracellular organelle'],
                           'go:0043231': ['intracellular membrane-bounded organelle']}
        wanted_classes = {'go:0043227': {'organelle', 'cellular anatomical entity',
                                         'membrane-bounded organelle', 'cellular_component', 'GO'},
                          'go:0043229': {'organelle', 'cellular anatomical entity',
                                         'intracellular organelle', 'cellular_component', 'GO'},
                          'go:0043231': {'organelle', 'cellular anatomical entity',
                                         'intracellular organelle', 'membrane-bounded organelle',
                                         'cellular_component',
                                         'intracellular membrane-bounded organelle', 'GO'}}

        self.assertEqual(output, 'No GO class found for : go:0044422')
        self.assertDictEqual(all_classes, wanted_classes)
        self.assertTrue(dicts_with_sorted_lists_equal(d_classes_ontology, wanted_ontology))