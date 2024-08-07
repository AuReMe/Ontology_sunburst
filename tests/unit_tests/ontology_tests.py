import unittest
from unittest.mock import patch
import io
import sys
from functools import wraps
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
KG_ONTO = {'a': ['ab'], 'b': ['ab'], 'c': ['cde', 'cf'], 'd': ['cde'], 'e': ['cde', 'eg'],
           'f': ['cf'], 'g': ['gh', 'eg'], 'h': ['gh'],
           'ab': [ROOTS[KEGG]], 'cde': ['cdecf', 'cdeeg'], 'cf': ['cdecf'],
           'eg': [ROOTS[KEGG], 'cdeeg'], 'gh': [ROOTS[KEGG]],
           'cdecf': [ROOTS[KEGG]], 'cdeeg': ['cdeeg+'], 'cdeeg+': [ROOTS[KEGG]]}

# EC
# --------------------------------------------------------------------------------------------------
EC_LST = ['1.4.5.6', '1.4.6.7', '2.1.2.3', '1.5.3', '1.6.9.-', '1.-.-.-', '1.4.-.-']
EC_ONTO = {'1.4.5.-': ['1.4.-.-'], '1.4.6.-': ['1.4.-.-'], '2.1.2.-': ['2.1.-.-'],
           '1.5.3.-': ['1.5.-.-'], '1.6.9.-': ['1.6.-.-'],
           '1.4.-.-': ['1.-.-.-'], '2.1.-.-': ['2.-.-.-'], '1.5.-.-': ['1.-.-.-'],
           '1.6.-.-': ['1.-.-.-'], '1.-.-.-': [ROOTS[EC]], '2.-.-.-': [ROOTS[EC]]}
EC_ONTO_FULL = {'1.4.5.-': ['1.4.-.-'], '1.4.6.-': ['1.4.-.-'], '2.1.2.-': ['2.1.-.-'],
                '1.5.3.-': ['1.5.-.-'], '1.6.9.-': ['1.6.-.-'], '1.4.-.-': ['1.-.-.-'],
                '2.1.-.-': ['2.-.-.-'], '1.5.-.-': ['1.-.-.-'], '1.6.-.-': ['1.-.-.-'],
                '1.-.-.-': [ROOTS[EC]], '2.-.-.-': [ROOTS[EC]], '1.4.5.6': ['1.4.5.-'],
                '1.4.6.7': ['1.4.6.-'], '2.1.2.3': ['2.1.2.-'], '1.5.3': ['1.5.-.-']}


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

# TEST CLASSES EXTRACTION
# --------------------------------------------------------------------------------------------------
class TestClassesExtraction(unittest.TestCase):
    @test_for(extract_met_classes)
    @patch('sys.stdout', new_callable=lambda: DualWriter(sys.stdout))
    def test_extract_met_classes_input_ok(self, mock_stdout):
        d_obj = extract_met_classes(MET_LST, MC_ONTO)
        output = mock_stdout.getvalue().strip()
        self.assertEqual(output, '3 metabolic objects to classify\n'
                                 '3/3 metabolic objects classified')
        self.assertEqual(d_obj, {'a': ['ab'], 'b': ['ab'], 'c': ['cde', 'cf']})

    @test_for(extract_met_classes)
    @patch('sys.stdout', new_callable=lambda: DualWriter(sys.stdout))
    def test_extract_met_classes_input_errors(self, mock_stdout):
        d_obj = extract_met_classes(MET_LST + ['x'], MC_ONTO)
        output = mock_stdout.getvalue().strip()
        self.assertEqual(output, '4 metabolic objects to classify\n'
                                 'x not classified.\n'
                                 '3/4 metabolic objects classified')
        self.assertEqual(d_obj, {'a': ['ab'], 'b': ['ab'], 'c': ['cde', 'cf']})

    @test_for(extract_ec_classes)
    @patch('sys.stdout', new_callable=lambda: DualWriter(sys.stdout))
    def test_extract_ec_classes_input_ok(self, mock_stdout):
        d_obj, d_onto = extract_ec_classes(EC_LST, EC_ONTO)
        output = mock_stdout.getvalue().strip()
        wanted_d_obj = {'1.4.5.6': ['1.4.5.-'], '1.4.6.7': ['1.4.6.-'], '2.1.2.3': ['2.1.2.-'],
                        '1.5.3': ['1.5.-.-'], '1.6.9.-': ['1.6.-.-'], '1.-.-.-': ['Enzyme'],
                        '1.4.-.-': ['1.-.-.-']}
        self.assertEqual(output, '7 EC numbers to classify\n'
                                 '7/7 EC numbers classified')
        self.assertEqual(d_obj, wanted_d_obj)
        self.assertEqual(d_onto, EC_ONTO_FULL)

    @test_for(extract_ec_classes)
    @patch('sys.stdout', new_callable=lambda: DualWriter(sys.stdout))
    def test_extract_ec_classes_input_errors(self, mock_stdout):
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
        self.assertEqual(d_onto, EC_ONTO_FULL)

    @test_for(get_parents)
    def test_get_parents_linear_path(self):
        # Simple linear direction
        parents = get_parents('a', {'ab'}, MC_ONTO, ROOTS[METACYC])
        self.assertEqual(parents, {'FRAMES', 'ab'})

    @test_for(get_parents)
    def test_get_parents_complex_path(self):
        # With multiple parents having multiple parents and different size of path until root
        parents = get_parents('c', {'cde', 'cf'}, MC_ONTO, ROOTS[METACYC])
        self.assertEqual(parents, {'cdeeg+', 'FRAMES', 'cf', 'cde', 'cdecf', 'cdeeg'})

    @test_for(get_parents)
    def test_get_parents_ec(self):
        # With EC (simple path)
        parents = get_parents('2.1.2.3', {'2.1.2.-'}, EC_ONTO_FULL, ROOTS[EC])
        self.assertEqual(parents, {'Enzyme', '2.1.-.-', '2.-.-.-', '2.1.2.-'})

    @test_for(get_parents)
    def test_get_parents_direct_root_child(self):
        # With direct child of root item
        parents = get_parents('1.-.-.-', {'Enzyme'}, EC_ONTO_FULL, ROOTS[EC])
        self.assertEqual(parents, {'Enzyme'})

    @test_for(get_all_classes)
    def test_get_all_classes(self):
        leaf_classes = {'a': ['ab'], 'b': ['ab'], 'c': ['cde', 'cf']}
        all_classes_met = get_all_classes(leaf_classes, MC_ONTO, ROOTS[METACYC])
        wanted_all_classes = {'a': {'FRAMES', 'ab'}, 'b': {'FRAMES', 'ab'},
                              'c': {'cdeeg+', 'cde', 'cdeeg', 'FRAMES', 'cdecf', 'cf'}}
        self.assertEqual(all_classes_met, wanted_all_classes)

    @test_for(get_all_classes)
    def test_get_all_classes_ec(self):
        ec_leaf_classes = {'1.4.5.6': ['1.4.5.-'], '1.4.6.7': ['1.4.6.-'], '2.1.2.3': ['2.1.2.-'],
                           '1.5.3': ['1.5.-.-'], '1.6.9.-': ['1.6.-.-'], '1.-.-.-': ['Enzyme'],
                           '1.4.-.-': ['1.-.-.-']}
        wanted_all_classes = {'1.4.5.6': {'1.4.5.-', '1.4.-.-', 'Enzyme', '1.-.-.-'},
                              '1.4.6.7': {'Enzyme', '1.4.-.-', '1.4.6.-', '1.-.-.-'},
                              '2.1.2.3': {'2.-.-.-', '2.1.-.-', '2.1.2.-', 'Enzyme'},
                              '1.5.3': {'1.5.-.-', 'Enzyme', '1.-.-.-'},
                              '1.6.9.-': {'1.6.-.-', 'Enzyme', '1.-.-.-'},
                              '1.-.-.-': {'Enzyme'}, '1.4.-.-': {'Enzyme', '1.-.-.-'}}
        all_classes_ec = get_all_classes(ec_leaf_classes, EC_ONTO, ROOTS[EC])
        self.assertEqual(all_classes_ec, wanted_all_classes)

    @test_for(extract_classes)
    def test_extract_classes_metacyc(self):
        mc_classes, d_classes_ontology, names = extract_classes(METACYC, MET_LST, ROOTS[METACYC],
                                                                MC_ONTO, None)
        wanted_mc_classes = {'a': {'FRAMES', 'ab'}, 'b': {'FRAMES', 'ab'},
                             'c': {'cdeeg+', 'cde', 'cdeeg', 'FRAMES', 'cdecf', 'cf'}}
        self.assertEqual(mc_classes, wanted_mc_classes)
        self.assertTrue(dicts_with_sorted_lists_equal(d_classes_ontology, MC_ONTO))

    @test_for(extract_classes)
    def test_extract_classes_kegg(self):
        kg_classes, d_classes_ontology, names = extract_classes(KEGG, MET_LST, ROOTS[KEGG], KG_ONTO,
                                                                None)
        wanted_kg_classes = {'a': {'ab', 'kegg'}, 'b': {'ab', 'kegg'},
                             'c': {'cde', 'cdeeg+', 'kegg', 'cdecf', 'cdeeg', 'cf'}}
        self.assertEqual(kg_classes, wanted_kg_classes)
        self.assertTrue(dicts_with_sorted_lists_equal(d_classes_ontology, KG_ONTO))

    @test_for(extract_classes)
    def test_extract_classes_ec(self):
        ec_classes, d_classes_ontology, names = extract_classes(EC, EC_LST, ROOTS[EC], EC_ONTO,
                                                                None)
        wanted_ec_classes = {'1.4.5.6': {'1.4.5.-', '1.4.-.-', 'Enzyme', '1.-.-.-'},
                             '1.4.6.7': {'Enzyme', '1.4.-.-', '1.4.6.-', '1.-.-.-'},
                             '2.1.2.3': {'2.-.-.-', '2.1.-.-', '2.1.2.-', 'Enzyme'},
                             '1.5.3': {'1.5.-.-', 'Enzyme', '1.-.-.-'},
                             '1.6.9.-': {'1.6.-.-', 'Enzyme', '1.-.-.-'},
                             '1.-.-.-': {'Enzyme'}, '1.4.-.-': {'Enzyme', '1.-.-.-'}}
        self.assertEqual(ec_classes, wanted_ec_classes)
        self.assertTrue(dicts_with_sorted_lists_equal(d_classes_ontology, EC_ONTO_FULL))


# TEST ABUNDANCES
# --------------------------------------------------------------------------------------------------
class TestAbundances(unittest.TestCase):
    @test_for(get_abundance_dict)
    def test_get_abundance_dict_abundances_not_ref(self):
        abundance_dict = get_abundance_dict(abundances=MET_LAB, metabolic_objects=MET_LST,
                                            ref=False)
        self.assertEqual(abundance_dict, {'a': 1, 'b': 2, 'c': 3})

    @test_for(get_abundance_dict)
    def test_get_abundance_dict_no_abundances_not_ref(self):
        abundance_dict = get_abundance_dict(abundances=None, metabolic_objects=MET_LST,
                                            ref=False)
        self.assertEqual(abundance_dict, {'a': 1, 'b': 1, 'c': 1})

    @test_for(get_abundance_dict)
    def test_get_abundance_dict_abundances_ref(self):
        abundance_dict = get_abundance_dict(abundances=MET_RAB, metabolic_objects=MET_REF,
                                            ref=True)
        self.assertEqual(abundance_dict, {'a': 1, 'b': 2, 'c': 3, 'd': 4,
                                          'e': 5, 'f': 6, 'g': 7, 'h': 8})

    @test_for(get_abundance_dict)
    def test_get_abundance_dict_no_abundances_ref(self):
        abundance_dict = get_abundance_dict(abundances=None, metabolic_objects=MET_REF,
                                            ref=True)
        self.assertEqual(abundance_dict, {'a': 1, 'b': 1, 'c': 1, 'd': 1,
                                          'e': 1, 'f': 1, 'g': 1, 'h': 1})

    @test_for(get_abundance_dict)
    def test_get_abundance_dict_errors_not_ref(self):
        with self.assertRaises(AttributeError) as e:
            get_abundance_dict(abundances=MET_LAB + [4], metabolic_objects=MET_LST, ref=False)
        self.assertEqual(str(e.exception), 'Length of "metabolic_objects" parameter must be '
                                           'equal to "abundances" parameter length : 3 != 4')

    @test_for(get_abundance_dict)
    def test_get_abundance_dict_errors_ref(self):
        with self.assertRaises(AttributeError) as e:
            get_abundance_dict(abundances=MET_RAB[:-1], metabolic_objects=MET_REF, ref=True)
        self.assertEqual(str(e.exception), 'Length of "reference_set" parameter must be '
                                           'equal to "ref_abundances" parameter length : 8 != 7')

    @test_for(get_classes_abundance)
    def test_get_classes_abundance_leaves(self):
        all_classes = {'a': {'FRAMES', 'ab'}, 'b': {'FRAMES', 'ab'},
                       'c': {'cdecf', 'cdeeg+', 'FRAMES', 'cde', 'cdeeg', 'cf'},
                       'd': {'cdecf', 'cdeeg+', 'FRAMES', 'cde', 'cdeeg'},
                       'e': {'cdeeg+', 'FRAMES', 'cde', 'cdecf', 'eg', 'cdeeg'},
                       'f': {'cdecf', 'FRAMES', 'cf'},
                       'g': {'cdeeg', 'cdeeg+', 'FRAMES', 'eg', 'gh'}, 'h': {'FRAMES', 'gh'}}
        abundances_dict = {'a': 1, 'b': 2, 'c': 3, 'd': 4, 'e': 5, 'f': 6, 'g': 7, 'h': 8}
        classes_abundances = get_classes_abundance(all_classes, abundances_dict, show_leaves=True)
        wanted_abundances = {'FRAMES': 36, 'cdeeg+': 19, 'cdeeg': 19, 'cdecf': 18, 'gh': 15,
                             'eg': 12, 'cde': 12, 'cf': 9, 'h': 8, 'g': 7, 'f': 6, 'e': 5,
                             'd': 4, 'c': 3, 'ab': 3, 'b': 2, 'a': 1}
        self.assertEqual(classes_abundances, wanted_abundances)

    @test_for(get_classes_abundance)
    def test_get_classes_abundance_leaves_sub(self):
        all_classes = {'a': {'FRAMES', 'ab'}, 'b': {'FRAMES', 'ab'},
                       'c': {'cdeeg+', 'cde', 'cdeeg', 'FRAMES', 'cdecf', 'cf'}}
        abundances_dict = {'a': 1, 'b': 2, 'c': 3}
        classes_abundances = get_classes_abundance(all_classes, abundances_dict, show_leaves=True)
        wanted_abundances = {'FRAMES': 6, 'cde': 3, 'cf': 3, 'cdecf': 3, 'cdeeg+': 3, 'cdeeg': 3,
                             'c': 3, 'ab': 3, 'b': 2, 'a': 1}
        self.assertEqual(classes_abundances, wanted_abundances)

    @test_for(get_classes_abundance)
    def test_get_classes_abundance_no_leaves(self):
        all_classes = {'a': {'FRAMES', 'ab'}, 'b': {'FRAMES', 'ab'},
                       'c': {'cdecf', 'cdeeg+', 'FRAMES', 'cde', 'cdeeg', 'cf'},
                       'd': {'cdecf', 'cdeeg+', 'FRAMES', 'cde', 'cdeeg'},
                       'e': {'cdeeg+', 'FRAMES', 'cde', 'cdecf', 'eg', 'cdeeg'},
                       'f': {'cdecf', 'FRAMES', 'cf'},
                       'g': {'cdeeg', 'cdeeg+', 'FRAMES', 'eg', 'gh'}, 'h': {'FRAMES', 'gh'}}
        abundances_dict = {'a': 1, 'b': 2, 'c': 3, 'd': 4, 'e': 5, 'f': 6, 'g': 7, 'h': 8}
        classes_abundances = get_classes_abundance(all_classes, abundances_dict, show_leaves=False)
        wanted_abundances = {'FRAMES': 36, 'cdeeg+': 19, 'cdeeg': 19, 'cdecf': 18, 'gh': 15,
                             'eg': 12, 'cde': 12, 'cf': 9, 'ab': 3}
        self.assertEqual(classes_abundances, wanted_abundances)

    @test_for(get_classes_abundance)
    def test_get_classes_abundance_no_leaves_sub(self):
        all_classes = {'a': {'FRAMES', 'ab'}, 'b': {'FRAMES', 'ab'},
                       'c': {'cdeeg+', 'cde', 'cdeeg', 'FRAMES', 'cdecf', 'cf'}}
        abundances_dict = {'a': 1, 'b': 2, 'c': 3}
        classes_abundances = get_classes_abundance(all_classes, abundances_dict, show_leaves=False)
        wanted_abundances = {'FRAMES': 6, 'cde': 3, 'cf': 3, 'cdecf': 3, 'cdeeg+': 3, 'cdeeg': 3,
                             'ab': 3}
        self.assertEqual(classes_abundances, wanted_abundances)

    @test_for(get_classes_abundance)
    def test_get_classes_abundance_different_level_abundances(self):
        all_classes = {'c': {'cdecf', 'cdeeg+', 'FRAMES', 'cde', 'cdeeg', 'cf'},
                       'd': {'cdecf', 'cdeeg+', 'FRAMES', 'cde', 'cdeeg'},
                       'e': {'cdeeg+', 'FRAMES', 'cde', 'cdecf', 'eg', 'cdeeg'},
                       'f': {'cdecf', 'FRAMES', 'cf'}, 'cf': {'cdecf', 'FRAMES'}}
        abundances_dict = {'c': 3, 'd': 4, 'e': 5, 'f': 2, 'cf': 2}
        classes_abundances = get_classes_abundance(all_classes, abundances_dict, show_leaves=True)
        wanted_abundances = {'FRAMES': 16, 'cdecf': 16, 'cde': 12, 'cdeeg+': 12, 'cdeeg': 12,
                             'cf': 7, 'eg': 5, 'e': 5, 'd': 4, 'c': 3, 'f': 2}
        self.assertEqual(classes_abundances, wanted_abundances)


# TEST UTILS
# --------------------------------------------------------------------------------------------------
class TestUtils(unittest.TestCase):
    @test_for(reduce_d_ontology)
    def test_reduce_d_ontology(self):
        classes_abundance = {'FRAMES': 6, 'cde': 3, 'cf': 3, 'cdecf': 3, 'cdeeg+': 3, 'cdeeg': 3,
                             'ab': 3}
        d_ontology_reduced = reduce_d_ontology(MC_ONTO, classes_abundance)
        wanted_d_ontology_reduced = {'ab': ['FRAMES'], 'cde': ['cdecf', 'cdeeg'],
                                     'cf': ['cdecf'], 'cdecf': ['FRAMES'], 'cdeeg': ['cdeeg+'],
                                     'cdeeg+': ['FRAMES']}
        self.assertEqual(d_ontology_reduced, wanted_d_ontology_reduced)


