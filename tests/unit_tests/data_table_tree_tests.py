import unittest
import io

from functools import wraps
from ontosunburst.data_table_tree import *
from ontosunburst.ontology import *

"""
Tests manually good file creation.
No automatic tests integrated.
"""

# ==================================================================================================
# GLOBAL
# ==================================================================================================

# --------------------------------------------------------------------------------------------------

CT_ONTO = {'a': ['ab'], 'b': ['ab'], 'c': ['cde', 'cf'], 'd': ['cde'], 'e': ['cde', 'eg'],
           'f': ['cf'], 'g': ['gh', 'eg'], 'h': ['gh'],
           'ab': [ROOTS[METACYC]], 'cde': ['cdecf', 'cdeeg'], 'cf': ['cdecf'],
           'eg': ['cdeeg', ROOTS[METACYC]], 'gh': [ROOTS[METACYC]],
           'cdecf': [ROOTS[METACYC]], 'cdeeg': ['cdeeg+'], 'cdeeg+': [ROOTS[METACYC]]}

CT_AB = {'FRAMES': 6, 'cde': 3, 'cf': 3, 'cdecf': 3, 'cdeeg+': 3, 'cdeeg': 3, 'c': 3, 'ab': 3,
         'b': 2, 'a': 1}
CT_REF_AB = {'FRAMES': 36, 'cdeeg+': 19, 'cdeeg': 19, 'cdecf': 18, 'gh': 15, 'eg': 12, 'cde': 12,
             'cf': 9, 'h': 8, 'g': 7, 'f': 6, 'e': 5, 'd': 4, 'c': 3, 'ab': 3, 'b': 2, 'a': 1}

CT_LAB = {'FRAMES': 'Root', 'cdeeg+': 'CDEEG+', 'cdeeg': 'CDEEG', 'cdecf': 'CDECF', 'gh': 'GH',
          'eg': 'EG', 'cde': 'CDE', 'cf': 'CF', 'h': 'H', 'g': 'G', 'f': 'F', 'e': 'E', 'd': 'D',
          'c': 'C', 'ab': 'AB', 'b': 'B'}

W_PROP = [1.0, 0.5, 0.5, 0.5, nan, nan, nan, 0.5, 0.5, 0.5, nan, nan, nan, nan, nan, nan, nan, nan,
          nan, nan, nan, 0.5, 0.5, 0.5, 0.5, 0.3333333333333333, 0.16666666666666666]

W_REF_PROP = [1.0, 0.5277777777777778, 0.5277777777777778, 0.5, 0.4166666666666667,
              0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.25,
              0.2222222222222222, 0.19444444444444445, 0.19444444444444445, 0.19444444444444445,
              0.16666666666666666, 0.1388888888888889, 0.1388888888888889, 0.1388888888888889,
              0.1388888888888889, 0.1111111111111111, 0.1111111111111111, 0.08333333333333333,
              0.08333333333333333, 0.08333333333333333, 0.08333333333333333, 0.05555555555555555,
              0.027777777777777776]

W_REL_PROP = {'FRAMES': 1000000, 'cdeeg+__FRAMES': 283582, 'cdeeg__cdeeg+__FRAMES': 283582,
              'cdecf__FRAMES': 268656, 'gh__FRAMES': 223880, 'eg__cdeeg__cdeeg+__FRAMES': 141791,
              'eg__FRAMES': 179104, 'cde__cdeeg__cdeeg+__FRAMES': 141791,
              'cde__cdecf__FRAMES': 153517, 'cf__cdecf__FRAMES': 115138, 'h__gh__FRAMES': 119402,
              'g__eg__cdeeg__cdeeg+__FRAMES': 82711, 'g__eg__FRAMES': 104477,
              'g__gh__FRAMES': 104477, 'f__cf__cdecf__FRAMES': 76758,
              'e__eg__cdeeg__cdeeg+__FRAMES': 59079, 'e__cde__cdeeg__cdeeg+__FRAMES': 59079,
              'e__eg__FRAMES': 74626, 'e__cde__cdecf__FRAMES': 63965,
              'd__cde__cdecf__FRAMES': 51172, 'd__cde__cdeeg__cdeeg+__FRAMES': 47263,
              'c__cde__cdecf__FRAMES': 38379, 'c__cf__cdecf__FRAMES': 38379,
              'c__cde__cdeeg__cdeeg+__FRAMES': 35447, 'ab__FRAMES': 44776, 'b__ab__FRAMES': 29850,
              'a__ab__FRAMES': 14925}


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


def data_to_lines(dico):
    lines = set()
    for i in range(len(dico[IDS])):
        line = (dico[IDS][i], dico[ONTO_ID][i], dico[PARENT][i], dico[LABEL][i], dico[COUNT][i],
                dico[REF_COUNT][i])
        if PROP in dico:
            line = line + (dico[PROP][i],)
        if REF_PROP in dico:
            line = line + (dico[REF_PROP][i],)
        if RELAT_PROP in dico:
            line = line + (dico[RELAT_PROP][i],)
        if PVAL in dico:
            line = line + (dico[PVAL][i],)
        lines.add(line)
    return lines


# ==================================================================================================
# UNIT TESTS
# ==================================================================================================

# TEST
# --------------------------------------------------------------------------------------------------

class TestGenerateDataTable(unittest.TestCase):

    @test_for(get_set2_abundance)
    def test_get_sub_abundances_exists_diff(self):
        sub_abu = get_set2_abundance(CT_AB, 'cf')
        self.assertEqual(sub_abu, 3)

    @test_for(get_set2_abundance)
    def test_get_sub_abundances_exists_equ(self):
        sub_abu = get_set2_abundance(CT_AB, 'a')
        self.assertEqual(sub_abu, 1)

    @test_for(get_set2_abundance)
    def test_get_sub_abundances_not_exists(self):
        sub_abu = get_set2_abundance(CT_AB, 'eg')
        self.assertTrue(np.isnan(sub_abu))

    @test_for(get_all_ids)
    def test_get_all_c_ids(self):
        all_ids = get_all_ids('c', 'c', CT_ONTO, ROOTS[METACYC], set())
        wanted_ids = {'c__cf__cdecf__FRAMES', 'c__cde__cdecf__FRAMES',
                      'c__cde__cdeeg__cdeeg+__FRAMES'}
        self.assertEqual(all_ids, wanted_ids)

    @test_for(get_all_ids)
    def test_get_all_e_ids(self):
        all_ids = get_all_ids('e', 'e', CT_ONTO, ROOTS[METACYC], set())
        wanted_ids = {'e__cde__cdeeg__cdeeg+__FRAMES', 'e__eg__FRAMES',
                      'e__eg__cdeeg__cdeeg+__FRAMES', 'e__cde__cdecf__FRAMES'}
        self.assertEqual(all_ids, wanted_ids)

    @test_for(get_all_ids)
    def test_get_all_eg_ids(self):
        all_ids = get_all_ids('eg', 'eg', CT_ONTO, ROOTS[METACYC], set())
        wanted_ids = {'eg__FRAMES', 'eg__cdeeg__cdeeg+__FRAMES'}
        self.assertEqual(all_ids, wanted_ids)

    @test_for(DataTable.add_value)
    def test_add_value_data(self):
        data = DataTable()
        data.add_value(m_id='bjr', onto_id='Bjr_0', label='bonjour', count=2, ref_count=8,
                       parent='salutations')
        data.add_value(m_id='slt', onto_id='sl_1', label='salut', count=0.5, ref_count=2.3,
                       parent='salutations')
        wanted_data = {IDS: ['bjr', 'slt'],
                       ONTO_ID: ['Bjr_0', 'sl_1'],
                       PARENT: ['salutations', 'salutations'],
                       LABEL: ['bonjour', 'salut'],
                       COUNT: [2, 0.5],
                       REF_COUNT: [8, 2.3],
                       PROP: [nan, nan], REF_PROP: [nan, nan], RELAT_PROP: [nan, nan],
                       PVAL: [nan, nan]}
        self.assertEqual(data.get_data_dict(), wanted_data)

    @test_for(DataTable.fill_parameters)
    def test_get_fig_parameters(self):
        data = DataTable()
        data.fill_parameters(ref_abundance=CT_REF_AB, parent_dict=CT_ONTO,
                             root_item=ROOTS[METACYC], set_abundance=CT_AB, names=None)
        lines = set(data.get_col())
        w_lines = {('cde__cdeeg__cdeeg+__FRAMES', 'cde', 'cde', 'cdeeg__cdeeg+__FRAMES', 3, 12, nan,
                    nan, nan, nan), (
                       'd__cde__cdeeg__cdeeg+__FRAMES', 'd', 'd', 'cde__cdeeg__cdeeg+__FRAMES', nan,
                       4,
                       nan, nan, nan, nan),
                   ('b__ab__FRAMES', 'b', 'b', 'ab__FRAMES', 2, 2, nan, nan, nan, nan),
                   ('g__eg__FRAMES', 'g', 'g', 'eg__FRAMES', nan, 7, nan, nan, nan, nan),
                   ('eg__FRAMES', 'eg', 'eg', 'FRAMES', nan, 12, nan, nan, nan, nan), (
                       'c__cde__cdecf__FRAMES', 'c', 'c', 'cde__cdecf__FRAMES', 3, 3, nan, nan, nan,
                       nan), (
                       'g__eg__cdeeg__cdeeg+__FRAMES', 'g', 'g', 'eg__cdeeg__cdeeg+__FRAMES', nan,
                       7,
                       nan, nan, nan, nan), (
                       'e__eg__cdeeg__cdeeg+__FRAMES', 'e', 'e', 'eg__cdeeg__cdeeg+__FRAMES', nan,
                       5,
                       nan, nan, nan, nan), (
                       'cdeeg__cdeeg+__FRAMES', 'cdeeg', 'cdeeg', 'cdeeg+__FRAMES', 3, 19, nan, nan,
                       nan, nan),
                   ('cdecf__FRAMES', 'cdecf', 'cdecf', 'FRAMES', 3, 18, nan, nan, nan, nan),
                   ('ab__FRAMES', 'ab', 'ab', 'FRAMES', 3, 3, nan, nan, nan, nan), (
                       'e__cde__cdeeg__cdeeg+__FRAMES', 'e', 'e', 'cde__cdeeg__cdeeg+__FRAMES', nan,
                       5,
                       nan, nan, nan, nan), (
                       'c__cde__cdeeg__cdeeg+__FRAMES', 'c', 'c', 'cde__cdeeg__cdeeg+__FRAMES', 3,
                       3,
                       nan, nan, nan, nan), (
                       'c__cf__cdecf__FRAMES', 'c', 'c', 'cf__cdecf__FRAMES', 3, 3, nan, nan, nan,
                       nan),
                   ('cde__cdecf__FRAMES', 'cde', 'cde', 'cdecf__FRAMES', 3, 12, nan, nan, nan, nan),
                   ('FRAMES', 'FRAMES', 'FRAMES', '', 6, 36, nan, nan, nan, nan), (
                       'f__cf__cdecf__FRAMES', 'f', 'f', 'cf__cdecf__FRAMES', nan, 6, nan, nan, nan,
                       nan), (
                       'eg__cdeeg__cdeeg+__FRAMES', 'eg', 'eg', 'cdeeg__cdeeg+__FRAMES', nan, 12,
                       nan,
                       nan, nan, nan),
                   ('e__eg__FRAMES', 'e', 'e', 'eg__FRAMES', nan, 5, nan, nan, nan, nan),
                   ('g__gh__FRAMES', 'g', 'g', 'gh__FRAMES', nan, 7, nan, nan, nan, nan),
                   ('h__gh__FRAMES', 'h', 'h', 'gh__FRAMES', nan, 8, nan, nan, nan, nan),
                   ('cdeeg+__FRAMES', 'cdeeg+', 'cdeeg+', 'FRAMES', 3, 19, nan, nan, nan, nan),
                   ('gh__FRAMES', 'gh', 'gh', 'FRAMES', nan, 15, nan, nan, nan, nan), (
                       'e__cde__cdecf__FRAMES', 'e', 'e', 'cde__cdecf__FRAMES', nan, 5, nan, nan,
                       nan,
                       nan),
                   ('cf__cdecf__FRAMES', 'cf', 'cf', 'cdecf__FRAMES', 3, 9, nan, nan, nan, nan),
                   ('a__ab__FRAMES', 'a', 'a', 'ab__FRAMES', 1, 1, nan, nan, nan, nan), (
                       'd__cde__cdecf__FRAMES', 'd', 'd', 'cde__cdecf__FRAMES', nan, 4, nan, nan,
                       nan,
                       nan)}
        self.assertEqual(lines, w_lines)

    @test_for(DataTable.fill_parameters)
    def test_get_fig_parameters_names(self):
        data = DataTable()
        data.fill_parameters(ref_abundance=CT_REF_AB, parent_dict=CT_ONTO,
                             root_item=ROOTS[METACYC], set_abundance=CT_AB, names=CT_LAB)
        lines = set(data.get_col())
        w_lines = {('ab__FRAMES', 'ab', 'AB', 'FRAMES', 3, 3, nan, nan, nan, nan), (
            'd__cde__cdecf__FRAMES', 'd', 'D', 'cde__cdecf__FRAMES', nan, 4, nan, nan, nan, nan),
                   ('b__ab__FRAMES', 'b', 'B', 'ab__FRAMES', 2, 2, nan, nan, nan, nan),
                   ('h__gh__FRAMES', 'h', 'H', 'gh__FRAMES', nan, 8, nan, nan, nan, nan), (
                       'g__eg__cdeeg__cdeeg+__FRAMES', 'g', 'G', 'eg__cdeeg__cdeeg+__FRAMES', nan,
                       7,
                       nan, nan, nan, nan),
                   ('cdeeg+__FRAMES', 'cdeeg+', 'CDEEG+', 'FRAMES', 3, 19, nan, nan, nan, nan),
                   ('FRAMES', 'FRAMES', 'FRAMES', '', 6, 36, nan, nan, nan, nan),
                   ('g__eg__FRAMES', 'g', 'G', 'eg__FRAMES', nan, 7, nan, nan, nan, nan), (
                       'cde__cdeeg__cdeeg+__FRAMES', 'cde', 'CDE', 'cdeeg__cdeeg+__FRAMES', 3, 12,
                       nan,
                       nan, nan, nan), (
                       'eg__cdeeg__cdeeg+__FRAMES', 'eg', 'EG', 'cdeeg__cdeeg+__FRAMES', nan, 12,
                       nan,
                       nan, nan, nan), (
                       'c__cde__cdeeg__cdeeg+__FRAMES', 'c', 'C', 'cde__cdeeg__cdeeg+__FRAMES', 3,
                       3,
                       nan, nan, nan, nan),
                   ('g__gh__FRAMES', 'g', 'G', 'gh__FRAMES', nan, 7, nan, nan, nan, nan),
                   ('e__eg__FRAMES', 'e', 'E', 'eg__FRAMES', nan, 5, nan, nan, nan, nan),
                   ('eg__FRAMES', 'eg', 'EG', 'FRAMES', nan, 12, nan, nan, nan, nan), (
                       'e__eg__cdeeg__cdeeg+__FRAMES', 'e', 'E', 'eg__cdeeg__cdeeg+__FRAMES', nan,
                       5,
                       nan, nan, nan, nan), (
                       'c__cf__cdecf__FRAMES', 'c', 'C', 'cf__cdecf__FRAMES', 3, 3, nan, nan, nan,
                       nan),
                   ('gh__FRAMES', 'gh', 'GH', 'FRAMES', nan, 15, nan, nan, nan, nan), (
                       'cdeeg__cdeeg+__FRAMES', 'cdeeg', 'CDEEG', 'cdeeg+__FRAMES', 3, 19, nan, nan,
                       nan, nan), (
                       'e__cde__cdeeg__cdeeg+__FRAMES', 'e', 'E', 'cde__cdeeg__cdeeg+__FRAMES', nan,
                       5,
                       nan, nan, nan, nan), (
                       'd__cde__cdeeg__cdeeg+__FRAMES', 'd', 'D', 'cde__cdeeg__cdeeg+__FRAMES', nan,
                       4,
                       nan, nan, nan, nan), (
                       'c__cde__cdecf__FRAMES', 'c', 'C', 'cde__cdecf__FRAMES', 3, 3, nan, nan, nan,
                       nan), ('a__ab__FRAMES', 'a', 'a', 'ab__FRAMES', 1, 1, nan, nan, nan, nan), (
                       'e__cde__cdecf__FRAMES', 'e', 'E', 'cde__cdecf__FRAMES', nan, 5, nan, nan,
                       nan,
                       nan),
                   ('cdecf__FRAMES', 'cdecf', 'CDECF', 'FRAMES', 3, 18, nan, nan, nan, nan),
                   ('cde__cdecf__FRAMES', 'cde', 'CDE', 'cdecf__FRAMES', 3, 12, nan, nan, nan, nan),
                   ('cf__cdecf__FRAMES', 'cf', 'CF', 'cdecf__FRAMES', 3, 9, nan, nan, nan, nan), (
                       'f__cf__cdecf__FRAMES', 'f', 'F', 'cf__cdecf__FRAMES', nan, 6, nan, nan, nan,
                       nan)}
        self.assertEqual(lines, w_lines)

    @test_for(DataTable.fill_parameters)
    def test_get_fig_parameters_no_ref_base(self):
        data = DataTable()
        data.fill_parameters(ref_abundance=CT_REF_AB, parent_dict=CT_ONTO,
                             root_item=ROOTS[METACYC], set_abundance=CT_AB, names=None,
                             ref_base=False)
        lines = set(data.get_col())
        w_lines = {('a__ab__FRAMES', 'a', 'a', 'ab__FRAMES', 1, 1, nan, nan, nan, nan),
                   ('cde__cdecf__FRAMES', 'cde', 'cde', 'cdecf__FRAMES', 3, 12, nan, nan, nan, nan),
                   ('c__cde__cdeeg__cdeeg+__FRAMES', 'c', 'c', 'cde__cdeeg__cdeeg+__FRAMES', 3, 3,
                    nan, nan, nan, nan),
                   ('c__cde__cdecf__FRAMES', 'c', 'c', 'cde__cdecf__FRAMES', 3, 3, nan, nan, nan,
                    nan),
                   ('cdeeg__cdeeg+__FRAMES', 'cdeeg', 'cdeeg', 'cdeeg+__FRAMES', 3, 19, nan, nan,
                    nan, nan),
                   ('FRAMES', 'FRAMES', 'FRAMES', '', 6, 36, nan, nan, nan, nan),
                   ('b__ab__FRAMES', 'b', 'b', 'ab__FRAMES', 2, 2, nan, nan, nan, nan),
                   ('cde__cdeeg__cdeeg+__FRAMES', 'cde', 'cde', 'cdeeg__cdeeg+__FRAMES', 3, 12, nan,
                    nan, nan, nan),
                   ('ab__FRAMES', 'ab', 'ab', 'FRAMES', 3, 3, nan, nan, nan, nan),
                   ('c__cf__cdecf__FRAMES', 'c', 'c', 'cf__cdecf__FRAMES', 3, 3, nan, nan, nan,
                    nan),
                   ('cf__cdecf__FRAMES', 'cf', 'cf', 'cdecf__FRAMES', 3, 9, nan, nan, nan, nan),
                   ('cdeeg+__FRAMES', 'cdeeg+', 'cdeeg+', 'FRAMES', 3, 19, nan, nan, nan, nan),
                   ('cdecf__FRAMES', 'cdecf', 'cdecf', 'FRAMES', 3, 18, nan, nan, nan, nan)}
        self.assertEqual(lines, w_lines)

    @test_for(DataTable.fill_parameters)
    def test_get_fig_parameters_no_ref(self):
        data = DataTable()
        data.fill_parameters(ref_abundance=CT_AB, parent_dict=CT_ONTO,
                             root_item=ROOTS[METACYC], set_abundance=CT_AB, names=None,
                             ref_base=False)
        lines = set(data.get_col())
        w_lines = {('b__ab__FRAMES', 'b', 'b', 'ab__FRAMES', 2, 2, nan, nan, nan, nan),
                   ('ab__FRAMES', 'ab', 'ab', 'FRAMES', 3, 3, nan, nan, nan, nan),
                   ('cde__cdeeg__cdeeg+__FRAMES', 'cde', 'cde', 'cdeeg__cdeeg+__FRAMES', 3, 3, nan,
                    nan, nan, nan),
                   ('c__cde__cdecf__FRAMES', 'c', 'c', 'cde__cdecf__FRAMES', 3, 3, nan, nan, nan,
                    nan),
                   ('cf__cdecf__FRAMES', 'cf', 'cf', 'cdecf__FRAMES', 3, 3, nan, nan, nan, nan),
                   ('cde__cdecf__FRAMES', 'cde', 'cde', 'cdecf__FRAMES', 3, 3, nan, nan, nan, nan),
                   ('cdeeg+__FRAMES', 'cdeeg+', 'cdeeg+', 'FRAMES', 3, 3, nan, nan, nan, nan),
                   ('cdecf__FRAMES', 'cdecf', 'cdecf', 'FRAMES', 3, 3, nan, nan, nan, nan),
                   ('cdeeg__cdeeg+__FRAMES', 'cdeeg', 'cdeeg', 'cdeeg+__FRAMES', 3, 3, nan, nan,
                    nan, nan),
                   ('FRAMES', 'FRAMES', 'FRAMES', '', 6, 6, nan, nan, nan, nan),
                   ('c__cf__cdecf__FRAMES', 'c', 'c', 'cf__cdecf__FRAMES', 3, 3, nan, nan, nan,
                    nan),
                   ('c__cde__cdeeg__cdeeg+__FRAMES', 'c', 'c', 'cde__cdeeg__cdeeg+__FRAMES', 3, 3,
                    nan, nan, nan, nan),
                   ('a__ab__FRAMES', 'a', 'a', 'ab__FRAMES', 1, 1, nan, nan, nan, nan)}
        self.assertEqual(lines, w_lines)


class TestAddProportionDataTable(unittest.TestCase):

    @test_for(DataTable.calculate_proportions)
    def test_get_data_proportion_no_relative(self):
        data = DataTable()
        data.fill_parameters(ref_abundance=CT_REF_AB, parent_dict=CT_ONTO,
                             root_item=ROOTS[METACYC], set_abundance=CT_AB, names=CT_LAB)
        data.calculate_proportions(True)
        for i in range(data.len):
            if np.isnan(data.prop[i]):
                self.assertTrue(np.isnan(W_PROP[i]))
            else:
                self.assertEqual(data.prop[i], W_PROP[i])

    @test_for(DataTable.calculate_proportions)
    def test_get_data_proportion_no_relative_ref(self):
        data = DataTable()
        data.fill_parameters(ref_abundance=CT_REF_AB, parent_dict=CT_ONTO,
                             root_item=ROOTS[METACYC], set_abundance=CT_AB, names=CT_LAB)
        data.calculate_proportions(True)
        for i in range(data.len):
            self.assertEqual(data.ref_prop[i], W_REF_PROP[i])

    @test_for(DataTable.calculate_proportions)
    def test_get_data_proportion_relative(self):
        data = DataTable()
        data.fill_parameters(ref_abundance=CT_REF_AB, parent_dict=CT_ONTO,
                             root_item=ROOTS[METACYC], set_abundance=CT_AB, names=CT_LAB)
        data.calculate_proportions(True)
        for k, v in W_REL_PROP.items():
            self.assertEqual(data.relative_prop[data.ids.index(k)], v)


# ENRICHMENT TESTS
# ==================================================================================================

ENRICH_AB = {'00': 50, '01': 5, '02': 25, '03': 20, '04': 1, '05': 5, '06': nan, '07': nan,
             '08': 1, '09': 1}
ENRICH_REF_AB = {'00': 100, '01': 40, '02': 30, '03': 20, '04': 10, '05': 20, '06': 5, '07': 1,
                 '08': 1, '09': 3}
E_LABElS = {'00': '0', '01': '1', '02': '2', '03': '3', '04': '4',
            '05': '5', '06': '6', '07': '7', '08': '8', '09': '9'}
E_ONTO = {'01': ['00'], '02': ['00'], '03': ['00'], '04': ['00'], '05': ['01'],
          '06': ['01'], '07': ['01'], '08': ['02'], '09': ['02']}


# Expected :
# Over : 2, 3 | Under : 1, 4, 5 | No diff : 0, 8, 9 | Nan : 6, 7


class TestEnrichmentAnalysis(unittest.TestCase):

    @test_for(DataTable.make_enrichment_analysis)
    def test_get_data_enrichment_analysis_single_value(self):
        data = DataTable()
        data.fill_parameters(ENRICH_AB, ENRICH_REF_AB, E_ONTO, '00', E_LABElS)
        data.calculate_proportions(True)
        data.make_enrichment_analysis(BINOMIAL_TEST)
        p_value_1 = [data.p_val[i] for i in range(data.len) if data.onto_ids[i] == '01'][0]
        M = 100
        N = 50
        m = 40
        n = 5
        exp_p_value_1 = stats.binomtest(n, N, m / M, alternative='two-sided').pvalue
        exp_p_value_1 = np.log10(exp_p_value_1)
        self.assertEqual(p_value_1, exp_p_value_1)

    @test_for(DataTable.make_enrichment_analysis)
    def test_get_data_enrichment_analysis_binomial(self):
        data = DataTable()
        data.fill_parameters(ENRICH_AB, ENRICH_REF_AB, E_ONTO, '00', E_LABElS)
        data.calculate_proportions(True)
        significant = data.make_enrichment_analysis(BINOMIAL_TEST)
        lines = set(data.get_col())
        exp_significant = {'01': 3.799562441228011e-06, '03': 0.001125114927936431,
                           '02': 0.003092409570144631}
        exp_lines = {('01__00', '01', '1', '00', 5, 40, 0.1, 0.4, 400000, -5.420266413988895),
                     ('02__00', '02', '2', '00', 25, 30, 0.5, 0.3, 300000, 2.509702991379166),
                     ('03__00', '03', '3', '00', 20, 20, 0.4, 0.2, 200000, 2.948803113091024),
                     ('05__01__00', '05', '5', '01__00', 5, 20, 0.1, 0.2, 200000,
                      -1.103304935668835),
                     ('00', '00', '00', '', 50, 100, 1.0, 1.0, 1000000, 0.0),
                     ('04__00', '04', '4', '00', 1, 10, 0.02, 0.1, 100000, -1.2341542222355069),
                     ('06__01__00', '06', '6', '01__00', nan, 5, nan, 0.05, 50000, nan),
                     ('08__02__00', '08', '8', '02__00', 1, 1, 0.02, 0.01, 10000,
                      0.4034095751193356),
                     ('07__01__00', '07', '7', '01__00', nan, 1, nan, 0.01, 10000, nan),
                     ('09__02__00', '09', '9', '02__00', 1, 3, 0.02, 0.03, 30000, 0.0)}
        for line in lines:
            line = tuple([nan if type(x) != str and np.isnan(x) else x for x in line])
            self.assertIn(line, exp_lines)
        self.assertEqual(len(lines), len(exp_lines))
        self.assertEqual(significant, exp_significant)

    @test_for(DataTable.make_enrichment_analysis)
    def test_get_data_enrichment_analysis_hypergeometric(self):
        data = DataTable()
        data.fill_parameters(ENRICH_AB, ENRICH_REF_AB, E_ONTO, '00', E_LABElS)
        data.calculate_proportions(True)
        significant = data.make_enrichment_analysis(HYPERGEO_TEST)
        lines = set(data.get_col())
        exp_lines = {('02__00', '02', '2', '00', 25, 30, 0.5, 0.3, 300000, 4.692610428021241),
                     ('00', '00', '00', '', 50, 100, 1.0, 1.0, 1000000, 0.3010299956639812),
                     ('08__02__00', '08', '8', '02__00', 1, 1, 0.02, 0.01, 10000, -0.0),
                     ('03__00', '03', '3', '00', 20, 20, 0.4, 0.2, 200000, 6.754831139005899),
                     ('05__01__00', '05', '5', '01__00', 5, 20, 0.1, 0.2, 200000,
                      -1.6413993451973743),
                     ('09__02__00', '09', '9', '02__00', 1, 3, 0.02, 0.03, 30000,
                      -1.4464911998299308e-16),
                     ('06__01__00', '06', '6', '01__00', nan, 5, nan, 0.05, 50000, nan),
                     ('01__00', '01', '1', '00', 5, 40, 0.1, 0.4, 400000, -9.138873998573988),
                     ('04__00', '04', '4', '00', 1, 10, 0.02, 0.1, 100000, -1.8051946563380086),
                     ('07__01__00', '07', '7', '01__00', nan, 1, nan, 0.01, 10000, nan)}
        exp_significant = {'01': 7.263166523971598e-10, '03': 1.7586072571039978e-07,
                           '02': 2.0295024128400847e-05}
        for line in lines:
            line = tuple([nan if type(x) != str and np.isnan(x) else x for x in line])
            self.assertIn(line, exp_lines)
        self.assertEqual(len(lines), len(exp_lines))
        self.assertEqual(significant, exp_significant)

    @test_for(DataTable.make_enrichment_analysis)
    def test_get_data_enrichment_analysis_scores(self):
        scores = {'00': 0.05, '01': 0.2, '02': 0.0004, '03': 0.5, '04': 0.000008, '05': 0.9,
                  '06': 0.01, '07': nan, '08': nan, '09': 0.000023}
        data = DataTable()
        data.fill_parameters(ENRICH_AB, ENRICH_REF_AB, E_ONTO, '00', E_LABElS)
        data.calculate_proportions(True)
        significant = data.make_enrichment_analysis(HYPERGEO_TEST, scores)
        lines = set(data.get_col())
        exp_lines = {('05__01__00', '05', '5', '01__00', 5, 20, 0.1, 0.2, 200000,
                      0.045757490560675115),
                     ('02__00', '02', '2', '00', 25, 30, 0.5, 0.3, 300000, 3.3979400086720375),
                     ('01__00', '01', '1', '00', 5, 40, 0.1, 0.4, 400000, 0.6989700043360187),
                     ('04__00', '04', '4', '00', 1, 10, 0.02, 0.1, 100000, 5.096910013008056),
                     ('07__01__00', '07', '7', '01__00', nan, 1, nan, 0.01, 10000, nan),
                     ('00', '00', '00', '', 50, 100, 1.0, 1.0, 1000000, 1.3010299956639813),
                     ('09__02__00', '09', '9', '02__00', 1, 3, 0.02, 0.03, 30000,
                      4.638272163982407),
                     ('03__00', '03', '3', '00', 20, 20, 0.4, 0.2, 200000, 0.3010299956639812),
                     ('06__01__00', '06', '6', '01__00', nan, 5, nan, 0.05, 50000, 2.0),
                     ('08__02__00', '08', '8', '02__00', 1, 1, 0.02, 0.01, 10000, nan)}
        exp_significant = {'04': 8e-06, '09': 2.3e-05, '02': 0.0004}
        for line in lines:
            line = tuple([nan if type(x) != str and np.isnan(x) else x for x in line])
            self.assertIn(line, exp_lines)
        self.assertEqual(len(lines), len(exp_lines))
        self.assertEqual(significant, exp_significant)


# TOPOLOGY MANAGEMENT TESTS
# ==================================================================================================

ROOT_AB = {'R': 50, 'R-1': 50, 'R-2': 50, '00': 50, '01': 5, '02': 25, '03': 20, '04': 1, '05': 5,
           '06': nan, '07': nan, '08': 1, '09': 1}
ROOT_REF_AB = {'R': 100, 'R-1': 100, 'R-2': 100, '00': 100, '01': 40, '02': 30, '03': 20, '04': 10,
               '05': 20, '06': 5, '07': 1, '08': 1, '09': 3}
ROOT_LABElS = {'R': 'r', '00': '0', '01': '1', '02': '2', '03': '3',
               '04': '4', '05': '5', '06': '6', '07': '7', '08': '8', '09': '9'}
ROOT_ONTO = {'01': ['00'], '02': ['00'], '03': ['00'], '04': ['00'], '05': ['01'], '00': ['R-2'],
             '06': ['01'], '07': ['01'], '08': ['02'], '09': ['02'], 'R-1': ['R'], 'R-2': ['R-1']}


PATH_ONTO = {'a': ['ab', 'cdeeg++++'], 'b': ['ab'], 'c': ['cde', 'cf'], 'd': ['cde'], 'e': ['cde', 'eg'],
             'f': ['cf'], 'g': ['gh', 'eg'], 'h': ['gh'], 'ab': [ROOTS[METACYC]],
             'cde': ['cde+'], 'cde+': ['cde++'], 'cde++': ['cde+++'], 'cde+++': ['cdecf', 'cdeeg'],
             'cf': ['cdecf'], 'eg': ['cdeeg', ROOTS[METACYC]], 'gh': [ROOTS[METACYC]],
             'cdecf': [ROOTS[METACYC]],
             'cdeeg': ['cdeeg+'], 'cdeeg+': ['cdeeg++'], 'cdeeg++': ['cdeeg+++'],
             'cdeeg+++': ['cdeeg++++'], 'cdeeg++++': [ROOTS[METACYC]]}
PATH_AB = {'FRAMES': 6, 'cde': 3, 'cde+': 3, 'cde++': 3, 'cde+++': 3, 'cf': 3, 'cdecf': 3,
           'cdeeg++++': 3, 'cdeeg+++': 3, 'cdeeg++': 3, 'cdeeg+': 3, 'cdeeg': 3, 'c': 3, 'ab': 3,
           'b': 2, 'a': 1}
PATH_REF_AB = {'FRAMES': 36, 'cdeeg++++': 19, 'cdeeg+++': 19, 'cdeeg++': 19, 'cdeeg+': 19,
               'cdeeg': 19, 'cdecf': 18, 'gh': 15, 'eg': 12, 'cde': 12, 'cde+': 12, 'cde++': 12,
               'cde+++': 12, 'cf': 9, 'h': 8, 'g': 7, 'f': 6, 'e': 5, 'd': 4, 'c': 3, 'ab': 3,
               'b': 2, 'a': 1}
PATH_LAB = {'FRAMES': 'Root', 'cdeeg+': 'CDEEG+', 'cdeeg': 'CDEEG', 'cdecf': 'CDECF', 'gh': 'GH',
            'eg': 'EG', 'cde': 'CDE', 'cf': 'CF', 'h': 'H', 'g': 'G', 'f': 'F', 'e': 'E', 'd': 'D',
            'c': 'C', 'ab': 'AB', 'b': 'B'}


class TestTopologyManagement(unittest.TestCase):

    @test_for(DataTable.cut_root)
    def test_data_cut_root_uncut(self):
        data = DataTable()
        data.fill_parameters(ROOT_AB, ROOT_REF_AB, ROOT_ONTO, 'R', ROOT_LABElS)
        data.calculate_proportions(True)
        exp_d = data.get_data_dict()
        data.cut_root(ROOT_UNCUT)
        self.assertEqual(data.get_data_dict(), exp_d)

    @test_for(DataTable.cut_root)
    def test_data_cut_root_cut(self):
        data = DataTable()
        data.fill_parameters(ROOT_AB, ROOT_REF_AB, ROOT_ONTO, 'R', ROOT_LABElS)
        data.calculate_proportions(True)
        data.cut_root(ROOT_CUT)
        lines = set(data.get_col())
        exp_lines = {('01__00__R-2__R-1__R', '01', '1', '00', 5, 40, 0.1, 0.4, 400000, nan),
                     (
                         '06__01__00__R-2__R-1__R', '06', '6', '01__00__R-2__R-1__R', nan, 5, nan,
                         0.05,
                         50000, nan),
                     (
                         '07__01__00__R-2__R-1__R', '07', '7', '01__00__R-2__R-1__R', nan, 1, nan,
                         0.01,
                         10000, nan),
                     ('05__01__00__R-2__R-1__R', '05', '5', '01__00__R-2__R-1__R', 5, 20, 0.1, 0.2,
                      200000, nan),
                     ('03__00__R-2__R-1__R', '03', '3', '00', 20, 20, 0.4, 0.2, 200000, nan),
                     ('04__00__R-2__R-1__R', '04', '4', '00', 1, 10, 0.02, 0.1, 100000, nan),
                     ('09__02__00__R-2__R-1__R', '09', '9', '02__00__R-2__R-1__R', 1, 3, 0.02, 0.03,
                      30000, nan),
                     ('02__00__R-2__R-1__R', '02', '2', '00', 25, 30, 0.5, 0.3, 300000, nan),
                     ('08__02__00__R-2__R-1__R', '08', '8', '02__00__R-2__R-1__R', 1, 1, 0.02, 0.01,
                      10000, nan)}
        self.assertEqual(data.len, 9)
        self.assertEqual(len(lines), len(exp_lines))
        for line in lines:
            line = tuple([nan if type(x) != str and np.isnan(x) else x for x in line])
            self.assertIn(line, exp_lines)

    @test_for(DataTable.cut_root)
    def test_data_cut_root_total_cut(self):
        data = DataTable()
        data.fill_parameters(ROOT_AB, ROOT_REF_AB, ROOT_ONTO, 'R', ROOT_LABElS)
        data.calculate_proportions(True)
        data.cut_root(ROOT_TOTAL_CUT)
        lines = set(data.get_col())
        exp_lines = {('01__00__R-2__R-1__R', '01', '1', '', 5, 40, 0.1, 0.4, 400000, nan),
                     (
                         '06__01__00__R-2__R-1__R', '06', '6', '01__00__R-2__R-1__R', nan, 5, nan,
                         0.05,
                         50000, nan),
                     (
                         '07__01__00__R-2__R-1__R', '07', '7', '01__00__R-2__R-1__R', nan, 1, nan,
                         0.01,
                         10000, nan),
                     ('05__01__00__R-2__R-1__R', '05', '5', '01__00__R-2__R-1__R', 5, 20, 0.1, 0.2,
                      200000, nan),
                     ('03__00__R-2__R-1__R', '03', '3', '', 20, 20, 0.4, 0.2, 200000, nan),
                     ('04__00__R-2__R-1__R', '04', '4', '', 1, 10, 0.02, 0.1, 100000, nan),
                     ('09__02__00__R-2__R-1__R', '09', '9', '02__00__R-2__R-1__R', 1, 3, 0.02, 0.03,
                      30000, nan),
                     ('02__00__R-2__R-1__R', '02', '2', '', 25, 30, 0.5, 0.3, 300000, nan),
                     ('08__02__00__R-2__R-1__R', '08', '8', '02__00__R-2__R-1__R', 1, 1, 0.02, 0.01,
                      10000, nan)}
        self.assertEqual(data.len, 9)
        self.assertEqual(len(lines), len(exp_lines))
        for line in lines:
            line = tuple([nan if type(x) != str and np.isnan(x) else x for x in line])
            self.assertIn(line, exp_lines)

    @test_for(DataTable.cut_nested_path)
    def test_cut_path_uncut(self):
        data = DataTable()
        data.fill_parameters(PATH_AB, PATH_REF_AB, PATH_ONTO, ROOTS[METACYC], PATH_LAB)
        data.calculate_proportions(True)
        data.cut_nested_path(PATH_UNCUT, False)
        exp_lines = {('f__cf__cdecf__FRAMES', 'f', 'F', 'cf__cdecf__FRAMES', nan, 6, nan, 0.16666666666666666, 76758, nan), ('cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'cdeeg', 'CDEEG', 'cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 3, 19, 0.5, 0.5277777777777778, 269402, nan), ('cde__cde+__cde++__cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'cde', 'CDE', 'cde+__cde++__cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 3, 12, 0.5, 0.3333333333333333, 134701, nan), ('ab__FRAMES', 'ab', 'AB', 'FRAMES', 3, 3, 0.5, 0.08333333333333333, 44776, nan), ('a__ab__FRAMES', 'a', 'a', 'ab__FRAMES', 1, 1, 0.16666666666666666, 0.027777777777777776, 14925, nan), ('cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'cdeeg+', 'CDEEG+', 'cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 3, 19, 0.5, 0.5277777777777778, 269402, nan), ('g__gh__FRAMES', 'g', 'G', 'gh__FRAMES', nan, 7, nan, 0.19444444444444445, 104477, nan), ('cdecf__FRAMES', 'cdecf', 'CDECF', 'FRAMES', 3, 18, 0.5, 0.5, 268656, nan), ('e__cde__cde+__cde++__cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'e', 'E', 'cde__cde+__cde++__cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', nan, 5, nan, 0.1388888888888889, 56125, nan), ('cde__cde+__cde++__cde+++__cdecf__FRAMES', 'cde', 'CDE', 'cde+__cde++__cde+++__cdecf__FRAMES', 3, 12, 0.5, 0.3333333333333333, 153517, nan), ('cde+++__cdecf__FRAMES', 'cde+++', 'cde+++', 'cdecf__FRAMES', 3, 12, 0.5, 0.3333333333333333, 153517, nan), ('b__ab__FRAMES', 'b', 'B', 'ab__FRAMES', 2, 2, 0.3333333333333333, 0.05555555555555555, 29850, nan), ('h__gh__FRAMES', 'h', 'H', 'gh__FRAMES', nan, 8, nan, 0.2222222222222222, 119402, nan), ('g__eg__FRAMES', 'g', 'G', 'eg__FRAMES', nan, 7, nan, 0.19444444444444445, 104477, nan), ('cde+__cde++__cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'cde+', 'cde+', 'cde++__cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 3, 12, 0.5, 0.3333333333333333, 134701, nan), ('cdeeg+++__cdeeg++++__FRAMES', 'cdeeg+++', 'cdeeg+++', 'cdeeg++++__FRAMES', 3, 19, 0.5, 0.5277777777777778, 269402, nan), ('cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'cdeeg++', 'cdeeg++', 'cdeeg+++__cdeeg++++__FRAMES', 3, 19, 0.5, 0.5277777777777778, 269402, nan), ('FRAMES', 'FRAMES', 'FRAMES', '', 6, 36, 1.0, 1.0, 1000000, nan), ('d__cde__cde+__cde++__cde+++__cdecf__FRAMES', 'd', 'D', 'cde__cde+__cde++__cde+++__cdecf__FRAMES', nan, 4, nan, 0.1111111111111111, 51172, nan), ('d__cde__cde+__cde++__cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'd', 'D', 'cde__cde+__cde++__cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', nan, 4, nan, 0.1111111111111111, 44900, nan), ('e__eg__FRAMES', 'e', 'E', 'eg__FRAMES', nan, 5, nan, 0.1388888888888889, 74626, nan), ('e__eg__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'e', 'E', 'eg__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', nan, 5, nan, 0.1388888888888889, 56125, nan), ('cde+__cde++__cde+++__cdecf__FRAMES', 'cde+', 'cde+', 'cde++__cde+++__cdecf__FRAMES', 3, 12, 0.5, 0.3333333333333333, 153517, nan), ('a__cdeeg++++__FRAMES', 'a', 'a', 'cdeeg++++__FRAMES', 1, 1, 0.16666666666666666, 0.027777777777777776, 14179, nan), ('eg__FRAMES', 'eg', 'EG', 'FRAMES', nan, 12, nan, 0.3333333333333333, 179104, nan), ('c__cde__cde+__cde++__cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'c', 'C', 'cde__cde+__cde++__cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 3, 3, 0.5, 0.08333333333333333, 33675, nan), ('c__cf__cdecf__FRAMES', 'c', 'C', 'cf__cdecf__FRAMES', 3, 3, 0.5, 0.08333333333333333, 38379, nan), ('gh__FRAMES', 'gh', 'GH', 'FRAMES', nan, 15, nan, 0.4166666666666667, 223880, nan), ('e__cde__cde+__cde++__cde+++__cdecf__FRAMES', 'e', 'E', 'cde__cde+__cde++__cde+++__cdecf__FRAMES', nan, 5, nan, 0.1388888888888889, 63965, nan), ('eg__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'eg', 'EG', 'cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', nan, 12, nan, 0.3333333333333333, 134701, nan), ('cdeeg++++__FRAMES', 'cdeeg++++', 'cdeeg++++', 'FRAMES', 3, 19, 0.5, 0.5277777777777778, 283582, nan), ('cde++__cde+++__cdecf__FRAMES', 'cde++', 'cde++', 'cde+++__cdecf__FRAMES', 3, 12, 0.5, 0.3333333333333333, 153517, nan), ('cde++__cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'cde++', 'cde++', 'cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 3, 12, 0.5, 0.3333333333333333, 134701, nan), ('g__eg__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'g', 'G', 'eg__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', nan, 7, nan, 0.19444444444444445, 78575, nan), ('c__cde__cde+__cde++__cde+++__cdecf__FRAMES', 'c', 'C', 'cde__cde+__cde++__cde+++__cdecf__FRAMES', 3, 3, 0.5, 0.08333333333333333, 38379, nan), ('cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'cde+++', 'cde+++', 'cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 3, 12, 0.5, 0.3333333333333333, 134701, nan), ('cf__cdecf__FRAMES', 'cf', 'CF', 'cdecf__FRAMES', 3, 9, 0.5, 0.25, 115138, nan)}
        lines = set(data.get_col())
        self.assertEqual(len(lines), len(exp_lines))
        for line in lines:
            line = tuple([nan if type(x) != str and np.isnan(x) else x for x in line])
            self.assertIn(line, exp_lines)
        # from ontosunburst.sunburst_fig import generate_sunburst_fig
        # generate_sunburst_fig(data, 'test', bg_color='black', font_color='white')

    @test_for(DataTable.cut_nested_path)
    def test_cut_path_cut_deeper(self):
        data = DataTable()
        data.fill_parameters(PATH_AB, PATH_REF_AB, PATH_ONTO, ROOTS[METACYC], PATH_LAB)
        data.calculate_proportions(True)
        data.cut_nested_path(PATH_DEEPER, False)
        # from ontosunburst.sunburst_fig import generate_sunburst_fig
        # generate_sunburst_fig(data, 'test')
        exp_lines = {('e__cde__cde+__cde++__cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'e', 'E', 'cde__cde+__cde++__cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', nan, 5, nan, 0.1388888888888889, 56125, nan), ('eg__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'eg', 'EG', 'cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', nan, 12, nan, 0.3333333333333333, 134701, nan), ('d__cde__cde+__cde++__cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'd', 'D', 'cde__cde+__cde++__cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', nan, 4, nan, 0.1111111111111111, 44900, nan), ('cf__cdecf__FRAMES', 'cf', 'CF', 'cdecf__FRAMES', 3, 9, 0.5, 0.25, 115138, nan), ('g__eg__FRAMES', 'g', 'G', 'eg__FRAMES', nan, 7, nan, 0.19444444444444445, 104477, nan), ('a__cdeeg++++__FRAMES', 'a', 'a', 'cdeeg++++__FRAMES', 1, 1, 0.16666666666666666, 0.027777777777777776, 14179, nan), ('b__ab__FRAMES', 'b', 'B', 'ab__FRAMES', 2, 2, 0.3333333333333333, 0.05555555555555555, 29850, nan), ('e__cde__cde+__cde++__cde+++__cdecf__FRAMES', 'e', 'E', 'cde__cde+__cde++__cde+++__cdecf__FRAMES', nan, 5, nan, 0.1388888888888889, 63965, nan), ('gh__FRAMES', 'gh', 'GH', 'FRAMES', nan, 15, nan, 0.4166666666666667, 223880, nan), ('cde__cde+__cde++__cde+++__cdecf__FRAMES', 'cde', '... CDE', 'cdecf__FRAMES', 3, 12, 0.5, 0.3333333333333333, 153517, nan), ('a__ab__FRAMES', 'a', 'a', 'ab__FRAMES', 1, 1, 0.16666666666666666, 0.027777777777777776, 14925, nan), ('g__eg__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'g', 'G', 'eg__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', nan, 7, nan, 0.19444444444444445, 78575, nan), ('FRAMES', 'FRAMES', 'FRAMES', '', 6, 36, 1.0, 1.0, 1000000, nan), ('f__cf__cdecf__FRAMES', 'f', 'F', 'cf__cdecf__FRAMES', nan, 6, nan, 0.16666666666666666, 76758, nan), ('e__eg__FRAMES', 'e', 'E', 'eg__FRAMES', nan, 5, nan, 0.1388888888888889, 74626, nan), ('c__cde__cde+__cde++__cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'c', 'C', 'cde__cde+__cde++__cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 3, 3, 0.5, 0.08333333333333333, 33675, nan), ('cdecf__FRAMES', 'cdecf', 'CDECF', 'FRAMES', 3, 18, 0.5, 0.5, 268656, nan), ('cde__cde+__cde++__cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'cde', '... CDE', 'cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 3, 12, 0.5, 0.3333333333333333, 134701, nan), ('e__eg__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'e', 'E', 'eg__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', nan, 5, nan, 0.1388888888888889, 56125, nan), ('g__gh__FRAMES', 'g', 'G', 'gh__FRAMES', nan, 7, nan, 0.19444444444444445, 104477, nan), ('cdeeg++++__FRAMES', 'cdeeg++++', 'cdeeg++++', 'FRAMES', 3, 19, 0.5, 0.5277777777777778, 283582, nan), ('c__cde__cde+__cde++__cde+++__cdecf__FRAMES', 'c', 'C', 'cde__cde+__cde++__cde+++__cdecf__FRAMES', 3, 3, 0.5, 0.08333333333333333, 38379, nan), ('d__cde__cde+__cde++__cde+++__cdecf__FRAMES', 'd', 'D', 'cde__cde+__cde++__cde+++__cdecf__FRAMES', nan, 4, nan, 0.1111111111111111, 51172, nan), ('cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'cdeeg', '... CDEEG', 'cdeeg++++__FRAMES', 3, 19, 0.5, 0.5277777777777778, 269402, nan), ('ab__FRAMES', 'ab', 'AB', 'FRAMES', 3, 3, 0.5, 0.08333333333333333, 44776, nan), ('eg__FRAMES', 'eg', 'EG', 'FRAMES', nan, 12, nan, 0.3333333333333333, 179104, nan), ('c__cf__cdecf__FRAMES', 'c', 'C', 'cf__cdecf__FRAMES', 3, 3, 0.5, 0.08333333333333333, 38379, nan), ('h__gh__FRAMES', 'h', 'H', 'gh__FRAMES', nan, 8, nan, 0.2222222222222222, 119402, nan)}
        lines = set(data.get_col())
        self.assertEqual(len(lines), len(exp_lines))
        for line in lines:
            line = tuple([nan if type(x) != str and np.isnan(x) else x for x in line])
            self.assertIn(line, exp_lines)

    @test_for(DataTable.cut_nested_path)
    def test_cut_path_cut_higher(self):
        data = DataTable()
        data.fill_parameters(PATH_AB, PATH_REF_AB, PATH_ONTO, ROOTS[METACYC], PATH_LAB)
        data.calculate_proportions(True)
        data.cut_root(ROOT_CUT)
        data.cut_nested_path(PATH_HIGHER, False)
        exp_lines = {('eg__FRAMES', 'eg', 'EG', 'FRAMES', nan, 12, nan, 0.3333333333333333, 179104, nan), ('d__cde__cde+__cde++__cde+++__cdecf__FRAMES', 'd', 'D', 'cde+++__cdecf__FRAMES', nan, 4, nan, 0.1111111111111111, 51172, nan), ('cdecf__FRAMES', 'cdecf', 'CDECF', 'FRAMES', 3, 18, 0.5, 0.5, 268656, nan), ('d__cde__cde+__cde++__cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'd', 'D', 'cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', nan, 4, nan, 0.1111111111111111, 44900, nan), ('c__cde__cde+__cde++__cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'c', 'C', 'cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 3, 3, 0.5, 0.08333333333333333, 33675, nan), ('cde+++__cdecf__FRAMES', 'cde+++', 'cde+++ ...', 'cdecf__FRAMES', 3, 12, 0.5, 0.3333333333333333, 153517, nan), ('cdeeg+++__cdeeg++++__FRAMES', 'cdeeg+++', 'cdeeg+++ ...', 'cdeeg++++__FRAMES', 3, 19, 0.5, 0.5277777777777778, 269402, nan), ('gh__FRAMES', 'gh', 'GH', 'FRAMES', nan, 15, nan, 0.4166666666666667, 223880, nan), ('cdeeg++++__FRAMES', 'cdeeg++++', 'cdeeg++++', 'FRAMES', 3, 19, 0.5, 0.5277777777777778, 283582, nan), ('cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'cde+++', 'cde+++ ...', 'cdeeg+++__cdeeg++++__FRAMES', 3, 12, 0.5, 0.3333333333333333, 134701, nan), ('g__eg__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'g', 'G', 'eg__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', nan, 7, nan, 0.19444444444444445, 78575, nan), ('a__ab__FRAMES', 'a', 'a', 'ab__FRAMES', 1, 1, 0.16666666666666666, 0.027777777777777776, 14925, nan), ('ab__FRAMES', 'ab', 'AB', 'FRAMES', 3, 3, 0.5, 0.08333333333333333, 44776, nan), ('c__cf__cdecf__FRAMES', 'c', 'C', 'cf__cdecf__FRAMES', 3, 3, 0.5, 0.08333333333333333, 38379, nan), ('c__cde__cde+__cde++__cde+++__cdecf__FRAMES', 'c', 'C', 'cde+++__cdecf__FRAMES', 3, 3, 0.5, 0.08333333333333333, 38379, nan), ('e__eg__FRAMES', 'e', 'E', 'eg__FRAMES', nan, 5, nan, 0.1388888888888889, 74626, nan), ('b__ab__FRAMES', 'b', 'B', 'ab__FRAMES', 2, 2, 0.3333333333333333, 0.05555555555555555, 29850, nan), ('e__cde__cde+__cde++__cde+++__cdecf__FRAMES', 'e', 'E', 'cde+++__cdecf__FRAMES', nan, 5, nan, 0.1388888888888889, 63965, nan), ('g__eg__FRAMES', 'g', 'G', 'eg__FRAMES', nan, 7, nan, 0.19444444444444445, 104477, nan), ('eg__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'eg', 'EG', 'cdeeg+++__cdeeg++++__FRAMES', nan, 12, nan, 0.3333333333333333, 134701, nan), ('a__cdeeg++++__FRAMES', 'a', 'a', 'cdeeg++++__FRAMES', 1, 1, 0.16666666666666666, 0.027777777777777776, 14179, nan), ('g__gh__FRAMES', 'g', 'G', 'gh__FRAMES', nan, 7, nan, 0.19444444444444445, 104477, nan), ('f__cf__cdecf__FRAMES', 'f', 'F', 'cf__cdecf__FRAMES', nan, 6, nan, 0.16666666666666666, 76758, nan), ('cf__cdecf__FRAMES', 'cf', 'CF', 'cdecf__FRAMES', 3, 9, 0.5, 0.25, 115138, nan), ('e__cde__cde+__cde++__cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'e', 'E', 'cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', nan, 5, nan, 0.1388888888888889, 56125, nan), ('h__gh__FRAMES', 'h', 'H', 'gh__FRAMES', nan, 8, nan, 0.2222222222222222, 119402, nan), ('e__eg__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'e', 'E', 'eg__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', nan, 5, nan, 0.1388888888888889, 56125, nan)}
        lines = set(data.get_col())
        self.assertEqual(len(lines), len(exp_lines))
        for line in lines:
            line = tuple([nan if type(x) != str and np.isnan(x) else x for x in line])
            self.assertIn(line, exp_lines)

    @test_for(DataTable.cut_nested_path)
    def test_cut_path_cut_bound(self):
        data = DataTable()
        data.fill_parameters(PATH_AB, PATH_REF_AB, PATH_ONTO, ROOTS[METACYC], PATH_LAB)
        data.calculate_proportions(True)
        data.cut_nested_path(PATH_BOUND, False)
        exp_lines = {('g__eg__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'g', 'G', 'eg__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', nan, 7, nan, 0.19444444444444445, 78575, nan), ('cf__cdecf__FRAMES', 'cf', 'CF', 'cdecf__FRAMES', 3, 9, 0.5, 0.25, 115138, nan), ('d__cde__cde+__cde++__cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'd', 'D', 'cde__cde+__cde++__cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', nan, 4, nan, 0.1111111111111111, 44900, nan), ('gh__FRAMES', 'gh', 'GH', 'FRAMES', nan, 15, nan, 0.4166666666666667, 223880, nan), ('h__gh__FRAMES', 'h', 'H', 'gh__FRAMES', nan, 8, nan, 0.2222222222222222, 119402, nan), ('c__cf__cdecf__FRAMES', 'c', 'C', 'cf__cdecf__FRAMES', 3, 3, 0.5, 0.08333333333333333, 38379, nan), ('cdeeg++++__FRAMES', 'cdeeg++++', 'cdeeg++++', 'FRAMES', 3, 19, 0.5, 0.5277777777777778, 283582, nan), ('eg__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'eg', 'EG', 'cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', nan, 12, nan, 0.3333333333333333, 134701, nan), ('c__cde__cde+__cde++__cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'c', 'C', 'cde__cde+__cde++__cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 3, 3, 0.5, 0.08333333333333333, 33675, nan), ('eg__FRAMES', 'eg', 'EG', 'FRAMES', nan, 12, nan, 0.3333333333333333, 179104, nan), ('cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'cde+++', 'cde+++ ...', 'cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 3, 12, 0.5, 0.3333333333333333, 134701, nan), ('b__ab__FRAMES', 'b', 'B', 'ab__FRAMES', 2, 2, 0.3333333333333333, 0.05555555555555555, 29850, nan), ('cde+++__cdecf__FRAMES', 'cde+++', 'cde+++ ...', 'cdecf__FRAMES', 3, 12, 0.5, 0.3333333333333333, 153517, nan), ('a__ab__FRAMES', 'a', 'a', 'ab__FRAMES', 1, 1, 0.16666666666666666, 0.027777777777777776, 14925, nan), ('ab__FRAMES', 'ab', 'AB', 'FRAMES', 3, 3, 0.5, 0.08333333333333333, 44776, nan), ('e__eg__FRAMES', 'e', 'E', 'eg__FRAMES', nan, 5, nan, 0.1388888888888889, 74626, nan), ('cde__cde+__cde++__cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'cde', '... CDE', 'cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 3, 12, 0.5, 0.3333333333333333, 134701, nan), ('e__cde__cde+__cde++__cde+++__cdecf__FRAMES', 'e', 'E', 'cde__cde+__cde++__cde+++__cdecf__FRAMES', nan, 5, nan, 0.1388888888888889, 63965, nan), ('FRAMES', 'FRAMES', 'FRAMES', '', 6, 36, 1.0, 1.0, 1000000, nan), ('cde__cde+__cde++__cde+++__cdecf__FRAMES', 'cde', '... CDE', 'cde+++__cdecf__FRAMES', 3, 12, 0.5, 0.3333333333333333, 153517, nan), ('cdecf__FRAMES', 'cdecf', 'CDECF', 'FRAMES', 3, 18, 0.5, 0.5, 268656, nan), ('e__cde__cde+__cde++__cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'e', 'E', 'cde__cde+__cde++__cde+++__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', nan, 5, nan, 0.1388888888888889, 56125, nan), ('f__cf__cdecf__FRAMES', 'f', 'F', 'cf__cdecf__FRAMES', nan, 6, nan, 0.16666666666666666, 76758, nan), ('a__cdeeg++++__FRAMES', 'a', 'a', 'cdeeg++++__FRAMES', 1, 1, 0.16666666666666666, 0.027777777777777776, 14179, nan), ('cdeeg+++__cdeeg++++__FRAMES', 'cdeeg+++', 'cdeeg+++ ...', 'cdeeg++++__FRAMES', 3, 19, 0.5, 0.5277777777777778, 269402, nan), ('e__eg__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'e', 'E', 'eg__cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', nan, 5, nan, 0.1388888888888889, 56125, nan), ('cdeeg__cdeeg+__cdeeg++__cdeeg+++__cdeeg++++__FRAMES', 'cdeeg', '... CDEEG', 'cdeeg+++__cdeeg++++__FRAMES', 3, 19, 0.5, 0.5277777777777778, 269402, nan), ('c__cde__cde+__cde++__cde+++__cdecf__FRAMES', 'c', 'C', 'cde__cde+__cde++__cde+++__cdecf__FRAMES', 3, 3, 0.5, 0.08333333333333333, 38379, nan), ('g__gh__FRAMES', 'g', 'G', 'gh__FRAMES', nan, 7, nan, 0.19444444444444445, 104477, nan), ('d__cde__cde+__cde++__cde+++__cdecf__FRAMES', 'd', 'D', 'cde__cde+__cde++__cde+++__cdecf__FRAMES', nan, 4, nan, 0.1111111111111111, 51172, nan), ('g__eg__FRAMES', 'g', 'G', 'eg__FRAMES', nan, 7, nan, 0.19444444444444445, 104477, nan)}
        lines = set(data.get_col())
        self.assertEqual(len(lines), len(exp_lines))
        for line in lines:
            line = tuple([nan if type(x) != str and np.isnan(x) else x for x in line])
            self.assertIn(line, exp_lines)


