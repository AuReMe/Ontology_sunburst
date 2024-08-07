import unittest
from ontosunburst.ontosunburst import *
"""
Tests manually good file creation.
No automatic tests integrated.
"""

# METACYC
# ==================================================================================================

MET_LST = ['CPD-24674', 'CPD-24687', 'CPD-24688']
MET_REF = ['CPD-24674', 'CPD-24687', 'CPD-24688',
           'CPD-12782', 'CPD-12784', 'CPD-12787',
           'CPD-12788', 'CPD-12789', 'CPD-12796',
           'CPD-12797', 'CPD-12798', 'CPD-12805',
           'CPD-12806', 'CPD-12812', 'CPD-12816',
           'CPD-1282', 'CPD-12824', 'CPD-1283']
MET_LAB = [1, 2, 3]
MET_RAB = [1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9]

RXN_LST = ['CROTCOALIG-RXN', 'CYSTHIOCYS-RXN', 'NQOR-RXN']
RXN_REF = ['CROTCOALIG-RXN', 'CYSTHIOCYS-RXN', 'NQOR-RXN',
           'RXN-14859', 'RXN-14873', 'RXN-14920',
           'RXN-14939', 'RXN-14975', 'RXN-21632',
           'RXN-21638', 'RXN-21652', 'RXN-8954']

PWY_LST = ['2ASDEG-PWY', '4AMINOBUTMETAB-PWY', 'ALLANTOINDEG-PWY']
PWY_REF = ['2ASDEG-PWY', '4AMINOBUTMETAB-PWY', 'ALLANTOINDEG-PWY',
           'CRNFORCAT-PWY', 'PWY-7195', 'PWY-7219',
           'PWY-7251', 'PWY-7351', 'PWY-7401',
           'PWY18C3-22', 'PWY0-1600', 'SERDEG-PWY']


class MetacycTest(unittest.TestCase):
    # Compounds
    def test_cpd_metacyc_proportion(self):
        fig = ontosunburst(ontology=METACYC, interest_set=MET_LST,
                           output='test_mc_cpd_prop', ref_base=False, abundances=MET_LAB,
                           show_leaves=True)
        fig.write_image('test_mc_cpd_prop.png', width=1900, height=1000, scale=1)

    def test_cpd_metacyc_comparison(self):
        fig = ontosunburst(ontology=METACYC, interest_set=MET_LST, reference_set=MET_REF,
                           analysis=ENRICHMENT_A, output='test_mc_cpd_comp', ref_base=True)
        fig.write_image('test_mc_cpd_comp.png', width=1900, height=1000, scale=1)

    # Reactions
    def test_rxn_metacyc_proportion(self):
        fig = ontosunburst(ontology=METACYC, interest_set=RXN_LST, reference_set=RXN_REF,
                           output='test_mc_rxn_prop', ref_base=True)
        fig.write_image('test_mc_rxn_prop.png', width=1900, height=1000, scale=1)

    def test_rxn_metacyc_comparison(self):
        fig = ontosunburst(ontology=METACYC, interest_set=RXN_LST, reference_set=RXN_REF,
                           analysis=ENRICHMENT_A, output='test_mc_rxn_comp', ref_base=True)
        fig.write_image('test_mc_rxn_comp.png', width=1900, height=1000, scale=1)

    # Pathways
    def test_pwy_metacyc_proportion(self):
        fig = ontosunburst(ontology=METACYC, interest_set=PWY_LST, reference_set=PWY_REF,
                           output='test_mc_pwy_prop', ref_base=True)
        fig.write_image('test_mc_pwy_prop.png', width=1900, height=1000, scale=1)

    def test_pwy_metacyc_comparison(self):
        fig = ontosunburst(ontology=METACYC, interest_set=PWY_LST, reference_set=PWY_REF,
                           analysis=ENRICHMENT_A, output='test_mc_pwy_comp', ref_base=True)
        fig.write_image('test_mc_pwy_comp.png', width=1900, height=1000, scale=1)


# EC
# ==================================================================================================

EC_LST = ['2.6.1.45', '1.1.1.25', '1.1.1.140', '1.1.2.-']
REF_EC = ['2.6.1.45', '1.1.1.25', '1.1.1.140',
          '1.14.14.52', '2.7.1.137', '7.1.1.8',
          '1.17.4.5', '2.3.1.165', '3.2.1.53',
          '3.2.1.91', '6.3.4.2', '5.4.99.8']


class EcTest(unittest.TestCase):

    def test_ec_proportion(self):
        fig = ontosunburst(ontology=EC, interest_set=EC_LST, reference_set=REF_EC,
                           output='test_ec_prop', ref_base=False, show_leaves=True)
        fig.write_image('test_ec_prop.png', width=1900, height=1000, scale=1)

    def test_ec_comparison(self):
        fig = ontosunburst(ontology=EC, interest_set=EC_LST, reference_set=REF_EC,
                           output='test_ec_comp', analysis=ENRICHMENT_A, ref_base=True)
        fig.write_image('test_ec_comp.png', width=1900, height=1000, scale=1)


# CHEBI
# ==================================================================================================

URL = 'http://localhost:3030/chebi/'
CH_LST = ['38028', '28604', '85146']
REF_CH = ['38028', '28604', '85146',
          '23066', '27803', '37565',
          '58215', '79983', '42639']


class ChEbiTest(unittest.TestCase):

    def test_chebi_proportion(self):
        fig = ontosunburst(ontology=CHEBI, interest_set=CH_LST, reference_set=REF_CH,
                           endpoint_url=URL, output='test_chebi_prop', ref_base=True)
        fig.write_image('test_chebi_prop.png', width=1900, height=1000, scale=1)

    def test_chebi_comparison(self):
        fig = ontosunburst(ontology=CHEBI, interest_set=CH_LST, reference_set=REF_CH,
                           endpoint_url=URL, output='test_chebi_comp', analysis=ENRICHMENT_A,
                           ref_base=True)
        fig.write_image('test_chebi_comp.png', width=1900, height=1000, scale=1)


# GO
# ==================================================================================================

GO_EX = ['go:0043226', 'go:0043227', 'go:0043229', 'go:0043231', 'go:0044422', 'go:0044424',
         'go:0044429', 'go:0044444']


class GOTest(unittest.TestCase):

    def test_go_proportion(self):
        fig = ontosunburst(ontology=GO, interest_set=GO_EX,
                           output='test_go_prop', ref_base=False, show_leaves=True)
        fig.write_image('test_go_prop.png', width=1900, height=1000, scale=1)


# KEGG
# ==================================================================================================

KEGG_EX = ['M00572', 'M00308', 'M00844', 'M00633', 'M00176', 'M00535', 'M00573', 'M00970', 'M00131',
           'M00620']


class KeggTest(unittest.TestCase):

    def test_kegg_proportion(self):
        fig = ontosunburst(ontology=GO, interest_set=GO_EX,
                           output='test_go_prop', ref_base=False, show_leaves=True)
        fig.write_image('test_kegg_prop.png', width=1900, height=1000, scale=1)
