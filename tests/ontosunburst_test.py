import unittest
from ontosunburst.ontosunburst import *

"""
Tests manually good file creation.
No automatic tests integrated.
"""

# METACYC
# ==================================================================================================

MET_SET = {'CPD-24674', 'CPD-24687', 'CPD-24688'}
REF_MET = {'CPD-24674', 'CPD-24687', 'CPD-24688',
           'CPD-12782', 'CPD-12784', 'CPD-12787',
           'CPD-12788', 'CPD-12789', 'CPD-12796',
           'CPD-12797', 'CPD-12798', 'CPD-12805',
           'CPD-12806', 'CPD-12812', 'CPD-12816',
           'CPD-1282', 'CPD-12824', 'CPD-1283'}

RXN_SET = {'CROTCOALIG-RXN', 'CYSTHIOCYS-RXN', 'NQOR-RXN'}
REF_RXN = {'CROTCOALIG-RXN', 'CYSTHIOCYS-RXN', 'NQOR-RXN',
           'RXN-14859', 'RXN-14873', 'RXN-14920',
           'RXN-14939', 'RXN-14975', 'RXN-21632',
           'RXN-21638', 'RXN-21652', 'RXN-8954'}

PWY_SET = {'2ASDEG-PWY', '4AMINOBUTMETAB-PWY', 'ALLANTOINDEG-PWY'}
REF_PWY = {'2ASDEG-PWY', '4AMINOBUTMETAB-PWY', 'ALLANTOINDEG-PWY',
           'CRNFORCAT-PWY', 'PWY-7195', 'PWY-7219',
           'PWY-7251', 'PWY-7351', 'PWY-7401',
           'PWY18C3-22', 'PWY0-1600', 'SERDEG-PWY'}


class MetacycTest(unittest.TestCase):
    # Compounds
    def test_cpd_metacyc_proportion(self):
        fig = metacyc_ontosunburst(metabolic_objects=REF_MET, output='test_mc_cpd_prop')
        fig.write_image('test_mc_cpd_prop.png', width=1900, height=1000, scale=1)

    def test_cpd_metacyc_comparison(self):
        fig = metacyc_ontosunburst(metabolic_objects=MET_SET, reference_set=REF_MET,
                                   output='test_mc_cpd_comp')
        fig.write_image('test_mc_cpd_comp.png', width=1900, height=1000, scale=1)
    # Reactions
    def test_rxn_metacyc_proportion(self):
        fig = metacyc_ontosunburst(metabolic_objects=REF_RXN, output='test_mc_rxn_prop')
        fig.write_image('test_mc_rxn_prop.png', width=1900, height=1000, scale=1)

    def test_rxn_metacyc_comparison(self):
        fig = metacyc_ontosunburst(metabolic_objects=RXN_SET, reference_set=REF_RXN,
                                   output='test_mc_rxn_comp')
        fig.write_image('test_mc_rxn_comp.png', width=1900, height=1000, scale=1)

    # Pathways
    def test_pwy_metacyc_proportion(self):
        fig = metacyc_ontosunburst(metabolic_objects=REF_PWY, output='test_mc_pwy_prop')
        fig.write_image('test_mc_pwy_prop.png', width=1900, height=1000, scale=1)

    def test_pwy_metacyc_comparison(self):
        fig = metacyc_ontosunburst(metabolic_objects=PWY_SET, reference_set=REF_PWY,
                                   output='test_mc_pwy_comp')
        fig.write_image('test_mc_pwy_comp.png', width=1900, height=1000, scale=1)


# EC
# ==================================================================================================

EC_SET = {'2.6.1.45', '1.1.1.25', '1.1.1.140'}
REF_EC = {'2.6.1.45', '1.1.1.25', '1.1.1.140',
          '1.14.14.52', '2.7.1.137', '7.1.1.8',
          '1.17.4.5', '2.3.1.165', '3.2.1.53',
          '3.2.1.91', '6.3.4.2', '5.4.99.8'}


class EcTest(unittest.TestCase):

    def test_ec_proportion(self):
        fig = ec_ontosunburst(ec_set=REF_EC, output='test_ec_prop')
        fig.write_image('test_ec_prop.png', width=1900, height=1000, scale=1)

    def test_ec_comparison(self):
        fig = ec_ontosunburst(ec_set=EC_SET, reference_set=REF_EC, output='test_ec_comp')
        fig.write_image('test_ec_comp.png', width=1900, height=1000, scale=1)


# CHEBI
# ==================================================================================================

URL = 'http://localhost:3030/chebi/'
CH_SET = {'38028', '28604', '85146'}
REF_CH = {'38028', '28604', '85146',
          '23066', '27803', '37565',
          '58215', '79983', '42639'}


class ChEbiTest(unittest.TestCase):

    def test_chebi_proportion(self):
        fig = chebi_ontosunburst(chebi_ids=REF_CH, endpoint_url=URL, output='test_chebi_prop')
        fig.write_image('test_chebi_prop.png', width=1900, height=1000, scale=1)

    def test_chebi_comparison(self):
        fig = chebi_ontosunburst(chebi_ids=CH_SET, reference_set=REF_CH, endpoint_url=URL,
                                 output='test_chebi_comp')
        fig.write_image('test_chebi_comp.png', width=1900, height=1000, scale=1)
