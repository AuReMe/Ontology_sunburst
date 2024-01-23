import csv
import json
import os.path

from padmet.utils.sbmlPlugin import convert_from_coded_id
from padmet.classes import PadmetSpec
from typing import List, Set, Tuple


# METABOLITES EXTRACTION
# ==================================================================================================

# From ClusterMap output
def extract_metabolites_clusters(input_file: str, clust: List[int]) -> Set[str]:
    """ Extract metabolites corresponding to clusters desired.

    Parameters
    ----------
    input_file: str
        Cluster tsv file from heatmap output
    clust: List[int]
        List of cluster numbers to extract

    Returns
    -------
    Set[str]
        Set of metabolites corresponding to cluster numbers.
    """
    with open(input_file, 'r') as f:
        rows = list(csv.reader(f, delimiter='\t'))
        del rows[0]
        metabolites = {row[1].replace('_C-BOUNDARY', '') for row in rows if int(row[0]) in clust}
    return metabolites


# From M2M MetaCom output
def extract_community_metabolites_m2m(input_directory: str) -> Tuple[Set[str], Set[str]]:
    """

    Parameters
    ----------
    input_directory: str
        Directory path to m2m results, must contain 'community_analysis' folder

    Returns
    -------
    Tuple[Set[str], Set[str]]
        Set of metabolites from Added Value
        Set of metabolites from Community Scopes
    """
    added_value_json = os.path.join(input_directory, 'community_analysis', 'addedvalue.json')
    comm_scopes_json = os.path.join(input_directory, 'community_analysis', 'comm_scopes.json')
    with open(added_value_json, 'r') as av_f, open(comm_scopes_json, 'r') as cs_f:
        av_set = {convert_from_coded_id(x)[0].replace('_C-BOUNDARY', '') for x in set(json.load(av_f)['addedvalue'])}
        cs_set = {convert_from_coded_id(x)[0].replace('_C-BOUNDARY', '') for x in set(json.load(cs_f)['com_scope'])}
    return av_set, cs_set


def extract_host_metabolites_m2m(input_directory: str) -> Set[str]:
    """

    Parameters
    ----------
    input_directory: str
        Directory path to m2m results, must contain 'community_analysis' folder

    Returns
    -------
    Set[str]
        Set of the host scope metabolites
    """
    comm_scopes_json = os.path.join(input_directory, 'community_analysis', 'comm_scopes.json')
    with open(comm_scopes_json, 'r') as h_f:
        h_set = {convert_from_coded_id(x)[0].replace('_C-BOUNDARY', '') for x in set(json.load(h_f)['host_scope'])}
    return h_set


# From Padmet Network
def extract_network_metabolites(network_file: str, consumes: bool = False) -> Set[str]:
    """

    Parameters
    ----------
    network_file: str
        Path to padmet network file.
    consumes: bool

    Returns
    -------
    Set[str]
        Set of metabolites from the metabolic network.
    """
    padmet = PadmetSpec(network_file)
    test_set = ['produces']
    if consumes:
        test_set.append('consumes')
    me = {rel.id_out for rel in padmet.getAllRelation() if rel.type in test_set}
    me.remove('Bio')
    return me


# From i_scope group selection
def extract_iscope_group_metabolites(scope_json, g_file, group_list):
    with open(scope_json, 'r') as f:
        scopes = json.load(f)
        for sp, met_l in scopes.items():
            scopes[sp] = set(met_l)
    # groups = extract_groups(g_file)
    union_scope = set()
    # for grp in group_list:
    #     sp_list = groups[grp]
    #     all_metabolites = [met for sp, met in scopes.items() if sp in sp_list]
    #     if all_metabolites:
    #         union_scope = union_scope.union(*all_metabolites)
    return set([convert_from_coded_id(met)[0].replace('_C-BOUNDARY', '') for met in union_scope])


# PATHWAYS EXTRACTION
# ==================================================================================================

# From padmet network
# def extract_pathway_classes(padmet_networks_list: List[str], completion_threshold: float = 0) -> Dict[str, List[str]]:
#     """
#
#     Parameters
#     ----------
#     padmet_networks_list: List[str]
#         List of padmet networks path
#     completion_threshold: float (default=0)
#         Minimal completion rate of the pathway to filter.
#
#     Returns
#     -------
#     Dict[str, List[str]]
#         Dict[pathway_ID, List[pathway_class]] : associate for each pathway, its classes associated.
#     """
#     pw_classes = dict()
    # for species_nw in padmet_networks_list:
    #     sp_nw = PadmetNetwork(species_nw)
    #     for pw_id, pw in sp_nw.pathways.items():
    #         if pw.completion_rate > completion_threshold:
    #             pw_classes[pw_id] = pw.is_class
    # return pw_classes


# GO TERMS EXTRACTION
# ==================================================================================================

def extract_go_esmecata(sp_annotations):
    go_abundance = dict()
    for sp_file in sp_annotations:
        with open(sp_file, 'r') as f:
            rows = csv.reader(f, delimiter='\t')
            rows.__next__()
            for row in rows:
                go_terms = row[3].split(',')
                for go in go_terms:
                    if go not in go_abundance:
                        go_abundance[go] = 0
                    go_abundance[go] += 1
    return go_abundance


# EC EXTRACTION
# ==================================================================================================

def extract_ec_esmecata(sp_annotations):
    ec_set = set()
    for sp_file in sp_annotations:
        with open(sp_file, 'r') as f:
            rows = csv.reader(f, delimiter='\t')
            rows.__next__()
            for row in rows:
                ec_number = row[4].split(',')
                for ec in ec_number:
                    if ec != '':
                        ec_set.add(ec)
    return ec_set
