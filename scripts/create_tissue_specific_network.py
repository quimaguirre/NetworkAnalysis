import argparse
import cPickle
import numpy as np
import pandas as pd
import time
import sys, os, re

from context import NetworkAnalysis
import NetworkAnalysis.network_analysis as NA
import NetworkAnalysis.tissue_specificity as TS




def main():

    options = parse_user_arguments()
    create_tissue_specific_network(options)

def parse_user_arguments(*args, **kwds):
    """
    Parses the arguments of the program
    """

    parser = argparse.ArgumentParser(
        description = "Create a tissue-specific network",
        epilog      = "@oliva's lab 2017")
    parser.add_argument('-n','--network_file',dest='network_file',action = 'store',
                        help = """" Input file with a protein-protein interaction network in SIF format. """)
    parser.add_argument('-nf','--network_format',dest='network_format',action = 'store',default='multi-fields',
                        help = '''Format of the edge file (network):\tsif, netscore, raw, multi-fields (default):\n
                                    'sif': <node1>\tscore\t<node2>\n
                                    'netscore': <node1>\t<node2>\t<score>\n
                                    'raw': <node1>\t<node2>\n
                                    'multi-fields' : <node1>\t<node2>\t<sources>\t<method_ids>\t<method_names>\t<pmids>\n''')
    parser.add_argument('-t','--tissue',dest='tissue',action = 'store',
                        help = """" Name of the tissue (liver, kidney, brain, heart or pancreas) """)
    parser.add_argument('-trans','--translation_file',dest='translation_file',action = 'store',
                        help = """" Input file with the translation file of biana codes to geneID """)
    parser.add_argument('-ws','--workspace',dest='workspace',action = 'store',default=os.path.join(os.path.join(os.path.dirname(__file__), '..'), 'workspace'),
                        help = """Define the workspace directory where the results will be created. """)
    parser.add_argument('-d','--data_dir',dest='data_dir',action = 'store',default=os.path.join(os.path.join(os.path.dirname(__file__), '..'), 'data'),
                        help = """Define the data directory where the databases are stored. """)

    options=parser.parse_args()

    return options


#################
#################
# MAIN FUNCTION #
#################
#################

def create_tissue_specific_network(options):
    """
    Generates the profiles of the input drug
    """

    # Start marker for time measure
    start = time.time()


    #----------------------#
    #   DEFINE THE PATHS   #
    #----------------------#

    # Get the program path
    main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    data_dir = workspace.data_dir


    #-----------------------------#
    #   PARSE THE INPUT NETWORK   #
    #-----------------------------#

    # Define the main network
    network_file = options.network_file
    type_id = 'biana'
    network_format = options.network_format
    network_instance = NA.Network(network_file, type_id, network_format)

    # Analyze the main network
    print('BIANA main network')
    print('Number of edges: {}'.format(len(network_instance.get_edges())))
    print('Number of nodes: {}\n'.format(len(network_instance.get_nodes())))

    methodid2interactions = network_instance.get_methodid2interactions() # Get the methods in the network and the number of interactions per method
    network_instance.get_numpmids2interactions() # Get the number of interactions in which we have a given number of pubmed IDs
    network_instance.get_database2interactions() # Get the database and the number of interactions


    #---------------------------#
    #   TRANSLATE THE NETWORK   #
    #---------------------------#

    translation_file = options.translation_file
    translation_id = 'geneid'
    translated_network = network_file+'.geneid'
    translated_network_instance = network_instance.translate_network(translation_file, translation_id, network_format, translated_network, None)


    #-----------------------------------#
    #   DEFINE THE HOUSEKEEPING GENES   #
    #-----------------------------------#

    # Define the housekeeping genes object
    translation_id = 'geneid'
    hk_genes = TS.HouseKeepingGenes(translation_id, pickles_dir)
    all_hk_genes = hk_genes.all_hk_genes
    hpa_hk_genes = hk_genes.hpa_hk_genes
    eisenberg_hk_genes = hk_genes.eisenberg_hk_genes


    #------------------------------------------------------------#
    #   CALCULATE THE COVERAGE OF HK GENES IN THE MAIN NETWORK   #
    #------------------------------------------------------------#

    # Calculate coverage of house-keeping genes in the main network
    num_all_hk_network = float(translated_network_instance.get_number_of_housekeeping_genes(all_hk_genes))
    num_hpa_hk_network = float(translated_network_instance.get_number_of_housekeeping_genes(hpa_hk_genes))
    num_eis_hk_network = float(translated_network_instance.get_number_of_housekeeping_genes(eisenberg_hk_genes))
    num_all_nonhk_network = float(len(set(translated_network_instance.get_nodes()))) - num_all_hk_network
    per_all_hk_network = num_all_hk_network / float(len(all_hk_genes)) * 100
    per_hpa_hk_network = num_hpa_hk_network / float(len(hpa_hk_genes)) * 100
    per_eis_hk_network = num_eis_hk_network / float(len(eisenberg_hk_genes)) * 100

    print('GeneID main network')
    print('Number of edges: {}'.format(len(translated_network_instance.get_edges())))
    print('Number of nodes: {}'.format(len(translated_network_instance.get_nodes())))
    print('Percentage of (all) housekeeping genes in the main network: {:.2f}%\t{:.0f} of {} genes'.format(per_all_hk_network, num_all_hk_network, len(all_hk_genes)))
    print('Percentage of (HPA) housekeeping genes in the main network: {:.2f}%\t{:.0f} of {} genes'.format(per_hpa_hk_network, num_hpa_hk_network, len(hpa_hk_genes)))
    print('Percentage of (Eisenberg) housekeeping genes in the main network: {:.2f}%\t{:.0f} of {} genes\n'.format(per_eis_hk_network, num_eis_hk_network, len(eisenberg_hk_genes)))


    #-----------------------#
    #   DEFINE THE TISSUE   #
    #-----------------------#

    tissue_instance = NA.Tissue(tissue_terms[options.tissue.lower()]['hpa'], tissue_terms[options.tissue.lower()]['jensen'], jensen_conf=3, hpa_level='low', hpa_rel='approved', pickles_path=pickles_dir)


    #----------------------------------------#
    #   CREATE THE TISSUE-SPECIFIC NETWORK   #
    #----------------------------------------#

    tissue_network_file = os.path.join(options.workspace, '{}_edges.biana.txt'.format(options.tissue))
    permission = 0
    tissue_network_instance = network_instance.filter_network_by_tissue(filtered_network_file=tissue_network_file, tissue_object=tissue_instance, permission=permission, filtered_nodes_file=None)

    # Analyze the network
    print('BIANA {} network'.format(options.tissue))
    print('Number of edges: {}'.format(len(tissue_network_instance.get_edges())))
    print('Number of nodes: {}\n'.format(len(tissue_network_instance.get_nodes())))

    # Translate the tissue-specific network
    tissue_network_file = os.path.join(options.workspace, '{}_edges.{}.txt'.format(options.tissue, translation_id))
    tissue_network_translated = tissue_network_instance.translate_network(translation_file, translation_id, network_format, tissue_network_file)


    #------------------------------------#
    #   DEFINE THE BENCHMARK FOR LIVER   #
    #------------------------------------#

    # Define the Wang liver network
    wang_liver_file = os.path.join(pickles_dir, 'wang_liver_network.pcl')
    wang_liver_network = cPickle.load(open(wang_liver_file))
    wang_edges_file = os.path.join(options.workspace, 'wang_liver_edges.geneid.txt')
    NA.write_network_file_from_networkx_graph(wang_liver_network, wang_edges_file, 'raw')
    wang_network = NA.Network(wang_edges_file, 'geneid', 'raw')


    #--------------------------------------------------------------#
    #   CALCULATE THE COVERAGE OF HK GENES IN THE TISSUE NETWORK   #
    #--------------------------------------------------------------#

    # Get info from tissue-specific network
    hpa_edges=tissue_network_translated.hpa_edges
    jensen_edges=tissue_network_translated.jensen_edges
    union=tissue_network_translated.get_union()
    intersection=tissue_network_translated.get_intersection()

    # Calculate coverage of housekeeping genes in the tissue-specific network
    num_all_hk_tissue = float(tissue_network_translated.get_number_of_housekeeping_genes(all_hk_genes))
    num_hpa_hk_tissue = float(tissue_network_translated.get_number_of_housekeeping_genes(hpa_hk_genes))
    num_eis_hk_tissue = float(tissue_network_translated.get_number_of_housekeeping_genes(eisenberg_hk_genes))
    num_all_nonhk_tissue = float(len(set(tissue_network_translated.get_nodes()))) - num_all_hk_tissue
    per_all_hk_tissue = num_all_hk_tissue / num_all_hk_network * 100
    per_hpa_hk_tissue = num_hpa_hk_tissue / num_hpa_hk_network * 100
    per_eis_hk_tissue = num_eis_hk_tissue / num_eis_hk_network * 100

    contingency_table = np.array([[num_all_hk_tissue, num_all_nonhk_tissue], [num_all_hk_network, num_all_nonhk_network]])
    chi2, pval, dof, expected = NA.calculate_contingency_table(contingency_table)

    print('GeneID {}-specific network'.format(options.tissue))
    print('Number of edges: {}'.format(len(tissue_network_translated.get_edges())))
    print('Number of nodes: {}'.format(len(tissue_network_translated.get_nodes())))
    print('Interactions using only Tissues (Jensen lab): {}'.format(len(jensen_edges)))
    print('Interactions using only Human Protein Atlas: {}'.format(len(hpa_edges)))
    print('Intersection: {}'.format(len(intersection)))
    print('Union: {}\n'.format(len(union)))
    print('Percentage of (all) HK genes in the tissue with respect to main network: {:.2f}%\t{:.0f} of {:.0f} genes'.format(per_all_hk_tissue, num_all_hk_tissue, num_all_hk_network))
    print('Percentage of (HPA) HK genes in the tissue with respect to main network: {:.2f}%\t{:.0f} of {:.0f} genes'.format(per_hpa_hk_tissue, num_hpa_hk_tissue, num_hpa_hk_network))
    print('Percentage of (Eisenberg) HK genes in the tissue with respect to main network: {:.2f}%\t{:.0f} of {:.0f} genes'.format(per_eis_hk_tissue, num_eis_hk_tissue, num_eis_hk_network))
    print('Contingency table: {}'.format(contingency_table))
    print('Contingency table result:\tchi2: {}\tp-value: {}\tdegrees of freedom: {}\n'.format(chi2, pval, dof))


    #--------------------------------------------------------------#
    #   CALCULATE THE COVERAGE OF LIVER GENES USING WANG NETWORK   #
    #--------------------------------------------------------------#

    if options.tissue == 'liver':

        # Calculate coverage of Wang liver network genes in the main network
        # and in our tissue-specific netwrok
        wang_nodes_intersection_network = NA.get_nodes_intersection_of_two_networks(translated_network_instance, wang_network)
        wang_nodes_intersection_liver = NA.get_nodes_intersection_of_two_networks(tissue_network_translated, wang_network)
        num_wang_network = float(len(wang_nodes_intersection_network))
        num_nonwang_network = float(len(set(translated_network_instance.get_nodes()))) - num_wang_network
        num_wang_liver = float(len(wang_nodes_intersection_liver))
        num_nonwang_liver = float(len(set(tissue_network_translated.get_nodes()))) - num_wang_liver
        per_wang_network = num_wang_network / float(len(wang_network.get_nodes())) * 100
        per_wang_tissue = num_wang_liver / num_wang_network * 100

        contingency_table = np.array([[num_wang_liver, num_nonwang_liver], [num_wang_network, num_nonwang_network]])
        chi2, pval, dof, expected = NA.calculate_contingency_table(contingency_table)
        print('Percentage of Wang liver genes in the main network: {:.2f}%\t{:.0f} of {:.0f} genes'.format(per_wang_network, num_wang_network, float(len(wang_network.get_nodes()))))
        print('Percentage of Wang liver genes in the liver-specific network with respect to main network: {:.2f}%\t{:.0f} of {:.0f} genes'.format(per_wang_tissue, num_wang_liver, num_wang_network))
        print('Contingency table: {}'.format(contingency_table))
        print('Contingency table result:\tchi2: {}\tp-value: {}\tdegrees of freedom: {}\n'.format(chi2, pval, dof))


        # Calculate coverage of Wang liver network interactions in the main
        # network and in our tissue-specific netwrok
        wang_intersection_network = NA.get_edges_intersection_of_two_networks(translated_network_instance, wang_network)
        wang_intersection_liver = NA.get_edges_intersection_of_two_networks(tissue_network_translated, wang_network)
        num_wang_network = float(len(wang_intersection_network))
        num_nonwang_network = float(len(set(translated_network_instance.get_edges()))) - num_wang_network
        num_wang_liver = float(len(wang_intersection_liver))
        num_nonwang_liver = float(len(set(tissue_network_translated.get_edges()))) - num_wang_liver
        per_wang_network = num_wang_network / float(len(wang_network.get_edges())) * 100
        per_wang_tissue = num_wang_liver / num_wang_network * 100

        contingency_table = np.array([[num_wang_liver, num_nonwang_liver], [num_wang_network, num_nonwang_network]])
        chi2, pval, dof, expected = NA.calculate_contingency_table(contingency_table)
        print('Percentage of Wang liver interactions in the main network: {:.2f}%\t{:.0f} of {:.0f} edges'.format(per_wang_network, num_wang_network, float(len(wang_network.get_edges()))))
        print('Percentage of Wang liver interactions in the liver-specific network with respect to main network: {:.2f}%\t{:.0f} of {:.0f} edges'.format(per_wang_tissue, num_wang_liver, num_wang_network))
        print('Contingency table: {}'.format(contingency_table))
        print('Contingency table result:\tchi2: {}\tp-value: {}\tdegrees of freedom: {}\n'.format(chi2, pval, dof))


    # End marker for time
    end = time.time()
    print('\nTIME OF EXECUTION: {:.3f} seconds or {:.3f} minutes.\n'.format(end - start, (end - start) / 60))


    return




#######################
#######################
# SECONDARY FUNCTIONS #
#######################
#######################



if  __name__ == "__main__":
    main()

