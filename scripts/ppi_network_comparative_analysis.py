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
    ppi_network_comparative_analysis(options)

def parse_user_arguments(*args, **kwds):
    """
    Parses the arguments of the program
    """

    parser = argparse.ArgumentParser(
        description = "Comparison of PPI networks",
        epilog      = "@oliva's lab 2017")
    parser.add_argument('-n','--network_file',dest='network_file',action = 'store',
                        help = """" Input file with a protein-protein interaction network in SIF format. """)
    parser.add_argument('-nf','--network_format',dest='network_format',action = 'store',default='multi-fields',
                        help = '''Format of the edge file (network):\tsif, netscore, raw, multi-fields (default):\n
                                    'sif': <node1>\tscore\t<node2>\n
                                    'netscore': <node1>\t<node2>\t<score>\n
                                    'raw': <node1>\t<node2>\n
                                    'multi-fields' : <node1>\t<node2>\t<sources>\t<method_ids>\t<method_names>\t<pmids>\n''')
    parser.add_argument('-trans','--translation_file',dest='translation_file',action = 'store',
                        help = """" Input file with the translation file of biana codes to geneID """)
    parser.add_argument('-ws','--workspace',dest='workspace',action = 'store',default=os.path.join(os.path.join(os.path.dirname(__file__), '..'), 'workspace'),
                        help = """Define the workspace directory where the results will be created. """)

    options=parser.parse_args()

    return options


#################
#################
# MAIN FUNCTION #
#################
#################

def ppi_network_comparative_analysis(options):
    """
    Comparative analysis of different PPI networks.
    """

    # Start marker for time measure
    start = time.time()


    #----------------------#
    #   DEFINE THE PATHS   #
    #----------------------#

    # Get the program path
    main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    pickles_dir = os.path.join(main_path, 'NetworkAnalysis/pickles')
    main_dir = os.path.join(options.workspace, 'main_networks') # Path to store the main networks
    create_directory(main_dir)
    tissue_spec_dir = os.path.join(options.workspace, 'tissue_specific_networks') # Path to store the tissue specific networks
    create_directory(tissue_spec_dir)

    # Define the results variable
    results_table = {}
    output_table_file = os.path.join(options.workspace, 'tissue_specific_results.tsv')



    ###################################
    ## CREATION OF THE MAIN NETWORKS ##
    ###################################

    # Define the main network
    network_file = options.network_file
    type_id = 'biana'
    network_format = options.network_format
    main_network = NA.Network(network_file, type_id, network_format)

    # Analyze the main network
    print('BIANA main network')
    print('Number of edges: {}'.format(len(main_network.get_edges())))
    print('Number of nodes: {}\n'.format(len(main_network.get_nodes())))
    methodid2interactions = main_network.get_methodid2interactions() # Get the methods in the network and the number of interactions per method
    main_network.get_numpmids2interactions() # Get the number of interactions in which we have a given number of pubmed IDs
    main_network.get_database2interactions() # Get the database and the number of interactions


    # Create a network filtered by methods
    psimi2score_file = os.path.join(pickles_dir, 'psimi2score.pcl')
    psimi2score = cPickle.load(open(psimi2score_file))
    method_ids_excluded = []
    method_ids_included = []

    # We exclude a method if it is not in the HIPPIE scoring system or if it is in 
    # the HIPPIE scoring system but with a score lower than 3
    for psimi in methodid2interactions:
        if psimi in psimi2score.keys():
            if psimi2score[psimi] <=3:
                method_ids_excluded.append(psimi)
            elif psimi2score[psimi] > 3:
                method_ids_included.append(psimi)
        else:
            method_ids_excluded.append(psimi)

    method_ex_network_file = os.path.join(main_dir, 'human_edges.eAFF.biana.nov17.method_ex.txt')

    # If using names to exclude, there may be errors because the excluded affinity methods are there in names (but not in ids)
    method_ex_network = main_network.filter_network_by_method(methods_excluded=None, method_ids_excluded=method_ids_excluded, methods_included=None, method_ids_included=None, output_network_file=method_ex_network_file)

    # Analyze the methods network
    print('Method-exclusive-filtered main network')
    print('Excluded methods: {}'.format(method_ids_excluded))
    print('Number of edges: {}'.format(len(method_ex_network.get_edges())))
    print('Number of nodes: {}\n'.format(len(method_ex_network.get_nodes())))

    # method_inc_network_file = os.path.join(main_dir, 'human_edges.eAFF.biana.nov17.method_inc.txt')

    # # If using names to exclude, there may be errors because the excluded affinity methods are there in names (but not in ids)
    # method_inc_network = main_network.filter_network_by_method(methods_excluded=None, method_ids_excluded=None, methods_included=None, method_ids_included=method_ids_included, output_network_file=method_inc_network_file)

    # # Analyze the methods network
    # print('Method-inclusive-filtered main network')
    # print('Included methods: {}'.format(method_ids_included))
    # print('Number of edges: {}'.format(len(method_inc_network.get_edges())))
    # print('Number of nodes: {}\n'.format(len(method_inc_network.get_nodes())))


    # Create a network filtered by number of pubmed IDs
    pmid_network_file = os.path.join(main_dir, 'human_edges.eAFF.biana.nov17.pmid.txt')
    min_num_pubmeds = 2
    pmid_network = main_network.filter_network_by_number_pubmeds(min_num_pubmeds, output_network_file=pmid_network_file)

    # Analyze the pubmed network
    print('Pubmed-filtered main network')
    print('Number of edges: {}'.format(len(pmid_network.get_edges())))
    print('Number of nodes: {}\n'.format(len(pmid_network.get_nodes())))


    # Create a network filtered by database
    databases = {
        'intact' : 'intact [release 2017_04 of 05-apr-2017]',
        'biogrid': 'biogrid [release 3.4.147 (31-mar-2017)]',
        'irefindex' : 'irefindex [14.0 (last edited in 2016-07-23)]',
        'hippie': 'hippie [v2.0 (06/24/2016)]',
        'hprd' : 'hprd [release 2010_04 of 13-apr-2010]',
        'dip': 'dip [release 2017_02 of 05-feb-2017]',
        'gpcr' : 'gpcr [1-apr-2017]'
    }
    intgrid_network_file = os.path.join(main_dir, 'human_edges.eAFF.biana.nov17.intgrid.txt')
    databases_included = [databases['intact'], databases['biogrid']]
    intgrid_network = main_network.filter_network_by_database(databases_included, intgrid_network_file)

    # Analyze the pubmed network

    print('Database-filtered main network')
    print('Number of edges: {}'.format(len(intgrid_network.get_edges())))
    print('Number of nodes: {}\n'.format(len(intgrid_network.get_nodes())))


    # We still have to implement the tissue-specificity to filter by microarray or RNAseq in the case of Human Protein Atlas data!!!

    ##############################################
    ## CREATION OF THE TISSUE-SPECIFIC NETWORKS ##
    ##############################################

    # Define the tissue liver
    tissue_terms_hpa = 'liver'
    tissue_terms_jensen = 'liver'
    liver = NA.Tissue(tissue_terms_hpa, tissue_terms_jensen, jensen_conf=3, hpa_level='low', hpa_rel='approved', pickles_path=pickles_dir)
    #print(liver.tissue_terms_hpa)
    #print(liver.tissue_terms_jensen)
    #print(liver.all_tissue_user_entities)
    #print(liver.all_tissue_proteins)

    # Define the tissue kidney
    tissue_terms_hpa = 'kidney'
    tissue_terms_jensen = 'kidney'
    kidney = NA.Tissue(tissue_terms_hpa, tissue_terms_jensen, jensen_conf=3, hpa_level='low', hpa_rel='approved', pickles_path=pickles_dir)

    # Define the tissue brain
    tissue_terms_hpa = ['cerebral cortex','cerebellum']
    tissue_terms_jensen = 'brain'
    brain = NA.Tissue(tissue_terms_hpa, tissue_terms_jensen, jensen_conf=3, hpa_level='low', hpa_rel='approved', pickles_path=pickles_dir)

    # Define the tissue heart
    tissue_terms_hpa = 'heart muscle'
    tissue_terms_jensen = 'heart'
    heart = NA.Tissue(tissue_terms_hpa, tissue_terms_jensen, jensen_conf=3, hpa_level='low', hpa_rel='approved', pickles_path=pickles_dir)

    # Define the tissue heart
    tissue_terms_hpa = 'pancreas'
    tissue_terms_jensen = 'pancreas'
    pancreas = NA.Tissue(tissue_terms_hpa, tissue_terms_jensen, jensen_conf=3, hpa_level='low', hpa_rel='approved', pickles_path=pickles_dir)


    # Define the housekeeping genes object
    translation_id = 'geneid'
    hk_genes = TS.HouseKeepingGenes(translation_id, pickles_dir)
    all_hk_genes = hk_genes.all_hk_genes
    hpa_hk_genes = hk_genes.hpa_hk_genes
    eisenberg_hk_genes = hk_genes.eisenberg_hk_genes

    # Define the Wang liver network
    wang_liver_file = os.path.join(pickles_dir, 'wang_liver_network.pcl')
    wang_liver_network = cPickle.load(open(wang_liver_file))
    wang_edges_file = os.path.join(main_dir, 'wang_liver_edges.geneid.txt')
    NA.write_network_file_from_networkx_graph(wang_liver_network, wang_edges_file, 'raw')
    wang_network = NA.Network(wang_edges_file, 'geneid', 'raw')


    for network, abbr in [ (main_network, 'complete'), (method_ex_network, 'method_ex'), (pmid_network, 'pmid'), (intgrid_network, 'intgrid') ]:

        results_table.setdefault(abbr, [])

        #### Translate the main network ####
        translation_file = options.translation_file
        translation_id = 'geneid'
        translated_network = os.path.join(main_dir, 'human_edges.eAFF.geneid.nov17.{}.txt'.format(abbr))
        geneid_net = network.translate_network(translation_file, translation_id, network_format, translated_network)

        # Calculate coverage of house-keeping genes in the main network
        num_all_hk_network = float(geneid_net.get_number_of_housekeeping_genes(all_hk_genes))
        num_hpa_hk_network = float(geneid_net.get_number_of_housekeeping_genes(hpa_hk_genes))
        num_eis_hk_network = float(geneid_net.get_number_of_housekeeping_genes(eisenberg_hk_genes))
        num_all_nonhk_network = float(len(set(geneid_net.get_nodes()))) - num_all_hk_network
        per_all_hk_network = num_all_hk_network / float(len(all_hk_genes)) * 100
        per_hpa_hk_network = num_hpa_hk_network / float(len(hpa_hk_genes)) * 100
        per_eis_hk_network = num_eis_hk_network / float(len(eisenberg_hk_genes)) * 100

        print('{} GeneID network'.format(abbr.upper()))
        print('Number of edges: {}'.format(len(geneid_net.get_edges())))
        print('Number of nodes: {}'.format(len(geneid_net.get_nodes())))
        print('Percentage of (all) housekeeping genes in the main network: {:.2f}%\t{:.0f} of {} genes'.format(per_all_hk_network, num_all_hk_network, len(all_hk_genes)))
        print('Percentage of (HPA) housekeeping genes in the main network: {:.2f}%\t{:.0f} of {} genes'.format(per_hpa_hk_network, num_hpa_hk_network, len(hpa_hk_genes)))
        print('Percentage of (Eisenberg) housekeeping genes in the main network: {:.2f}%\t{:.0f} of {} genes\n'.format(per_eis_hk_network, num_eis_hk_network, len(eisenberg_hk_genes)))

        results_table[abbr].append(len(geneid_net.get_edges()))
        results_table[abbr].append(len(geneid_net.get_nodes()))

        for tissue, abbr_tis in [ (liver, 'liver'), (kidney, 'kidney'), (brain, 'brain'), (heart, 'heart'), (pancreas, 'pancreas') ]:

            tissue_dir = os.path.join(tissue_spec_dir,abbr_tis)
            if not os.path.exists(tissue_dir):
                os.makedirs(tissue_dir)

            #### Create tissue-specific network ####
            tissue_network_file = os.path.join(tissue_dir, '{}_edges.biana.nov17.{}.txt'.format(abbr_tis,abbr))
            permission = 0
            tissue_network = network.filter_network_by_tissue(tissue_network_file, tissue, permission)

            #### Translate the tissue-specific network to geneID ####
            translated_network = os.path.join(tissue_dir, '{}_edges.geneid.nov17.{}.txt'.format(abbr_tis,abbr))
            geneid_tissue_network = tissue_network.translate_network(translation_file, translation_id, network_format, translated_network)

            # Get info from tissue-specific network
            hpa_edges=geneid_tissue_network.hpa_edges
            jensen_edges=geneid_tissue_network.jensen_edges
            union=geneid_tissue_network.get_union()
            intersection=geneid_tissue_network.get_intersection()

            # Calculate coverage of housekeeping genes in the tissue-specific network
            num_all_hk_tissue = float(geneid_tissue_network.get_number_of_housekeeping_genes(all_hk_genes))
            num_hpa_hk_tissue = float(geneid_tissue_network.get_number_of_housekeeping_genes(hpa_hk_genes))
            num_eis_hk_tissue = float(geneid_tissue_network.get_number_of_housekeeping_genes(eisenberg_hk_genes))
            num_all_nonhk_tissue = float(len(set(geneid_tissue_network.get_nodes()))) - num_all_hk_tissue
            per_all_hk_tissue = num_all_hk_tissue / num_all_hk_network * 100
            per_hpa_hk_tissue = num_hpa_hk_tissue / num_hpa_hk_network * 100
            per_eis_hk_tissue = num_eis_hk_tissue / num_eis_hk_network * 100

            contingency_table = np.array([[num_all_hk_tissue, num_all_nonhk_tissue], [num_all_hk_network, num_all_nonhk_network]])
            chi2, pval, dof, expected = NA.calculate_contingency_table(contingency_table)

            print('{} {}-specific network'.format(abbr.upper(), abbr_tis.upper()))
            print('Number of edges: {}'.format(len(geneid_tissue_network.get_edges())))
            print('Number of nodes: {}'.format(len(geneid_tissue_network.get_nodes())))
            print('Interactions using only Tissues (Jensen lab): {}'.format(len(jensen_edges)))
            print('Interactions using only Human Protein Atlas: {}'.format(len(hpa_edges)))
            print('Intersection: {}'.format(len(intersection)))
            print('Union: {}\n'.format(len(union)))
            print('Percentage of (all) HK genes in the tissue with respect to main network: {:.2f}%\t{:.0f} of {:.0f} genes'.format(per_all_hk_tissue, num_all_hk_tissue, num_all_hk_network))
            print('Percentage of (HPA) HK genes in the tissue with respect to main network: {:.2f}%\t{:.0f} of {:.0f} genes'.format(per_hpa_hk_tissue, num_hpa_hk_tissue, num_hpa_hk_network))
            print('Percentage of (Eisenberg) HK genes in the tissue with respect to main network: {:.2f}%\t{:.0f} of {:.0f} genes'.format(per_eis_hk_tissue, num_eis_hk_tissue, num_eis_hk_network))
            print('Contingency table: {}'.format(contingency_table))
            print('Contingency table result:\tchi2: {}\tp-value: {}\tdegrees of freedom: {}\n'.format(chi2, pval, dof))

            results_table[abbr].append(len(geneid_tissue_network.get_edges()))
            results_table[abbr].append(len(geneid_tissue_network.get_nodes()))
            results_table[abbr].append(per_all_hk_tissue)
            results_table[abbr].append(per_hpa_hk_tissue)
            results_table[abbr].append(per_eis_hk_tissue)
            results_table[abbr].append(pval)


            if abbr_tis == 'liver':

                # Calculate coverage of Wang liver network genes in the main network
                # and in our tissue-specific netwrok
                wang_nodes_intersection_network = NA.get_nodes_intersection_of_two_networks(geneid_net, wang_network)
                wang_nodes_intersection_liver = NA.get_nodes_intersection_of_two_networks(geneid_tissue_network, wang_network)
                num_wang_network = float(len(wang_nodes_intersection_network))
                num_nonwang_network = float(len(set(geneid_net.get_nodes()))) - num_wang_network
                num_wang_liver = float(len(wang_nodes_intersection_liver))
                num_nonwang_liver = float(len(set(geneid_tissue_network.get_nodes()))) - num_wang_liver
                per_wang_network = num_wang_network / float(len(wang_network.get_nodes())) * 100
                per_wang_tissue = num_wang_liver / num_wang_network * 100

                contingency_table = np.array([[num_wang_liver, num_nonwang_liver], [num_wang_network, num_nonwang_network]])
                chi2, pval, dof, expected = NA.calculate_contingency_table(contingency_table)
                print('Percentage of Wang liver genes in the main network: {:.2f}%\t{:.0f} of {:.0f} genes'.format(per_wang_network, num_wang_network, float(len(wang_network.get_nodes()))))
                print('Percentage of Wang liver genes in the liver-specific network with respect to main network: {:.2f}%\t{:.0f} of {:.0f} genes'.format(per_wang_tissue, num_wang_liver, num_wang_network))
                print('Contingency table: {}'.format(contingency_table))
                print('Contingency table result:\tchi2: {}\tp-value: {}\tdegrees of freedom: {}\n'.format(chi2, pval, dof))

                results_table[abbr].append(per_wang_tissue)
                results_table[abbr].append(pval)

                # Calculate coverage of Wang liver network interactions in the main
                # network and in our tissue-specific netwrok
                wang_intersection_network = NA.get_edges_intersection_of_two_networks(geneid_net, wang_network)
                wang_intersection_liver = NA.get_edges_intersection_of_two_networks(geneid_tissue_network, wang_network)
                num_wang_network = float(len(wang_intersection_network))
                num_nonwang_network = float(len(set(geneid_net.get_edges()))) - num_wang_network
                num_wang_liver = float(len(wang_intersection_liver))
                num_nonwang_liver = float(len(set(geneid_tissue_network.get_edges()))) - num_wang_liver
                per_wang_network = num_wang_network / float(len(wang_network.get_edges())) * 100
                per_wang_tissue = num_wang_liver / num_wang_network * 100

                contingency_table = np.array([[num_wang_liver, num_nonwang_liver], [num_wang_network, num_nonwang_network]])
                chi2, pval, dof, expected = NA.calculate_contingency_table(contingency_table)
                print('Percentage of Wang liver interactions in the main network: {:.2f}%\t{:.0f} of {:.0f} edges'.format(per_wang_network, num_wang_network, float(len(wang_network.get_edges()))))
                print('Percentage of Wang liver interactions in the liver-specific network with respect to main network: {:.2f}%\t{:.0f} of {:.0f} edges'.format(per_wang_tissue, num_wang_liver, num_wang_network))
                print('Contingency table: {}'.format(contingency_table))
                print('Contingency table result:\tchi2: {}\tp-value: {}\tdegrees of freedom: {}\n'.format(chi2, pval, dof))

                results_table[abbr].append(per_wang_tissue)
                results_table[abbr].append(pval)

    output_table_fd = open(output_table_file, 'w')

    for method in results_table:

        output_table_fd.write('{}'.format(method))
        print(method)
        for result in results_table[method]:
            output_table_fd.write('\t{}'.format(result))
            print(result)
        output_table_fd.write('\n')

    output_table_fd.close()
    print(results_table)



    # End marker for time
    end = time.time()
    print('\nTIME OF EXECUTION: {:.3f} seconds or {:.3f} minutes.\n'.format(end - start, (end - start) / 60))

    return



#######################
#######################
# SECONDARY FUNCTIONS #
#######################
#######################

def create_directory(directory):
    """
    Checks if a directory exists and if not, creates it
    """
    try:
        os.stat(directory)
    except:
        os.mkdir(directory)
    return



if  __name__ == "__main__":
    main()

