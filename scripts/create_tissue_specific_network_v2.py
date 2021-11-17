import NetworkAnalysis.network_analysis as NA
import NetworkAnalysis.tissue_specificity as TS
import cPickle
import numpy as np
import sys, os
import pandas as pd

###################################
## CREATION OF THE MAIN NETWORKS ##
###################################

# Define the paths
main_dir = '/home/quim/project/tissue_specificity/main' # Path to store the main networks
tissue_spec_dir = '/home/quim/project/tissue_specificity/tissue_spec' # Path to store the tissue specific networks
scripts_dir = '/home/quim/project/tissue_specificity/scripts/pickles' # Path to execute the pickles

# Define the results variable
results_table = {}
output_table_file = '/home/quim/project/tissue_specificity/' + 'tissue_specific_results.tsv'


# Define the main network

network_file = os.path.join(main_dir, 'human_edges.eAFF.biana.010617.txt')
node_file = os.path.join(main_dir, 'human_nodes.eAFF.biana.010617.txt')
type_id = 'biana'
network_format = 'multi-fields'

main_network = NA.Network(network_file, node_file, type_id, network_format)

# Analyze the main network

print('BIANA main network')
print('Number of edges: {}'.format(len(main_network.get_edges())))
print('Number of nodes: {}\n'.format(len(main_network.get_nodes())))
methodid2interactions = main_network.get_methodid2interactions() # Get the methods in the network and the number of interactions per method
main_network.get_numpmids2interactions() # Get the number of interactions in which we have a given number of pubmed IDs
main_network.get_database2interactions() # Get the database and the number of interactions


# Create a network filtered by methods

psimi2score_file = os.path.join(scripts_dir, 'psimi2score.pcl')
psimi2score = cPickle.load(open(psimi2score_file))
method_ids_excluded = []

# We exclude a method if it is not in the HIPPIE scoring system or if it is in 
# the HIPPIE scoring system but with a score lower than 3
for psimi in methodid2interactions:
    if psimi in psimi2score.keys():
        if psimi2score[psimi] <=3:
            method_ids_excluded.append(psimi)
    else:
        method_ids_excluded.append(psimi)

method_network_file = os.path.join(main_dir, 'human_edges.eAFF.biana.method.txt')
method_nodes_file = os.path.join(main_dir, 'human_nodes.eAFF.biana.method.txt')

# If using names to exclude, there may be errors because the excluded affinity methods are there in names (but not in ids)
method_network = main_network.filter_network_by_method(methods_excluded=None, method_ids_excluded=method_ids_excluded, methods_included=None, method_ids_included=None, output_network_file=method_network_file, output_nodes_file=method_nodes_file)

# Analyze the methods network

print('Methods-filtered main network')
print('Excluded methods: {}'.format(method_ids_excluded))
print('Number of edges: {}'.format(len(method_network.get_edges())))
print('Number of nodes: {}\n'.format(len(method_network.get_nodes())))


# Create a network filtered by number of pubmed IDs

pmid_network_file = os.path.join(main_dir, 'human_edges.eAFF.biana.pmid.txt')
pmid_nodes_file = os.path.join(main_dir, 'human_nodes.eAFF.biana.pmid.txt')
min_num_pubmeds = 2
pmid_network = main_network.filter_network_by_number_pubmeds(min_num_pubmeds, output_network_file=pmid_network_file, output_nodes_file=pmid_nodes_file)

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

intgrid_network_file = os.path.join(main_dir, 'human_edges.eAFF.biana.intgrid.txt')
intgrid_nodes_file = os.path.join(main_dir, 'human_nodes.eAFF.biana.intgrid.txt')
databases_included = [databases['intact'], databases['biogrid']]
intgrid_network = main_network.filter_network_by_database(databases_included, intgrid_network_file, intgrid_nodes_file)

# Analyze the pubmed network

print('Database-filtered main network')
print('Number of edges: {}'.format(len(intgrid_network.get_edges())))
print('Number of nodes: {}\n'.format(len(intgrid_network.get_nodes())))


##############################################
## CREATION OF THE TISSUE-SPECIFIC NETWORKS ##
##############################################

# Define the tissue liver
tissue_terms_hpa = 'liver'
tissue_terms_jensen = 'liver'
liver = NA.Tissue(tissue_terms_hpa, tissue_terms_jensen, jensen_conf=3, hpa_level='low', hpa_rel='approved')
#print(liver.tissue_terms_hpa)
#print(liver.tissue_terms_jensen)
#print(liver.all_tissue_user_entities)
#print(liver.all_tissue_proteins)

# Define the tissue kidney
tissue_terms_hpa = 'kidney'
tissue_terms_jensen = 'kidney'
kidney = NA.Tissue(tissue_terms_hpa, tissue_terms_jensen, jensen_conf=3, hpa_level='low', hpa_rel='approved')

# Define the tissue brain
tissue_terms_hpa = ['cerebral cortex','cerebellum']
tissue_terms_jensen = 'brain'
brain = NA.Tissue(tissue_terms_hpa, tissue_terms_jensen, jensen_conf=3, hpa_level='low', hpa_rel='approved')

# Define the tissue heart
tissue_terms_hpa = 'heart muscle'
tissue_terms_jensen = 'heart'
heart = NA.Tissue(tissue_terms_hpa, tissue_terms_jensen, jensen_conf=3, hpa_level='low', hpa_rel='approved')

# Define the tissue heart
tissue_terms_hpa = 'pancreas'
tissue_terms_jensen = 'pancreas'
pancreas = NA.Tissue(tissue_terms_hpa, tissue_terms_jensen, jensen_conf=3, hpa_level='low', hpa_rel='approved')


# Define the housekeeping genes object
translation_id = 'geneid'
hk_genes = TS.HouseKeepingGenes(translation_id)
all_hk_genes = hk_genes.all_hk_genes
hpa_hk_genes = hk_genes.hpa_hk_genes
eisenberg_hk_genes = hk_genes.eisenberg_hk_genes

# Define the Wang liver network
wang_liver_file = os.path.join(scripts_dir, 'wang_liver_network.pcl')
wang_liver_network = cPickle.load(open(wang_liver_file))
wang_edges_file = os.path.join(main_dir, 'wang_liver_edges.geneid.txt')
wang_nodes_file = os.path.join(main_dir, 'wang_liver_nodes.geneid.txt')
NA.write_network_file_from_networkx_graph(wang_liver_network, wang_edges_file, wang_nodes_file, 'raw')
wang_network = NA.Network(wang_edges_file, wang_nodes_file, 'geneid', 'raw')


for network, abbr in [ (main_network, '010617'), (method_network, 'method'), (pmid_network, 'pmid'), (intgrid_network, 'intgrid') ]:

    results_table.setdefault(abbr, [])

    #### Translate the main network ####
    translation_file = os.path.join(main_dir, 'human_network.eAFF.010617.geneid.trans')
    translation_id = 'geneid'
    translated_network = os.path.join(main_dir, 'human_edges.eAFF.geneid.{}.txt'.format(abbr))
    translated_nodes = os.path.join(main_dir, 'human_nodes.eAFF.geneid.{}.txt'.format(abbr))
    geneid_net = network.translate_network(translation_file, translation_id, network_format, translated_network, translated_nodes)

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
        tissue_network_file = os.path.join(tissue_dir, '{}_edges.biana.{}.txt'.format(abbr_tis,abbr))
        tissue_nodes_file = os.path.join(tissue_dir, '{}_nodes.biana.{}.txt'.format(abbr_tis,abbr))
        permission = 0
        tissue_network = network.filter_network_by_tissue(tissue_network_file, tissue_nodes_file, tissue, permission)

        #### Translate the tissue-specific network to geneID ####
        translated_network = os.path.join(tissue_dir, '{}_edges.geneid.{}.txt'.format(abbr_tis,abbr))
        translated_nodes = os.path.join(tissue_dir, '{}_nodes.geneid.{}.txt'.format(abbr_tis,abbr))
        geneid_tissue_network = tissue_network.translate_network(translation_file, translation_id, network_format, translated_network, translated_nodes)

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
