import argparse
import sys, os
import cPickle

from context import NetworkAnalysis
import NetworkAnalysis.network_analysis as NA



def main():

    options = parse_user_arguments()
    create_perfect_filtered_network(options)

def parse_user_arguments(*args, **kwds):
    """
    Parses the arguments of the program
    """

    parser = argparse.ArgumentParser(
        description = "Create a tissue-specific network",
        epilog      = "@oliva's lab 2017")
    parser.add_argument('-n','--network_file',dest='network_file',action = 'store',
                        help = """ Input file with a protein-protein interaction network in multi-fields format. """)
    parser.add_argument('-ws','--workspace',dest='workspace',action = 'store',default=os.path.join(os.path.join(os.path.dirname(__file__), '..'), 'workspace'),
                        help = """Define the workspace directory where the results will be created. """)

    options=parser.parse_args()

    return options


#################
#################
# MAIN FUNCTION #
#################
#################

def create_perfect_filtered_network(options):
    """
    Filters in various ways
    """

    # Get the program path
    main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    pickles_dir = os.path.join(main_path, 'NetworkAnalysis/pickles')

    # Define the parameters
    network_file = options.network_file
    workspace_dir = options.workspace
    type_id = 'geneID'
    network_format = 'multi-fields'
    score_conf = 5
    score_semiconf = 1
    min_num_methods = 3
    min_num_pubmeds = 5

    # Load HIPPIE scores and print them
    psimi2score_file = os.path.join(pickles_dir, 'psimi2score.pcl')
    psimi2score = cPickle.load(open(psimi2score_file))
    psimi2method_file = os.path.join(pickles_dir,'psimi2method.pcl')
    psimi2method = cPickle.load(open(psimi2method_file))
    print_str = 0
    print('\nHIPPIE SCORES:\n')
    for psimi, score in sorted(psimi2score.iteritems(), key=lambda (x, y): y, reverse=True):
        if score >= score_conf and print_str == 0:
            print('CONFIDENCE METHODS:')
            print_str+=1
        elif score < score_conf and score >= score_semiconf and print_str == 1:
            print('SEMICONFIDENCE METHODS:')
            print_str+=1
        elif score < score_semiconf and print_str == 2:
            print('LOW CONFIDENCE METHODS')
            print_str+=1
        if psimi in psimi2method:
            name = psimi2method[psimi]
        else:
            name = '-'
        print('Method ID: {}\tName: {}\tScore: {}'.format(psimi, name, score))
    print('')

    # Define the network instance
    network_instance = NA.Network(network_file, type_id, network_format, node_file=None)

    # Analyze the main network
    print('Input network')
    print('Number of edges: {}'.format(len(network_instance.get_edges())))
    print('Number of nodes: {}\n'.format(len(network_instance.get_nodes())))
    methodid2interactions = network_instance.get_methodid2interactions() # Get the methods in the network and the number of interactions per method

    # Classify method IDs by their confidence using the HIPPIE scores
    # If using names, there may be errors because the excluded affinity methods are there in names (but not in ids)
    method_ids_confident = []
    method_ids_semiconfident = []
    method_ids_excluded = []

    for psimi in methodid2interactions:
        if psimi in psimi2score.keys():
            if psimi2score[psimi] >= score_conf:
                method_ids_confident.append(psimi)
            if psimi2score[psimi] >= score_semiconf:
                method_ids_semiconfident.append(psimi)
        else:
            # Exclude the methods for which we do not have information
            method_ids_excluded.append(psimi)

    # Create a confidence network
    confidence_nework_file = os.path.join(workspace_dir, 'network.confidence.txt')
    confidence_network = network_instance.filter_network_by_method(methods_excluded=None, method_ids_excluded=None, methods_included=None, method_ids_included=method_ids_confident, output_network_file=confidence_nework_file)

    print('Confidence network')
    print('Included methods: {}'.format(method_ids_confident))
    print('Number of edges: {}'.format(len(confidence_network.get_edges())))
    print('Number of nodes: {}\n'.format(len(confidence_network.get_nodes())))

    # Create a network of interactions with at least 2 methods
    nmet_network_file = os.path.join(workspace_dir, 'network.{}met.txt'.format(str(min_num_methods)))
    nmet_network = network_instance.filter_network_by_number_methods(min_num_methods, output_network_file=nmet_network_file, output_nodes_file=None)
    print('>= {} methods network'.format(str(min_num_methods)))
    print('Number of edges: {}'.format(len(nmet_network.get_edges())))
    print('Number of nodes: {}\n'.format(len(nmet_network.get_nodes())))

    # Filter the >= 2 methods network by semiconfidence methods
    semiconfidence_nework_file = os.path.join(workspace_dir, 'network.semiconfidence.txt')
    semiconfidence_network = nmet_network.filter_network_by_method(methods_excluded=None, method_ids_excluded=None, methods_included=None, method_ids_included=method_ids_semiconfident, output_network_file=semiconfidence_nework_file)

    print('Semiconfidence network')
    print('Included methods: {}'.format(method_ids_semiconfident))
    print('Number of edges: {}'.format(len(semiconfidence_network.get_edges())))
    print('Number of nodes: {}\n'.format(len(semiconfidence_network.get_nodes())))

    # Create a network filtered by number of pubmed IDs
    pmid_network_file = os.path.join(workspace_dir, 'network.{}pmid.txt'.format(str(min_num_pubmeds)))
    pmid_network = network_instance.filter_network_by_number_pubmeds(min_num_pubmeds, output_network_file=pmid_network_file, output_nodes_file=None)

    print('>= {} pubmeds network'.format(str(min_num_pubmeds)))
    print('Number of edges: {}'.format(len(pmid_network.get_edges())))
    print('Number of nodes: {}\n'.format(len(pmid_network.get_nodes())))

    # Create a union of semiconfidence and pubmed networks
    pmidormet_network_file = os.path.join(workspace_dir, 'network.semiconfidence.{}pmid.txt'.format(str(min_num_pubmeds)))
    pmidormet_network = NA.get_union_network_from_two_networks(semiconfidence_network, pmid_network, pmidormet_network_file)

    print('Semiconfidence-Pubmed network')
    print('Number of edges: {}'.format(len(pmidormet_network.get_edges())))
    print('Number of nodes: {}\n'.format(len(pmidormet_network.get_nodes())))

    return


#######################
#######################
# SECONDARY FUNCTIONS #
#######################
#######################


if  __name__ == "__main__":
    main()