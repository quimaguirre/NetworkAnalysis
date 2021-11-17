import argparse
import sys, os
import cPickle

from context import NetworkAnalysis
import NetworkAnalysis.network_analysis as NA



def main():

    options = parse_user_arguments()
    filter_network_by_hippie_methods(options)

def parse_user_arguments(*args, **kwds):
    """
    Parses the arguments of the program
    """

    parser = argparse.ArgumentParser(
        description = "Create a tissue-specific network",
        epilog      = "@oliva's lab 2017")
    parser.add_argument('-n','--network_file',dest='network_file',action = 'store',
                        help = """ Input file with a protein-protein interaction network in multi-fields format. """)
    parser.add_argument('-o','--output_file',dest='output_file',action = 'store',
                        help = """ Output file of the filtered protein-protein interaction network in multi-fields format. """)
    parser.add_argument('-s','--score',dest='score',action = 'store', default = 3,
                        help = """ Input the cut-off score from which you want to exclude/include methods """)
    parser.add_argument('-c','--criteria',dest='criteria',action = 'store', default = 'exclude',
                        help = """ Introduce the criteria of the filtering: 'include' or 'exclude' the interactions with the methods """)

    options=parser.parse_args()

    return options


#################
#################
# MAIN FUNCTION #
#################
#################

def filter_network_by_hippie_methods(options):
    """
    Filters by methods by Hippie score
    """

    # Get the program path
    main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    pickles_dir = os.path.join(main_path, 'NetworkAnalysis/pickles')

    # Define the parameters
    network_file = options.network_file
    output_file = options.output_file
    score = float(options.score)
    criteria = options.criteria
    type_id = 'biana'
    network_format = 'multi-fields'

    # Load HIPPIE scores and print them
    psimi2score_file = os.path.join(pickles_dir, 'psimi2score.pcl')
    psimi2score = cPickle.load(open(psimi2score_file))
    psimi2method_file = os.path.join(pickles_dir,'psimi2method.pcl')
    psimi2method = cPickle.load(open(psimi2method_file))
    print_str = 0
    print('\nHIPPIE SCORES:\n')
    for psimi, hippie_score in sorted(psimi2score.iteritems(), key=lambda (x, y): y, reverse=True):
        if hippie_score >= score and print_str == 0:
            print('METHODS HIGHER THAN {}:'.format(score))
            print_str+=1
        if hippie_score < score and print_str == 1:
            print('METHODS LOWER THAN {}:'.format(score))
            print_str+=1
        if psimi in psimi2method:
            name = psimi2method[psimi]
        else:
            name = '-'
        print('Method ID: {}\tName: {}\tScore: {}'.format(psimi, name, hippie_score))
    print('')

    # Define the network instance
    network_instance = NA.Network(network_file, type_id, network_format, node_file=None)

    # Analyze the main network
    print('Input network')
    print('Number of edges: {}'.format(len(network_instance.get_edges())))
    print('Number of nodes: {}\n'.format(len(network_instance.get_nodes())))
    methodid2interactions=network_instance.get_methodid2interactions(verbose=True) # Get the methods in the network and the number of interactions per method

    # Create a network filtered by methods
    method_ids_excluded = []
    method_ids_included = []

    # We exclude a method if it is not in the HIPPIE scoring system or if it is in 
    # the HIPPIE scoring system but with a score lower than 3 when criterium 'exclude'
    for psimi in methodid2interactions:
        if psimi in psimi2score.keys():
            if criteria == 'exclude':
                if psimi2score[psimi] <= score:
                    method_ids_excluded.append(psimi)
            elif criteria == 'include':
                if psimi2score[psimi] >= score:
                    method_ids_included.append(psimi)
        else:
            # Exclude the methods for which we do not have information
            method_ids_excluded.append(psimi)


    # If using names to exclude, there may be errors because the excluded affinity methods are there in names (but not in ids)
    method_network = network_instance.filter_network_by_method(methods_excluded=None, method_ids_excluded=method_ids_excluded, methods_included=None, method_ids_included=method_ids_included, output_network_file=output_file)

    # Analyze the methods network
    print('Method-({})-filtered main network'.format(criteria))
    print('Excluded methods: {}'.format(method_ids_excluded))
    print('Included methods: {}'.format(method_ids_included))
    print('Number of edges: {}'.format(len(method_network.get_edges())))
    print('Number of nodes: {}\n'.format(len(method_network.get_nodes())))
    method_network.get_methodid2interactions(verbose=True) # Get the methods in the network and the number of interactions per method

    workspace_dir = os.path.join(main_path, 'workspace')
    output_edge_file = os.path.join(workspace_dir, 'edge_scores.txt')
    output_edge_file_diamond = os.path.join(workspace_dir, 'edge_diamond_file.txt')
    NA.write_network_file_from_networkx_graph(method_network.network, output_edge_file, 'sif', output_nodes_file=None, tissue_specific=False)
    NA.write_network_file_from_networkx_graph(method_network.network, output_edge_file_diamond, 'raw', output_nodes_file=None, tissue_specific=False)

    return


#######################
#######################
# SECONDARY FUNCTIONS #
#######################
#######################


if  __name__ == "__main__":
    main()