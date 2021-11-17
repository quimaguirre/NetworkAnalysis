import argparse
import sys, os
import cPickle

from context import NetworkAnalysis
import NetworkAnalysis.network_analysis as NA



def main():

    options = parse_user_arguments()
    analyze_network(options)

def parse_user_arguments(*args, **kwds):
    """
    Parses the arguments of the program
    """

    parser = argparse.ArgumentParser(
        description = "Analyze a network",
        epilog      = "@oliva's lab 2017")
    parser.add_argument('-n','--network_file',dest='network_file',action = 'store',
                        help = """ Input file with a protein-protein interaction network in multi-fields format. """)
    parser.add_argument('-f','--network_format',dest='network_format',action = 'store',
                        help = '''Format file of the network file:\tsif, netscore, raw, multi-fields:\n
                                    'sif': <node1>\tscore\t<node2>\n
                                    'netscore': <node1>\t<node2>\t<score>\n
                                    'raw': <node1>\t<node2>\n
                                    'multi-fields' : <node1>\t<node2>\t<sources>\t<method_ids>\t<method_names>\t<pmids>\n''')

    options=parser.parse_args()

    return options


#################
#################
# MAIN FUNCTION #
#################
#################

def analyze_network(options):
    """
    Filters by methods by Hippie score
    """

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

    methodid2interactions = network_instance.get_methodid2interactions(verbose=True) # Get the methods in the network and the number of interactions per method
    network_instance.get_nummethodids2interactions(verbose=True) # Get the number of interactions in which we have a given number of method IDs
    network_instance.get_numpmids2interactions(verbose=True) # Get the number of interactions in which we have a given number of pubmed IDs
    network_instance.get_database2interactions(verbose=True) # Get the database and the number of interactions


#######################
#######################
# SECONDARY FUNCTIONS #
#######################
#######################


if  __name__ == "__main__":
    main()