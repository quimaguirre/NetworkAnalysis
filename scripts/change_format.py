import argparse
import time
import sys, os, re

from context import NetworkAnalysis
import NetworkAnalysis.network_analysis as NA




def main():

    options = parse_user_arguments()
    filter_network(options)


def parse_user_arguments(*args, **kwds):
    """
    Parses the arguments of the program
    """

    parser = argparse.ArgumentParser(
        description = "Filter network",
        epilog      = "@oliva's lab 2017")
    parser.add_argument('-i','--input_file',dest='input_file',action = 'store',
                        help = """ Input file with a protein-protein interaction network. """)
    parser.add_argument('-o','--output_file',dest='output_file',action = 'store',
                        help = """ Output file with a protein-protein interaction network. """)
    parser.add_argument('-iformat','--input_format',dest='input_format',action = 'store',default='multi-fields',
                        help = '''Format file of the edge file:\tsif (default), netscore, raw, multi-fields:\n
                                    'sif': <node1>\tscore\t<node2>\n
                                    'netscore': <node1>\t<node2>\t<score>\n
                                    'raw': <node1>\t<node2>\n
                                    'multi-fields' : <node1>\t<node2>\t<sources>\t<method_ids>\t<method_names>\t<pmids>\n''')
    parser.add_argument('-oformat','--output_format',dest='output_format',action = 'store',default='sif',
                        help = '''Format file of the edge file:\tsif (default), netscore, raw, multi-fields:\n
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

def filter_network(options):
    """
    Generates the profiles of the input drug
    """

    # Start marker for time measure
    start = time.time()


    #-----------------------#
    #   PARSE THE NETWORK   #
    #-----------------------#

    input_file = options.input_file
    output_file = options.output_file
    type_id = 'biana'
    input_format = options.input_format
    output_format = options.output_format
    network_instance = NA.Network(input_file, type_id, input_format)

    # Analyze the main network
    print('Network:')
    print('Number of edges: {}'.format(len(network_instance.get_edges())))
    print('Number of nodes: {}\n'.format(len(network_instance.get_nodes())))


    # Write the output network
    NA.write_network_file_from_networkx_graph(network_instance.network, output_file, output_format, output_nodes_file=None, tissue_specific=False)


    # End marker for time
    end = time.time()
    print('\nTIME OF EXECUTION: {:.3f} seconds or {:.3f} minutes.\n'.format(end - start, (end - start) / 60))


    return




#######################
#######################
# SECONDARY FUNCTIONS #
#######################
#######################

def fileExist(file):
    """
    Checks if a file exists AND is a file
    """
    return os.path.exists(file) and os.path.isfile(file)


if  __name__ == "__main__":
    main()
