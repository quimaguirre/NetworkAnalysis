import argparse
import sys, os

from context import NetworkAnalysis
import NetworkAnalysis.network_analysis as NA



def main():

    options = parse_user_arguments()
    create_guildify_files(options)

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
    parser.add_argument('-ts','--tissue_specific',dest='tissue_specific',action = 'store_true',
                        help = """" If the input network is tissue_specific, introduce this argument """)
    parser.add_argument('-ws','--workspace',dest='workspace',action = 'store',default=os.path.join(os.path.join(os.path.dirname(__file__), '..'), 'workspace'),
                        help = """Define the workspace directory where the results will be created. """)

    options=parser.parse_args()

    return options


#################
#################
# MAIN FUNCTION #
#################
#################

def create_guildify_files(options):
    """
    Generates the profiles of the input drug
    """

    # Define the parameters
    output_node_file = os.path.join(options.workspace, 'node_info.txt')
    output_edge_file = os.path.join(options.workspace, 'edge_scores.txt')
    output_edge_file_diamond = os.path.join(options.workspace, 'edge_diamond_file.txt')
    output_multifields_file = os.path.join(options.workspace, 'multi-fields_file.txt')
    type_id = 'biana'
    network_format = options.network_format

    if options.tissue_specific:
        # Define the network instance
        network_instance = NA.TissueSpecificNetwork(options.network_file, type_id, network_format, None, None, node_file=None)
    else:
        # Define the network instance
        network_instance = NA.Network(options.network_file, type_id, network_format, node_file=None)

    # Define the networkX graph
    nx_graph = network_instance.network

    # Create the GUILDify files
    new_network = NA.write_node_and_edge_files_for_guildify_from_networkx_graph(nx_graph, type_id, output_node_file, output_edge_file, output_edge_file_diamond)

    # Create multi-fields file
    if network_format == 'multi-fields':
        # Define the network instance
        NA.write_network_file_from_networkx_graph(new_network, output_multifields_file, 'multi-fields', output_nodes_file=None, tissue_specific=False)

    return


#######################
#######################
# SECONDARY FUNCTIONS #
#######################
#######################


if  __name__ == "__main__":
    main()

