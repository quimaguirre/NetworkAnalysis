import argparse
import time
import sys, os, re

from context import NetworkAnalysis
import NetworkAnalysis.network_analysis as NA
import NetworkAnalysis.validation as VA


def main():
    """
    Validation of a PPI network by comparing to other networks of reference.
    """

    # Start marker for time measure
    start = time.time()

    #------------#
    #   INPUTS   #
    #------------#

    # Our input network 1
    input_network_file1 = '/home/quim/data/networks/human_guildify_jan18/multifields/human_network_guildify_jan18.geneID.multifields'
    type_id = 'geneID'
    network_format = 'multi-fields'
    name_input_network1 = 'BIANA MAIN'

    # Our input network 2
    input_network_file2 = '/home/quim/data/networks/human_guildify_jan18/multifields/network.method.txt'
    network_format2 = 'multi-fields'
    name_input_network2 = 'BIANA TAP FILTERED'

    # Our input network 3
    input_network_file3 = '/home/quim/data/networks/human_guildify_jan18/multifields/network.2db.txt'
    network_format3 = 'multi-fields'
    name_input_network3 = 'BIANA 2 >= NUMBER DATABASES'

    # HIPPIE network
    hippie_file = '/home/quim/Databases/hippie/HIPPIE-current.mitab.txt'
    hippie_type_id = 'geneID' # It can be geneID or UniprotEntry
    output_hippie_file = '/home/quim/data/networks/HIPPIE/HIPPIE.{}.multifields'.format(hippie_type_id)
    output_hippie_newID_file = '/home/quim/data/networks/HIPPIE/HIPPIE.{}.multifields'.format(type_id)
    hippie_network_format = 'multi-fields'

    # ConsensusPathDB network
    ConsensusPathDB_file = '/home/quim/data/networks/ConsensusPathDB/ConsensusPathDB_human_PPI'
    output_ConsensusPath_file = '/home/quim/data/networks/ConsensusPathDB/ConsensusPathDB_human_PPI.multifields'
    output_ConsensusPath_newID_file = '/home/quim/data/networks/ConsensusPathDB/ConsensusPathDB_human_PPI.{}.multifields'.format(type_id)
    consensus_network_format = 'multi-fields'

    # I2D network
    I2D_file = '/home/quim/data/networks/I2D/i2d.2_9.Public.HUMAN.tab'
    output_I2D_file = '/home/quim/data/networks/I2D/i2d.2_9.Public.HUMAN.multifields'
    output_I2D_newID_file = '/home/quim/data/networks/I2D/i2d.2_9.Public.HUMAN.{}.multifields'.format(type_id)
    I2D_network_format = 'multi-fields'


    #------------------------#
    #   DEFINE OUR NETWORK   #
    #------------------------#

    # Define the input network 1
    network1 = NA.Network(input_network_file1, type_id, network_format)

    print('{} network'.format(name_input_network1))
    print('Number of edges: {}'.format(len(network1.get_edges())))
    print('Number of nodes: {}\n'.format(len(network1.get_nodes())))

    # Define the input network 2
    network2 = NA.Network(input_network_file2, type_id, network_format2)

    print('{} network'.format(name_input_network2))
    print('Number of edges: {}'.format(len(network2.get_edges())))
    print('Number of nodes: {}\n'.format(len(network2.get_nodes())))

    # Define the input network 3
    network3 = NA.Network(input_network_file3, type_id, network_format3)

    print('{} network'.format(name_input_network3))
    print('Number of edges: {}'.format(len(network3.get_edges())))
    print('Number of nodes: {}\n'.format(len(network3.get_nodes())))


    #---------------------------#
    #   DEFINE HIPPIE NETWORK   #
    #---------------------------#

    if not fileExist(output_hippie_file):
        hippie_instance = VA.HippieParser(hippie_file)
        hippie_instance.parse()
        hippie_network = hippie_instance.write_network_file(output_hippie_file, hippie_network_format, hippie_type_id)
    else:
        hippie_network = NA.Network(output_hippie_file, hippie_type_id, hippie_network_format)

    # Translate HIPPIE to 'type_id'
    if type_id.lower() != hippie_type_id.lower():
        if not fileExist(output_hippie_newID_file):
            hippie_network = VA.translate_network_from_BIANA(hippie_network, hippie_type_id, type_id, output_hippie_newID_file)
        else:
            hippie_network = NA.Network(output_hippie_newID_file, type_id, hippie_network_format)

    print('HIPPIE network')
    print('Number of edges: {}'.format(len(hippie_network.get_edges())))
    print('Number of nodes: {}\n'.format(len(hippie_network.get_nodes())))


    #------------------------------------#
    #   DEFINE CONSENSUSPATHDB NETWORK   #
    #------------------------------------#

    if not fileExist(output_ConsensusPath_file):
        consensus_instance = VA.ConsensusPathDBParser(ConsensusPathDB_file)
        consensus_network_uniprot = consensus_instance.parse(output_ConsensusPath_file, consensus_network_format)
    else:
        consensus_network_uniprot = NA.Network(output_ConsensusPath_file, 'uniprotentry', consensus_network_format)

    # Translate ConsensusPathDB to 'type_id'
    if type_id.lower() != 'uniprotentry':
        if not fileExist(output_ConsensusPath_newID_file):
            consensus_network = VA.translate_network_from_BIANA(consensus_network_uniprot, 'uniprotentry', type_id, output_ConsensusPath_newID_file)
        else:
            consensus_network = NA.Network(output_ConsensusPath_newID_file, type_id, consensus_network_format)
    else:
        consensus_network = consensus_network_uniprot

    print('ConsensusPath (uniprotentry) network')
    print('Number of edges: {}'.format(len(consensus_network_uniprot.get_edges())))
    print('Number of nodes: {}\n'.format(len(consensus_network_uniprot.get_nodes())))

    print('ConsensusPath network')
    print('Number of edges: {}'.format(len(consensus_network.get_edges())))
    print('Number of nodes: {}\n'.format(len(consensus_network.get_nodes())))


    #------------------------#
    #   DEFINE I2D NETWORK   #
    #------------------------#

    if not fileExist(output_I2D_file):
        I2D_instance = VA.I2DParser(I2D_file)
        I2D_network_uniprot = I2D_instance.parse(output_I2D_file, I2D_network_format)
    else:
        I2D_network_uniprot = NA.Network(output_I2D_file, 'uniprotaccession', I2D_network_format)

    # Translate I2D to 'type_id'
    if type_id.lower() != 'uniprotaccession':
        if not fileExist(output_I2D_newID_file):
            I2D_network = VA.translate_network_from_BIANA(I2D_network_uniprot, 'uniprotaccession', type_id, output_I2D_newID_file)
        else:
            I2D_network = NA.Network(output_I2D_newID_file, type_id, I2D_network_format)
    else:
        I2D_network = I2D_network_uniprot

    print('I2D (uniprotaccession) network')
    print('Number of edges: {}'.format(len(I2D_network_uniprot.get_edges())))
    print('Number of nodes: {}\n'.format(len(I2D_network_uniprot.get_nodes())))

    print('I2D network')
    print('Number of edges: {}'.format(len(I2D_network.get_edges())))
    print('Number of nodes: {}\n'.format(len(I2D_network.get_nodes())))


    #----------------------------------#
    #   CHECK OVERLAP OF NODES/EDGES   #
    #----------------------------------#

    print_summary_overlap(network1, hippie_network, name_input_network1, 'HIPPIE')
    print_summary_overlap(network1, consensus_network, name_input_network1, 'CONSENSUSPATHDB')
    print_summary_overlap(network1, I2D_network, name_input_network1, 'I2D')

    print_summary_overlap(network2, hippie_network, name_input_network2, 'HIPPIE')
    print_summary_overlap(network2, consensus_network, name_input_network2, 'CONSENSUSPATHDB')
    print_summary_overlap(network2, I2D_network, name_input_network2, 'I2D')

    print_summary_overlap(network3, hippie_network, name_input_network3, 'HIPPIE')
    print_summary_overlap(network3, consensus_network, name_input_network3, 'CONSENSUSPATHDB')
    print_summary_overlap(network3, I2D_network, name_input_network3, 'I2D')

    print_summary_overlap(hippie_network, consensus_network, 'HIPPIE', 'CONSENSUSPATHDB')
    print_summary_overlap(hippie_network, I2D_network, 'HIPPIE', 'I2D')
    print_summary_overlap(consensus_network, I2D_network, 'CONSENSUSPATHDB', 'I2D')


    # End marker for time
    end = time.time()
    print('\nTIME OF EXECUTION: {:.3f} seconds or {:.3f} minutes.\n'.format(end - start, (end - start) / 60))



#######################
#######################
# SECONDARY FUNCTIONS #
#######################
#######################

def fileExist(file):
    """
    Checks if a file exists AND is a file.
    """
    return os.path.exists(file) and os.path.isfile(file)

def print_summary_overlap(network1, network2, name_network1, name_network2):
    """
    Prints the information of the overlap of two networks.
    """
    common_edges = NA.get_edges_intersection_of_two_networks(network1, network2)
    common_nodes = NA.get_nodes_intersection_of_two_networks(network1, network2)
    percent_net1 = float(len(common_edges)) / float(len(network1.get_edges())) * 100
    percent_net2 = float(len(common_edges)) / float(len(network2.get_edges())) * 100
    print('Common edges between {} and {}: {}'.format(name_network1, name_network2, len(common_edges)))
    print('Common nodes between {} and {}: {}'.format(name_network1, name_network2, len(common_nodes)))
    print('{:.2f}% of {} covered and {:.2f}% of {} covered\n'.format(percent_net1, name_network1, percent_net2, name_network2))
    return


if  __name__ == "__main__":
    main()

