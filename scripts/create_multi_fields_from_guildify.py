import argparse
import sys, os

from context import NetworkAnalysis
import NetworkAnalysis.network_analysis as NA



def main():

    options = parse_user_arguments()
    create_multifields(options)

def parse_user_arguments(*args, **kwds):
    """
    Parses the arguments of the program
    """

    parser = argparse.ArgumentParser(
        description = "Create a tissue-specific network",
        epilog      = "@oliva's lab 2017")
    parser.add_argument('-i','--input_file',dest='input_file',action = 'store',
                        help = """ File of the input network. """)
    parser.add_argument('-nod','--node_info',dest='node_info_file',action = 'store',
                        help = """ Node info file. """)
    parser.add_argument('-met','--method_file',dest='method_file',action = 'store',
                        help = """ File of the input network with method ids. """)
    parser.add_argument('-pub','--pubmed_file',dest='pubmed_file',action = 'store',
                        help = """ File of the input network with pubmed ids. """)
    parser.add_argument('-src','--source_file',dest='source_file',action = 'store',
                        help = """ File of the input network with sources. """)
    parser.add_argument('-o','--output_file',dest='output_file',action = 'store',
                        help = """ File of the output network in multi-fields format. """)
    parser.add_argument('-ot','--output_trans',dest='output_trans_file',action = 'store',
                        help = """ Output translation file """)

    options=parser.parse_args()

    return options


#################
#################
# MAIN FUNCTION #
#################
#################

def create_multifields(options):
    """
    Creates a multi-fields network file from GUILDify output files
    """

    # Define the parameters
    input_file = options.input_file
    node_info_file = options.node_info_file
    method_file = options.method_file
    pubmed_file = options.pubmed_file
    source_file = options.source_file
    output_file = options.output_file
    output_trans_file = options.output_trans_file

    # Parse network, methods, pubmeds and sources files
    guildify_network = parse_guildify_network(input_file)
    interaction_to_method = parse_guildify_information_file(method_file)
    interaction_to_pubmed = parse_guildify_information_file(pubmed_file)
    interaction_to_source = parse_guildify_information_file(source_file, source=True)

    # Write multi-fields
    with open(output_file, 'w') as output_fd:
        for edge in guildify_network:
            node1, node2 = edge
            if edge not in interaction_to_method or edge not in interaction_to_source or edge not in interaction_to_pubmed:
                continue
            method_ids = ';'.join(interaction_to_method[edge])
            sources = ';'.join(interaction_to_source[edge])
            pubmeds = ';'.join(interaction_to_pubmed[edge])
            method_names = '-'
            output_fd.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format( node1, node2, sources, method_ids, method_names, pubmeds ))

    # Set memory free
    guildify_network = set()
    interaction_to_method = {}
    interaction_to_pubmed = {}
    interaction_to_source = {}

    # Parse node info file
    node_to_geneID, node_to_uniprotacc, node_to_genesymbol = parse_node_info_file(node_info_file)

    # Write translation files
    write_translation_file(node_to_geneID, output_trans_file+'.geneID')
    write_translation_file(node_to_uniprotacc, output_trans_file+'.uniprotacc')
    write_translation_file(node_to_genesymbol, output_trans_file+'.genesymbol')

    return


#######################
#######################
# SECONDARY FUNCTIONS #
#######################
#######################

def parse_guildify_network(input_file):
    """
    Parses GUILDify network
    """
    guildify_network = set()
    with open(input_file, 'r') as input_fd:
        for line in input_fd:
            fields = line.strip().split()
            # Skip first lines of network with only nodes
            if len(fields) < 3:
                continue
            node1, _, node2 = line.strip().split()
            edge = frozenset([node1, node2])
            guildify_network.add(edge)
    return guildify_network

def parse_guildify_information_file(input_file, source=False):
    """
    Parses GUILDify files of information (pubmed, method_id, source (if source=True))
    """
    interaction_to_content = {}
    with open(input_file, 'r') as input_fd:
        first_line = input_fd.readline()
        for line in input_fd:
            node1, _, node2, _, content = line.strip().split()
            if source:
                content = content.split('(')[0]
            edge = frozenset([node1, node2])
            interaction_to_content.setdefault(edge, set())
            interaction_to_content[edge].add(content)
    return interaction_to_content

def parse_node_info_file(input_file):
    """
    Parses the node info file and returns the translation of the nodes to 
    GeneID, UniprotAccession and GeneSymbol
    """
    node_to_geneID = {}
    node_to_uniprotacc = {}
    node_to_genesymbol = {}
    with open(input_file, 'r') as input_fd:
        first_line = input_fd.readline()
        for line in input_fd:
            user_entity, uniprotaccs, gene_symbols, name, geneIDs, eq_entries = line.strip().split('\t')
            if geneIDs != '-':
                geneIDs = geneIDs.split('; ')
                for geneID in geneIDs:
                    node_to_geneID.setdefault(user_entity, set())
                    node_to_geneID[user_entity].add(geneID)
            else:
                node_to_geneID.setdefault(user_entity, set())
                node_to_geneID[user_entity].add('')
            if uniprotaccs != '-':
                uniprotaccs = uniprotaccs.split('; ')
                for uacc in uniprotaccs:
                    node_to_uniprotacc.setdefault(user_entity, set())
                    node_to_uniprotacc[user_entity].add(uacc)
            else:
                node_to_uniprotacc.setdefault(user_entity, set())
                node_to_uniprotacc[user_entity].add('')
            if gene_symbols != '-':
                gene_symbols = gene_symbols.split('; ')
                for gene_symbol in gene_symbols:
                    node_to_genesymbol.setdefault(user_entity, set())
                    node_to_genesymbol[user_entity].add(gene_symbol)
            else:
                node_to_genesymbol.setdefault(user_entity, set())
                node_to_genesymbol[user_entity].add('')
    return node_to_geneID, node_to_uniprotacc, node_to_genesymbol

def write_translation_file(node_to_trans, output_file):
    """
    Writes the translation files
    """
    with open(output_file, 'w') as output_fd:
        for node in node_to_trans:
            translated = []
            for translation in node_to_trans[node]:
                translated.append("'{}'".format(translation))
            output_fd.write("{0}\t{1}\n".format(node,','.join(translated)))
    return


if  __name__ == "__main__":
    main()

