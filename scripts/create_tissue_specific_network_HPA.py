import argparse
import ConfigParser
import cPickle
import mysql.connector
import networkx as nx
import time
import sys, os, re

from context import NetworkAnalysis
import NetworkAnalysis.network_analysis as NA
import NetworkAnalysis.tissue_specificity as TS
import NetworkAnalysis.drug as NA_drug


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
    parser.add_argument('-ws','--workspace',dest='workspace',action = 'store',default=os.path.join(os.path.join(os.path.dirname(__file__), '..'), 'workspace'),
                        help = """Define the workspace directory where the results will be created. """)

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

    # Get the program path
    main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    pickles_dir = os.path.join(main_path, 'NetworkAnalysis/pickles')

    # Read the config file
    config_file = os.path.join(main_path, 'config.ini')
    config = ConfigParser.ConfigParser()
    config.read(config_file)

    # Connection to BIANA
    cnx = mysql.connector.connect( user=config.get('BIANA', 'user'),
                                   password=config.get('BIANA', 'password'),
                                   host=config.get('BIANA', 'host'),
                                   database=config.get('BIANA', 'database') )


    # Get Human Protein Atlas data
    expression_data_file = os.path.join(options.workspace, 'HPA_rnaseq_BIANA.tsv')

    if not fileExist(expression_data_file):
        tissue_to_UEprot_to_value = create_BIANA_expression_file(cnx, config.get('BIANA', 'unification_protocol'), expression_data_file)
    else:
        tissue_to_UEprot_to_value = read_BIANA_expression_file(expression_data_file)
    #print(tissue_to_UEprot_to_value)


    # Define the main network
    network_file = options.network_file
    type_id = 'biana'
    network_format = options.network_format
    network_instance = NA.Network(network_file, type_id, network_format)


    for tissue in tissue_to_UEprot_to_value:
        print(tissue)
        output_network_file = os.path.join(options.workspace, '{}.txt'.format(tissue))
        tissue_prots = [ UEprot for UEprot in tissue_to_UEprot_to_value[tissue] if tissue_to_UEprot_to_value[tissue][UEprot] >= 1 ]
        tissue_network = filter_network_by_tissue(network_instance, tissue_prots)
        NA.write_network_file_from_networkx_graph(tissue_network, output_network_file, 'sif', output_nodes_file=None, tissue_specific=False)


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

def create_BIANA_expression_file(cnx, unification_protocol, output_file):
    """
    Create a file containing the tissue and the corresponding expression value for BIANA user entities.
    """     

    print('\n.....Obtaining dictionary of tissue to user entity proteins to expression value (RNAseq).....\n')

    up_table = NA_drug.return_unification_protocol_table(cnx, unification_protocol)
    cursor = cnx.cursor()

    query = (''' SELECT UP.userEntityID, HT.value, RNA.value, RNA.unit 
                 FROM {} UP, externalEntityRelationParticipant RP, externalEntityRelationParticipant RT, externalEntityHumanProteinAtlas_RNAseq_value RNA, {} UT, externalEntity ET, externalEntityHumanProteinAtlas_tissue HT
                 WHERE UP.externalEntityID = RP.externalEntityID AND RP.externalEntityID != RT.externalEntityID AND RP.externalEntityRelationID = RT.externalEntityRelationID AND RT.externalEntityRelationID = RNA.externalEntityID AND RT.externalEntityID = UT.externalEntityID AND RT.externalEntityID = ET.externalEntityID AND ET.type = 'tissue' AND ET.externalEntityID = HT.externalEntityID
             '''.format(up_table, up_table))

    cursor.execute(query)

    tissue_to_UEprot_to_value = {}

    for items in cursor:

        uEprot, tissue, value, unit = items

        value = float(value)
        unit = unit.lower()

        if unit != 'tpm':
            print('Incorrect RNAseq unit for tissue {} and uE {}: {}'.format(tissue, uEprot, unit))
            sys.exit(10)

        tissue_to_UEprot_to_value.setdefault(tissue, {})
        tissue_to_UEprot_to_value[tissue][uEprot] = value

    cursor.close()

    with open(output_file, 'w') as output_fd:
        output_fd.write('Tissue\tBIANA ID\tExpression (tpm)\n')
        for tissue in tissue_to_UEprot_to_value:
            for UEprot in tissue_to_UEprot_to_value[tissue]:
                expression_value = float(tissue_to_UEprot_to_value[tissue][UEprot])
                output_fd.write('{}\t{}\t{}\n'.format(tissue, UEprot, expression_value))

    print('\nProtein to HUMAN PROTEIN ATLAS (RNAseq) dictionary obtained!\n')

    return tissue_to_UEprot_to_value

def read_BIANA_expression_file(input_file):
    """
    Read a the file and obtain the dictionary.
    """
    tissue_to_UEprot_to_value = {}
    with open(input_file, 'r') as input_fd:
        first_line = input_fd.readline()
        for line in input_fd:
            tissue, UEprot, expression_value = line.strip().split('\t')
            tissue_to_UEprot_to_value.setdefault(tissue, {})
            tissue_to_UEprot_to_value[tissue][UEprot] = float(expression_value)
    return tissue_to_UEprot_to_value

def filter_network_by_tissue(network_instance, tissue_proteins):
    """
    Read a the file and obtain the dictionary.
    """
    tissue_network=nx.Graph()
    for u,v in network_instance.network.edges():
        if u in tissue_proteins and v in tissue_proteins:
            tissue_network.add_edge(u,v)
    return tissue_network


if  __name__ == "__main__":
    main()

