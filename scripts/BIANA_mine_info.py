import argparse
import ConfigParser
import networkx as nx
import pandas as pd
import sys, os, re

import biana
try: from biana import *
except: sys.exit(10)


def main():
    """
    Requires:
    - BIANA package installed
    - BIANA database available

    Usage:
    python BIANA_mine_interactions.py -i <input_file> -t <type_of_protein_identifier> -o <output_file>

    Example of command:
    python /home/quim/PHD/Projects/BIANA/mining_scripts/BIANA_mine_info.py -i /home/quim/PHD/Projects/BIANA/data/COVID_test.txt -t genesymbol -o /home/quim/PHD/Projects/BIANA/outputs/COVID_info.txt -v
    """

    options = parse_user_arguments()
    mine_biana_info(options)

def parse_user_arguments(*args, **kwds):

    parser = argparse.ArgumentParser(
        description = "Generate a protein-protein interaction network",
        epilog      = "@oliva's lab 2017")
    parser.add_argument('-i','--input_file',dest='input_file',action = 'store',
                        help = 'Input file')
    parser.add_argument('-t','--id_type',dest='id_type',action = 'store',default='uniprotaccession',
                        help = '''Type of identifier of the proteins 
                        e.g. uniprotaccession, uniprotentry, geneid, genesymbol 
                        (default is uniprotaccession)''')
    parser.add_argument('-o','--output_file',dest='output_file',action = 'store',
                        help = 'Output file')
    parser.add_argument('-v','--verbose',dest='verbose',action = 'store_true',
                        help = 'Flag to use verbose mode')
    options=parser.parse_args()

    return options


def mine_biana_info(options):
    """
    Generates an interaction network extracting information from BIANA.
    """
    # Load config file
    scripts_path = os.path.abspath(os.path.dirname(__file__))
    config_file = os.path.join(scripts_path, 'config.ini')
    config = ConfigParser.ConfigParser()
    config.read(config_file)

    # Check ID type
    id_type = options.id_type.lower()
    accepted_id_types = ['uniprotaccession', 'uniprotentry', 'genesymbol', 'geneid']
    if id_type not in accepted_id_types:
        print('Identifier type "{}" not accepted! Please, introduce one of the followings: {}'.format(id_type, ', '.join(accepted_id_types)))
        sys.exit(10)
    else:
        node_attributes = ['taxid', 'genesymbol', 'geneid', 'uniprotaccession', 'uniprotentry']
        if id_type not in node_attributes:
            node_attributes.append(id_type)
        node_attributes.append('proteinsequence')
        if options.verbose:
            print('Node attributes to search: {}'.format(', '.join(node_attributes)))

    # Read input file
    if fileExist(options.input_file):
        biological_entities = read_input_file(options.input_file)
        if options.verbose:
            print('Biological entities introduced: {}'.format(', '.join(biological_entities)))
    else:
        print('Input file {} does not exist!'.format(options.input_file))
        sys.exit(10)


    # START BIANA SESSION
    session = create_new_session( sessionID="biana_session",
                                  dbname=config.get('BIANA', 'database'),
                                  dbhost=config.get('BIANA', 'host'),
                                  dbuser=config.get('BIANA', 'user'),
                                  dbpassword=config.get('BIANA', 'password'),
                                  unification_protocol=config.get('BIANA', 'unification_protocol') )


    # CREATE A LIST WITH ALL THE SEED IDENTIFIERS
    # Example: list_input_identifiers = [("uniprotentry","ACE_YEAST"),("uniprotentry","PGH2_HUMAN"),("uniprotentry","RND3_HUMAN")]
    list_input_identifiers = [ (id_type, biological_entity) for biological_entity in biological_entities ]
    list_input_restriction_identifiers = []
    #list_input_restriction_identifiers = [("taxid",options.taxid)]
    list_input_negative_restriction_identifiers = []


    # CREATE THE SET
    proteome = session.create_new_user_entity_set(  identifier_description_list =list_input_identifiers,
                                                    attribute_restriction_list=list_input_restriction_identifiers,
                                                    negative_attribute_restriction_list = list_input_negative_restriction_identifiers,
                                                    id_type='embedded',
                                                    only_uniques=True,
                                                    new_user_entity_set_id='proteome'  )

    user_entity_ids = proteome.get_user_entity_ids()

    if options.verbose:
        print('User entity set created.')
        print('User entities selected: {}'.format(', '.join([str(user_entity) for user_entity in user_entity_ids])))


    # OUTPUT COMMANDS
    session.output_user_entity_details( user_entity_set = proteome, 
                                        user_entity_id_list = user_entity_ids, 
                                        out_method = open(options.output_file,'w').write, 
                                        attributes = node_attributes, 
                                        include_level_info = False, 
                                        include_degree_info=False, 
                                        include_tags_info=False, 
                                        include_tags_linkage_degree_info=[], 
                                        substitute_node_attribute_if_not_exists=False, 
                                        output_1_value_per_attribute=True, 
                                        output_format="tabulated", 
                                        include_command_in_rows=False, 
                                        output_only_unique_values=True
                                        )

    if options.verbose:
        print('Output file created.')


    sys.exit(0)

    # READ OUTPUT FILE
    # Parse the BIANA network output file and obtain the network (TFs and interactors)
    network = nx.Graph()
    interactors = set()
    interactors_to_tfs = {}
    interactions_df = pd.read_csv(options.output_file, sep='\t')
    id1_field = 'Participant 1 '+id_type
    id2_field = 'Participant 2 '+id_type
    for index, row in interactions_df.iterrows():
        id1 = str(row[id1_field])
        id2 = str(row[id2_field])
        if id1 != '-' and id2 != '-':
            for sub_id1 in id1.split(', '):
                for sub_id2 in id2.split(', '):
                    network.add_edge(sub_id1, sub_id2)
                    # Obtain the interactors and the pairs interactor-tf
                    if sub_id1 not in all_tfs:
                        interactors.add(sub_id1)
                    if sub_id2 not in all_tfs:
                        interactors.add(sub_id2)
                    if sub_id1 in interactors and sub_id2 in all_tfs:
                        interactors_to_tfs.setdefault(sub_id1, set()).add(sub_id2)
                    if sub_id2 in interactors and sub_id1 in all_tfs:
                        interactors_to_tfs.setdefault(sub_id2, set()).add(sub_id1)

    # Add the connections between binding sites (ordering them from 5' to 3')
    bs1 = None
    for bs2, info in sorted(binding_site_to_info.iteritems(), key=lambda (x, y): y[2], reverse = False):
        if bs1:
            network.add_edge(bs1, bs2)
        else:
            #network.add_edge("5'", bs2)
            pass
        bs1 = (bs2 + '.')[:-1] # copy the variable
    #network.add_edge(bs1, "3'")

    # Add the connections between binding sites and transcription factors
    for bs in binding_site_to_tfs:
        for tf in binding_site_to_tfs[bs]:
            network.add_edge(bs, tf)


    # Filter interactors that are interacting with at least 2 transcription factors
    filtered_interactors = set([interactor for interactor in interactors_to_tfs if len(interactors_to_tfs[interactor])>1])

    # Filter interactors that only interact with one transcription factor but interact with another interactor connected with a different transcription factor
    one_link_interactors = set()
    for u,v in network.edges():
        if u in interactors and v in interactors:
            u_tfs = interactors_to_tfs[u]
            v_tfs = interactors_to_tfs[v]
            if u_tfs != v_tfs:
                one_link_interactors.add(u)
                one_link_interactors.add(v)


    new_network = nx.Graph()
    for u,v in network.edges():
        # Skip the interactors that are not filtered
        if u in interactors and u not in filtered_interactors and u not in one_link_interactors:
            continue
        if v in interactors and v not in filtered_interactors and v not in one_link_interactors:
            continue
        new_network.add_edge(u,v)


    # CREATE JSON FILE

    output = []

    # First the nodes
    for bs in binding_site_to_tfs:
        [bs_id, family, start, end] = binding_site_to_info[bs]
        output.append('{ "data": { "id": "%s", "label":"%s", "type": "binding site" } }' % (bs, bs_id))
    for tf in all_tfs:
        output.append('{ "data": { "id": "%s", "label":"%s", "type": "transcription factor" } }' % (tf, tf))
    for interactor in filtered_interactors:
        output.append('{ "data": { "id": "%s", "label":"%s", "type": "interactor" } }' % (interactor, interactor))
    #output.append('{ "data": { "id": "5\'", "label":"5\'", "type": "indication" } }')
    #output.append('{ "data": { "id": "3\'", "label":"3\'", "type": "indication" } }')

    # Then the interactions
    for u,v in new_network.edges():
        if v in binding_site_to_info and u not in binding_site_to_info:
            # Put always the binding site as source if the other molecule is not a binding site
            output.append('{ "data": { "id": "%s-%s", "source": "%s", "target": "%s" }  }' % (v, u, v, u))
        # Put always the tfs before the interactors
        elif u in interactors and v in all_tfs:
            output.append('{ "data": { "id": "%s-%s", "source": "%s", "target": "%s" }  }' % (v, u, v, u))
        else:
            output.append('{ "data": { "id": "%s-%s", "source": "%s", "target": "%s" }  }' % (u, v, u, v))

    # Write them in the output file
    with open(options.output_file, 'w') as output_f:
        output_f.write('[{}]'.format(','.join(output)))

    return


def read_input_file(input_file):
    """
    Reads the input file.
    """
    biological_entities = set()
    with open(input_file, 'r') as input_f:
        for line in input_f:
            entity = line.strip()
            biological_entities.add(entity)
    return biological_entities


def fileExist(file):
    """
    Checks if a file exists AND is a file.
    """
    return os.path.exists(file) and os.path.isfile(file)


def create_directory(directory):
    """
    Checks if a directory exists and if not, creates it.
    """
    try:
        os.stat(directory)
    except:
        os.mkdir(directory)
    return


if  __name__ == "__main__":
    main()



