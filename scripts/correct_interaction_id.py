import argparse
import networkx as nx
import time
import sys, os



def main():

    options = parse_user_arguments()
    get_interaction_ids_of_network(options)

def parse_user_arguments(*args, **kwds):
    """
    Parses the arguments of the program
    """

    parser = argparse.ArgumentParser(
        description = "Filter network",
        epilog      = "@oliva's lab 2017")
    parser.add_argument('-n','--network_file',dest='network_file',action = 'store',
                        help = """ Input file with a protein-protein interaction network. """)
    parser.add_argument('-o','--output_file',dest='output_file',action = 'store',
                        help = """ Output file of a protein-protein interaction network. """)
    parser.add_argument('-ws','--workspace',dest='workspace',action = 'store',default=os.path.join(os.path.join(os.path.dirname(__file__), '..'), 'workspace'),
                        help = """Define the workspace directory where the results will be created. """)

    options=parser.parse_args()

    return options


#################
#################
# MAIN FUNCTION #
#################
#################

def get_interaction_ids_of_network(options):
    """
    Gets the interaction IDs of all the interactions of an input network.
    """

    # Start marker for time measure
    start = time.time()


    # #--------------------------#
    # #   PARSE THE INTACT OBO   #
    # #--------------------------#

    # # Get the program path
    # workspace_dir = options.workspace
    # intact_obo_file = os.path.join(workspace_dir, 'intact.obo')

    # with open(intact_obo_file, 'r') as intact_fd:
    #     for line in intact_fd:
    #         term = False
    #         if line.startswith('[TERM]'):
    #             term = True
    #         elif line.startswith('id')



    #-----------------------------#
    #   PARSE THE INPUT NETWORK   #
    #-----------------------------#

    # Define the main network
    network_file = options.network_file
    output_file = options.output_file
    type_id = 'biana'
    network_format = 'multi-fields'
    network_instance = parse_network(network_file, network_format, tissue_specific=False, verbose=False)

    # Analyze the main network
    print('Number of edges: {}'.format(len(network_instance.edges())))
    print('Number of nodes: {}\n'.format(len(network_instance.nodes())))


    #-----------------------------#
    #   CORRECT INTERACTION IDS   #
    #-----------------------------#

    intact_id_to_num = {}
    dip_id_to_num = {}
    for id1,id2,d in network_instance.edges(data=True):
        sources = ';'.join(d['sources'])
        method_ids = ';'.join(d['method_ids'])
        method_names = ';'.join(d['method_names'])
        pmids = ';'.join(d['pmids'])
        interactions = d['interaction_ids']
        if interactions == '-':
            continue
        for interaction_group in interactions:
            database, interaction_ids = interaction_group.split(':')
            interaction_ids = interaction_ids.split(',')
            if database == 'intact':
                for interaction_id in interaction_ids:
                    intact_id_to_num.setdefault(interaction_id, 0)
                    intact_id_to_num[interaction_id] += 1
            if database == 'dip':
                for interaction_id in interaction_ids:
                    dip_id_to_num.setdefault(interaction_id, 0)
                    dip_id_to_num[interaction_id] += 1



    #------------------------------#
    #   WRITE THE OUTPUT NETWORK   #
    #------------------------------#

    with open(output_file, 'w') as out_fd:
        for id1,id2,d in network_instance.edges(data=True):
            sources = ';'.join(d['sources'])
            method_ids = ';'.join(d['method_ids'])
            method_names = ';'.join(d['method_names'])
            pmids = ';'.join(d['pmids'])
            interactions = d['interaction_ids']
            interactions_accepted = {}
            if interactions == '-':
                dbids = '-'
            else:
                for interaction_group in interactions:
                    database, interaction_ids = interaction_group.split(':')
                    interaction_ids = interaction_ids.split(',')
                    if database == 'dip':
                        for interaction_id in interaction_ids:
                            if interaction_id.endswith('E'):
                                interactions_accepted.setdefault(database, []).append(interaction_id)
                                continue
                            elif interaction_id.endswith('X'):
                                continue
                            if dip_id_to_num[interaction_id] == 1:
                                interactions_accepted.setdefault(database, []).append(interaction_id)
                    elif database == 'intact':
                        for interaction_id in interaction_ids:
                            if intact_id_to_num[interaction_id] == 1:
                                interactions_accepted.setdefault(database, []).append('EBI-'+str(interaction_id))
                    else:
                        for interaction_id in interaction_ids:
                            interactions_accepted.setdefault(database, []).append(interaction_id)
                if len(interactions_accepted) > 0:
                    dbids = []
                    for database in interactions_accepted:
                        ids_str = database + ':' +','.join([str(dbid) for dbid in interactions_accepted[database]])
                        dbids.append(ids_str)
                    dbids = ';'.join(dbids)
                else:
                    dbids = '-'

            out_fd.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format( id1,id2,sources,method_ids,method_names,pmids,dbids ))


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

def parse_network(network_file, network_format, tissue_specific=False, verbose=False):
    """
    Parse the network file using the Python module NetworkX.
    It is possible to parse the network in two formats:
        - 'raw' : <node1>\t<node2>
        - 'sif' : <node1>\t<score>\t<node2>
        - 'multi-fields' : <node1>\t<node2>\t<sources>\t<method_ids>\t<method_names>\t<pmids>
        - 'multi-fields' + tissue_specific=True : <node1>\t<node2>\t<sources>\t<method_ids>\t<method_names>\t<pmids>\t<tissue_db>\t<tissue_additional>
    """

    if verbose:
        print('Parsing network...\n')

    G=nx.Graph()

    network_fd = open(network_file, 'r')

    for line in network_fd:

        if network_format == 'sif':
            (node1, score, node2) = line.strip().split()
            G.add_edge(node1,node2,score=score)

        elif network_format == 'multi-fields':
            if tissue_specific:
                (node1, node2, sources, method_ids, method_names, pubmeds, databases, additional) = line.strip().split('\t')
                sources = sources.split(';')
                method_ids = method_ids.split(';')
                method_names = method_names.split(';')
                pubmeds = pubmeds.split(';')
                databases = databases.split(';')
                additional = additional.split(';')
                G.add_edge(node1, node2, sources=sources, method_ids=method_ids, method_names=method_names, pmids=pubmeds, tissue_db=databases, tissue_additional=additional)
            else:
                fields = line.strip().split('\t')
                if len(fields) == 6:
                    (node1, node2, sources, method_ids, method_names, pubmeds) = fields
                    sources = sources.split(';')
                    method_ids = method_ids.split(';')
                    method_names = method_names.split(';')
                    pubmeds = pubmeds.split(';')
                    interaction_ids = '-'
                    G.add_edge(node1, node2, sources=sources, method_ids=method_ids, method_names=method_names, pmids=pubmeds, interaction_ids=interaction_ids)
                elif len(fields) == 7:
                    (node1, node2, sources, method_ids, method_names, pubmeds, interaction_ids) = fields
                    sources = sources.split(';')
                    method_ids = method_ids.split(';')
                    method_names = method_names.split(';')
                    pubmeds = pubmeds.split(';')
                    interaction_ids = interaction_ids.split(';')
                    G.add_edge(node1, node2, sources=sources, method_ids=method_ids, method_names=method_names, pmids=pubmeds, interaction_ids=interaction_ids)
                else:
                    print(line.strip())
                    print(len(fields))
                    print('Multi-fields not recognized! Check your file')
                    sys.exit(10)

        elif network_format == 'raw':
            (node1, node2) = line.strip().split('\t')
            G.add_edge(node1,node2)
        else:
            print('Incorrect format!')
            sys.exit(10)

    network_fd.close()

    if verbose:
        print('Parsing network... finished!\n')

    return G



if  __name__ == "__main__":
    main()
