import argparse
import time
import sys, os, re

from context import NetworkAnalysis
import NetworkAnalysis.network_analysis as NA
import NetworkAnalysis.methods_dictionaries as methods_dicts




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
    parser.add_argument('-n','--network_file',dest='network_file',action = 'store',
                        help = """ Input file with a protein-protein interaction network in SIF format. """)
    parser.add_argument('-rmethod','--rmethod',dest='restrict_method',action = 'store',
                        help = """ Input file containing a list of methods.
                                   The interactions reported at least by methods in this file will be included. """)
    parser.add_argument('-emethod','--emethod',dest='except_method',action = 'store',
                        help = """ Input file containing a list of methods.
                                   The interactions only reported by methods in this file will be excluded. """)
    parser.add_argument('-rY2H','--restricted_to_Y2H',dest='restricted_to_Y2H',action = 'store_true',
                        help = 'Flag to use interactions at least described by yeast two hybrid methods (Y2H)')
    parser.add_argument('-eAFF','--except_TAP',dest='except_TAP',action = 'store_true',
                        help = 'Flag to use all interactions except those described by affinity methods (i.e. Tandem Affinity Purification)')
    parser.add_argument('-npmid','--npmid',dest='npmid',action = 'store',
                        help = """ Filter the interactions by minimum number of Pubmeds. """)
    parser.add_argument('-ndbs','--ndbs',dest='ndbs',action = 'store',
                        help = """ Filter the interactions by minimum number of databases. """)
    parser.add_argument('-nmet','--nmethods',dest='nmethods',action = 'store',
                        help = """ Filter the interactions by minimum number of methods. """)
    parser.add_argument('-wh','--withouthippie',dest='withouthippie',action = 'store_true',
                        help = """ Exclude HIPPIE from the network. """)
    parser.add_argument('-ws','--workspace',dest='workspace',action = 'store',default=os.path.join(os.path.join(os.path.dirname(__file__), '..'), 'workspace'),
                        help = """Define the workspace directory where the results will be created. """)

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


    #----------------------#
    #   DEFINE THE PATHS   #
    #----------------------#

    # Get the program path
    main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    pickles_dir = os.path.join(main_path, 'NetworkAnalysis/pickles')



    #---------------------------#
    #   GET METHODS TO FILTER   #
    #---------------------------#

    # Get the affinity dictionary
    affinity_dict = methods_dicts.affinity_dict
    affinity=set(affinity_dict.keys())

    # Get the complementation dictionary
    complementation_dict = methods_dicts.complementation_dict
    complementation=set(complementation_dict.keys())

    method_ids_included = set()
    method_ids_excluded = set()

    if options.restrict_method:
        if fileExist(options.restrict_method):
            with open(options.restrict_method, 'r') as restrict_fd:
                for line in restrict_fd:
                    method = line.strip().split('\t')[0]
                    method_ids_included.add(method)

    if options.except_method:
        if fileExist(options.except_method):
            with open(options.except_method, 'r') as except_fd:
                for line in except_fd:
                    method = line.strip().split('\t')[0]
                    method_ids_excluded.add(method)

    if options.restricted_to_Y2H:
        for method in affinity:
            method_name = affinity[method]
            if method_name == 'two hybrid':
                method_ids_included.add(method)

    if options.except_TAP:
        for method in affinity:
            method_ids_excluded.add(method)

    print('Included methods: {}'.format(', '.join(method_ids_included)))
    print('Excluded methods: {}'.format(', '.join(method_ids_excluded)))



    #-----------------------------#
    #   PARSE THE INPUT NETWORK   #
    #-----------------------------#

    # Define the main network
    network_file = options.network_file
    type_id = 'biana'
    network_format = 'multi-fields'
    network_instance = NA.Network(network_file, type_id, network_format)

    # Analyze the main network
    print('BIANA main network')
    print('Number of edges: {}'.format(len(network_instance.get_edges())))
    print('Number of nodes: {}\n'.format(len(network_instance.get_nodes())))

    methodid2interactions = network_instance.get_methodid2interactions() # Get the methods in the network and the number of interactions per method
    numpmids2interactions = network_instance.get_numpmids2interactions() # Get the number of interactions in which we have a given number of pubmed IDs
    database2interactions = network_instance.get_database2interactions() # Get the database and the number of interactions



    #----------------------------------------#
    #   CREATE NETWORK FILTERED BY METHODS   #
    #---------------------------------------.#

    if len(method_ids_included) > 0 or len(method_ids_excluded) > 0:

        method_network_file = os.path.join(options.workspace, 'network.method.txt')

        # If using names to exclude, there may be errors because the excluded affinity methods are there in names (but not in ids)
        method_network = network_instance.filter_network_by_method(methods_excluded=None, method_ids_excluded=method_ids_excluded, methods_included=None, method_ids_included=method_ids_included, output_network_file=method_network_file, output_nodes_file=None)

        # Analyze the methods network

        print('Methods-filtered main network')
        print('Included methods: {}'.format(method_ids_included))
        print('Excluded methods: {}'.format(method_ids_excluded))
        print('Number of edges: {}'.format(len(method_network.get_edges())))
        print('Number of nodes: {}\n'.format(len(method_network.get_nodes())))



    #--------------------------------------------------#
    #   CREATE NETWORK FILTERED BY NUMBER OF PUBMEDS   #
    #--------------------------------------------------#

    if options.npmid:

        # Create a network filtered by number of pubmed IDs
        pmid_network_file = os.path.join(options.workspace, 'network.{}pmid.txt'.format(str(options.npmid)))
        min_num_pubmeds = int(options.npmid)
        pmid_network = network_instance.filter_network_by_number_pubmeds(min_num_pubmeds, output_network_file=pmid_network_file, output_nodes_file=None)

        # Analyze the pubmed network

        print('Pubmed-filtered main network')
        print('Number of edges: {}'.format(len(pmid_network.get_edges())))
        print('Number of nodes: {}\n'.format(len(pmid_network.get_nodes())))



    #----------------------------------------------------#
    #   CREATE NETWORK FILTERED BY NUMBER OF DATABASES   #
    #----------------------------------------------------#

    if options.ndbs:

        # Create a network filtered by number of pubmed IDs
        db_network_file = os.path.join(options.workspace, 'network.{}db.txt'.format(str(options.ndbs)))
        min_num_databases = int(options.ndbs)
        db_network = network_instance.filter_network_by_number_databases(min_num_databases, output_network_file=db_network_file, output_nodes_file=None)

        # Analyze the pubmed network

        print('Database-filtered main network')
        print('Number of edges: {}'.format(len(db_network.get_edges())))
        print('Number of nodes: {}\n'.format(len(db_network.get_nodes())))



    #--------------------------------------------------#
    #   CREATE NETWORK FILTERED BY NUMBER OF METHODS   #
    #--------------------------------------------------#

    if options.nmethods:

        # Create a network filtered by number of methods
        nmet_network_file = os.path.join(options.workspace, 'network.{}met.txt'.format(str(options.nmethods)))
        min_num_methods = int(options.nmethods)
        nmet_network = network_instance.filter_network_by_number_methods(min_num_methods, output_network_file=nmet_network_file, output_nodes_file=None)

        # Analyze the method network

        print('Number of methods-filtered main network')
        print('Number of edges: {}'.format(len(nmet_network.get_edges())))
        print('Number of nodes: {}\n'.format(len(nmet_network.get_nodes())))



    #-------------------------------------#
    #   CREATE NETWORK EXCLUDING HIPPIE   #
    #-------------------------------------#

    if options.withouthippie:

        # Get all databases except HIPPIE
        hippie_name = 'hippie'
        databases = [x.lower() for x in database2interactions.keys()]
        for db in databases:
            if 'hippie' in db:
                hippie_name = db
                break
        # Create network excluding HIPPIE
        withouthippie_network_file = os.path.join(options.workspace, 'network.withouthippie.txt')
        withouthippie_network = network_instance.remove_database_from_network(hippie_name, withouthippie_network_file, output_nodes_file=None, verbose=False)
        # Analyze the method network

        print('HIPPIE-excluded network')
        print('HIPPIE name: {}'.format(hippie_name))
        print('Number of edges: {}'.format(len(withouthippie_network.get_edges())))
        print('Number of nodes: {}\n'.format(len(withouthippie_network.get_nodes())))



    #---------------------------------------------------------------------------#
    #   CREATE THE UNION NETWORK OF THE FILTERED ONES BY PUBMED AND DATABASES   #
    #---------------------------------------------------------------------------#

    if options.npmid and options.nmethods:
        pmidormet_network_file = os.path.join(options.workspace, 'network.{}pmidor{}met.txt'.format(str(options.npmid), str(options.nmethods)))
        pmidormet_network = NA.get_union_network_from_two_networks(pmid_network, nmet_network, pmidormet_network_file)

        # Analyze the network
        print('Pubmed-Database-Union-filtered main network')
        print('Number of edges: {}'.format(len(pmidormet_network.get_edges())))
        print('Number of nodes: {}\n'.format(len(pmidormet_network.get_nodes())))

    if options.npmid and options.ndbs:
        pmidordb_network_file = os.path.join(options.workspace, 'network.{}pmidor{}db.txt'.format(str(options.npmid), str(options.ndbs)))
        pmidordb_network = NA.get_union_network_from_two_networks(pmid_network, db_network, pmidordb_network_file)

        # Analyze the network
        print('Pubmed-Database-Union-filtered main network')
        print('Number of edges: {}'.format(len(pmidordb_network.get_edges())))
        print('Number of nodes: {}\n'.format(len(pmidordb_network.get_nodes())))

    if options.npmid and options.ndbs and options.nmethods:
        pmidordbormet_network_file = os.path.join(options.workspace, 'network.{}pmidor{}dbor{}met.txt'.format(str(options.npmid), str(options.ndbs), str(options.nmethods)))
        pmidordbormet_network = NA.get_union_network_from_two_networks(pmidordb_network, nmet_network, pmidordbormet_network_file)

        # Analyze the network
        print('Pubmed-Database-NumMethods-Union-filtered main network')
        print('Number of edges: {}'.format(len(pmidordbormet_network.get_edges())))
        print('Number of nodes: {}\n'.format(len(pmidordbormet_network.get_nodes())))



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



