import argparse
import ConfigParser
import mysql.connector
import time
import sys, os, re

from context import NetworkAnalysis
import NetworkAnalysis.network_analysis as NA




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


    #----------------------#
    #   DEFINE THE PATHS   #
    #----------------------#

    # Get the program path
    main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))



    #-----------------------------#
    #   PARSE THE INPUT NETWORK   #
    #-----------------------------#

    # Define the main network
    network_file = options.network_file
    output_file = options.output_file
    type_id = 'biana'
    network_format = 'multi-fields'
    network_instance = NA.Network(network_file, type_id, network_format)

    # Analyze the main network
    print('Number of edges: {}'.format(len(network_instance.get_edges())))
    print('Number of nodes: {}\n'.format(len(network_instance.get_nodes())))


    # Get the program path
    main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

    # Read the config file
    config_file = os.path.join(main_path, 'config.ini')
    config = ConfigParser.ConfigParser()
    config.read(config_file)

    cnx = mysql.connector.connect( user=config.get('BIANA', 'user'),
                                   password=config.get('BIANA', 'password'),
                                   host=config.get('BIANA', 'host'),
                                   database=config.get('BIANA', 'database') )

    up_table = NA.return_unification_protocol_table(cnx, config.get('BIANA', 'unification_protocol'))

    with open(output_file, 'w') as out_fd:
        for id1,id2,d in network_instance.network.edges(data=True):
            sources = ';'.join(d['sources'])
            method_ids = ';'.join(d['method_ids'])
            method_names = ';'.join(d['method_names'])
            pmids = ';'.join(d['pmids'])
            database_to_databaseIDs = get_id_of_interaction(cnx, up_table, id1, id2)
            dbids = []
            for database in database_to_databaseIDs:
                ids_str = database + ':' +','.join([str(dbid) for dbid in database_to_databaseIDs[database]])
                dbids.append(ids_str)
            dbids = ';'.join(dbids)
            out_fd.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format( id1,id2,sources,method_ids,method_names,pmids,dbids ))

    cnx.close()

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


def get_id_of_interaction(cnx, up_table, id1, id2):
    """
    Get the ID of an interaction.
    """

    database_to_query = {
        'intact' : 'SELECT value FROM externalEntityIntAct WHERE externalEntityID = %s',
        'dip' : 'SELECT value FROM externalEntityDIP WHERE externalEntityID = %s',
        'innatedb' : 'SELECT value FROM externalEntityInnateDB WHERE externalEntityID = %s',
    }

    cursor = cnx.cursor()

    query = (''' SELECT R1.externalEntityRelationID, DB.databaseName 
                 FROM {} U1, {} U2, externalEntityRelationParticipant R1, externalEntityRelationParticipant R2, externalEntity E, externalDatabase DB 
                 WHERE U1.externalEntityID = R1.externalEntityID AND U2.externalEntityID = R2.externalEntityID AND R1.externalEntityRelationID = R2.externalEntityRelationID AND R1.externalEntityRelationID = E.externalEntityID AND E.externalDatabaseID = DB.externalDatabaseID AND U1.userEntityID = %s AND U2.userEntityID = %s;
             '''.format(up_table, up_table))

    cursor.execute(query, (id1, id2,))

    database_to_relationIDs = {}
    for items in cursor:
        relationID, database = items
        database_to_relationIDs.setdefault(database, set()).add(relationID)

    database_to_databaseIDs = {}
    for database in database_to_relationIDs:
        if database in database_to_query:
            query = ('{}'.format(database_to_query[database]))
            for relationID in database_to_relationIDs[database]:
                cursor.execute(query, (relationID,))
                for items in cursor:
                    databaseID = items[0]
                    database_to_databaseIDs.setdefault(database, set()).add(databaseID)

    cursor.close()

    return database_to_databaseIDs


if  __name__ == "__main__":
    main()
