import ConfigParser
import mysql.connector
import sys, os, re

from context import NetworkAnalysis
import NetworkAnalysis.network_analysis as NA

def main():
    """
    module load Python/2.7.11
    python /home/quim/PHD/Projects/NetworkAnalysis/scripts/design_file_to_filter_by_methods.py
    """

    main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    hippie_scores_file = os.path.join(main_path, 'NetworkAnalysis/hippie_scores/experimental_scores.tsv')
    network_file = '/home/quim/PHD/Projects/BIANA/outputs/geneID_seqtax_PROUS/human_network_prous_2020.txt'
    physical_methods_file = '/home/quim/PHD/Projects/NetworkAnalysis/workspace/restrict_method_phy_prous2020.txt'
    functional_methods_file = '/home/quim/PHD/Projects/NetworkAnalysis/workspace/restrict_method_func_prous2020.txt'


    # Read the config file
    config_file = os.path.join(main_path, 'config.ini')
    config = ConfigParser.ConfigParser()
    config.read(config_file)

    # Connect to BIANA
    cnx = mysql.connector.connect( user=config.get('BIANA', 'user'),
                                   password=config.get('BIANA', 'password'),
                                   host=config.get('BIANA', 'host'),
                                   database=config.get('BIANA', 'database') )


    # Parse psimiobo
    key_attribute_table = return_key_attribute_table(cnx, ontology_name='psimiobo')
    psimi2method = obtain_psimi_to_method(cnx, key_attribute_table)
    # Parse hippie scores file
    psimi2score = parse_hippie_scores(hippie_scores_file)


    # Parse network
    # Define the main network
    #type_id = 'biana'
    #network_format = 'multi-fields'
    #network_instance = NA.Network(network_file, type_id, network_format)


    # Analyze the main network
    #print('Number of edges: {}'.format(len(network_instance.get_edges())))
    #print('Number of nodes: {}\n'.format(len(network_instance.get_nodes())))


    # Get the methods in the network and the number of interactions per method
    y2h_methods = set()
    skip_methods = set([80, 432, 437, 438, 588, 655, 726, 2277])
    #methodid2interactions = network_instance.get_methodid2interactions()
    #for psimi, interactions in sorted(methodid2interactions.iteritems(), key=lambda (x, y): y, reverse=True):
    for psimi in psimi2method:
        name = '-'
        if psimi in psimi2method:
            name = psimi2method[psimi]
            if 'hybrid' in name and psimi not in skip_methods:
                y2h_methods.add(psimi)


    # Get methods from HIPPIE
    hippie_methods = set()
    for psimi, score in sorted(psimi2score.iteritems(), key=lambda (x, y): y, reverse=True):
        if score >= 7:
            hippie_methods.add(psimi)


    # Write PHYSICAL methods file
    phy_methods = hippie_methods | y2h_methods
    with open(physical_methods_file, 'w') as phy_fd:
        for psimi in sorted(phy_methods):
            name = '-'
            if psimi in psimi2method:
                name = psimi2method[psimi]
            phy_fd.write('{}\t{}\n'.format(psimi, name))


    # Write FUNCTIONAL methods file
    coimmuno_methods = set([6])
    func_methods = hippie_methods | y2h_methods | coimmuno_methods
    with open(functional_methods_file, 'w') as func_fd:
        for psimi in sorted(func_methods):
            name = '-'
            if psimi in psimi2method:
                name = psimi2method[psimi]
            func_fd.write('{}\t{}\n'.format(psimi, name))


    # Close BIANA connection
    cnx.close()

    return


def return_key_attribute_table(cnx, ontology_name):
    """
    Returns the table that contains the Unification Protocol
    introduced as query
    """

    query = (''' SELECT key_id FROM ExternalEntityOntology WHERE name = %s ''')

    cursor = cnx.cursor()
    cursor.execute(query, (ontology_name,))
    key_ids = []
    for items in cursor:
        for key in items:
            key_ids.append(key)
    key_id = key_ids[0]
    key_attribute_table = 'key_attribute_{}'.format(str(key_id))
    cursor.close()

    return key_attribute_table


def obtain_psimi_to_method(cnx, key_attribute_table):
    """
    Obtain dictionary uE_prot : {'psi-mi code' : 'method_name' }
    """     

    print('\n.....Obtaining dictionary of PSI-MI codes to Method names.....\n')

    cursor = cnx.cursor()

    query = (''' SELECT K.value, P.value FROM externalEntitypsimi_name P, {} K where P.externalEntityID = K.externalEntityID
             '''.format(key_attribute_table))

    cursor.execute(query)

    psimi2method = {}

    for items in cursor:

        psimi = items[0]
        method = items[1]
        psimi2method[psimi] = method

    cursor.close()

    print('\nPSI-MI 2 METHOD dictionary obtained!\n')

    return psimi2method


def parse_hippie_scores(hippie_scores_file):
    """
    Obtain dictionary uE_prot : {'psi-mi code' : 'method_name' }
    """     

    print('\n.....Parsing HIPPIE scores.....\n')

    psimi2score = {}

    mi_regex = re.compile('MI:([0-9]{4})')
    hippie_scores_fd = open(hippie_scores_file, 'r')

    for line in hippie_scores_fd:
        method_name, psimi, score = line.strip().split('\t')
        m = mi_regex.search(psimi)
        if m:
            psimi = int(m.group(1))
            psimi2score[psimi]=float(score)

    hippie_scores_fd.close()

    print('\n.....Parsing of HIPPIE scores done!.....\n')

    return psimi2score


def fileExist(file):
    """
    Checks if a file exists AND is a file
    """
    return os.path.exists(file) and os.path.isfile(file)


if  __name__ == "__main__":
    main()
