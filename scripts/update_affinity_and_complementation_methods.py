import ConfigParser
import mysql.connector
import sys, os, re

def main():

    #----------------------#
    #   DEFINE THE PATHS   #
    #----------------------#

    # Get the program path
    scripts_path = os.path.abspath(os.path.dirname(__file__))
    main_path = os.path.abspath(os.path.join(scripts_path, '..'))

    # Define the outputs
    output_file = os.path.join(scripts_path, 'methods_dictionaries.py')


    #--------------------------------------#
    #   GET INFORMATION FROM CONFIG FILE   #
    #--------------------------------------#

    # Read the config file
    config_file = os.path.join(scripts_path, 'config.ini')
    config = ConfigParser.ConfigParser()
    config.read(config_file)

    cnx = mysql.connector.connect( user = config.get('BIANA', 'user'), 
                                   password = config.get('BIANA', 'password'),
                                   host = config.get('BIANA', 'host'),
                                   database = config.get('BIANA', 'database') )

    cursor = cnx.cursor()


    #----------------------------------------------------#
    #   QUERY BIANA TO OBTAIN THE CODES OF THE METHODS   #
    #----------------------------------------------------#

    # Get the predefined affinity and complementation dictionaries
    affinity = list(affinity_dict.values())
    complementation = list(complementation_dict.values())

    # Define the query
    key_attribute_table = return_key_attribute_table(cnx, ontology_name='psimiobo')
    query = ("SELECT k.value, e.value FROM {} AS k JOIN externalEntitypsimi_name AS e ON e.externalEntityID = k.externalEntityID ".format(key_attribute_table))

    cursor.execute(query)


    # Pattern of affinity
    pa = re.compile('affinity|precipitation')
    # Pattern of complementation
    pc = re.compile('complementation|hybrid')

    affinity_str = ''
    complementation_str = ''

    affinity_str += 'affinity_dict={\n'
    complementation_str += 'complementation_dict={\n'

    for (externalEntityID, value) in cursor:
        # Write the method if it is in the predefined dicts or matches the patterns of affinity/complementation
        if value in affinity or pa.search(value) != None:
            affinity_str += "\t'{}':\t'{}',\n".format(externalEntityID, value)
            #fa.write("'{}':\t'{}',\n".format(externalEntityID, value))
        if value in complementation or pc.search(value) != None:
            complementation_str += "\t'{}':\t'{}',\n".format(externalEntityID, value)
            #fc.write("'{}':\t'{}',\n".format(externalEntityID, value))

    affinity_str += '}\n'
    complementation_str += '}\n'

    cursor.close()
    cnx.close()


    #---------------------------#
    #   WRITE THE OUTPUT FILE   #
    #---------------------------#

    with open(output_file, 'w') as output_file_fd:

        for line in affinity_str:
            output_file_fd.write(line)
        output_file_fd.write('\n')
        for line in complementation_str:
            output_file_fd.write(line)

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

if  __name__ == "__main__":

    affinity_dict={
    '0':    'molecular interaction',
    '4':    'affinity chromatography technology',
    '6':    'anti bait coimmunoprecipitation',
    '7':    'anti tag coimmunoprecipitation',
    '8':    'array technology',
    '9':    'bacterial display',
    '19':   'coimmunoprecipitation',
    # '27':   'cosedimentation',
    # '28':   'cosedimentation in solution',
    # '29':   'cosedimentation through density gradient',
    '34':   'display technology',
    '47':   'far western blotting',
    '48':   'filamentous phage display',
    '49':   'filter binding',
    '50':   'flag tag coimmunoprecipitation',
    '60':   'ha tag coimmunoprecipitation',
    '62':   'his tag coimmunoprecipitation',
    '66':   'lambda phage display',
    '71':   'molecular sieving',
    '73':   'mrna display',
    '75':   'myc tag coimmunoprecipitation',
    '81':   'peptide array',
    '84':   'phage display',
    '89':   'protein array',
    '92':   'protein in situ array',
    '95':   'proteinchip(r) on a surface-enhanced laser desorption/ionization',
    '96':   'pull down',
    '98':   'ribosome display',
    '108':  't7 phage display',
    '109':  'tap tag coimmunoprecipitation',
    '115':  'yeast display',
    '225':  'chromatin immunoprecipitation array',
    '400':  'affinity technology',
    '402':  'chromatin immunoprecipitation assay',
    '405':  'competition binding',
    '411':  'enzyme linked immunosorbent assay',
    '412':  'electrophoretic mobility supershift assay',
    '413':  'electrophoretic mobility shift assay',
    '440':  'saturation binding',
    '463':  'biogrid', 
    '469':  'intact', 
    '471':  'mint', 
    '492':  'in vitro',
    '493':  'in vivo',
    '657':  'systematic evolution of ligands by exponential enrichment',
    '676':  'tandem affinity purification',
    '678':  'antibody array',
    '686':  'unspecified method',
    '695':  'sandwich immunoassay',
    '729':  'luminescence based mammalian interactome mapping',
    '858':  'immunodepleted coimmunoprecipitation',
    '892':  'solid phase assay',
    '899':  'p3 filamentous phage display',
    '900':  'p8 filamentous phage display',
    '921':  'surface plasmon resonance array',
    '946':  'miniaturized immunoprecipitation',
    '947':  'bead aggregation assay',
    '963':  'interactome parallel affinity capture',
    '1017': 'rna immunoprecipitation',
    '1028': 'modified chromatin immunoprecipitation',
    '1029': 'proteomics of isolated chromatin segments',
    '1031': 'protein folding/unfolding',
    '1087': 'monoclonal antibody blockade',
    }

    complementation_dict={
    '0':    'molecular interaction',
    '10':   'beta galactosidase complementation',
    '11':   'beta lactamase complementation',
    '14':   'adenylate cyclase complementation',
    '18':   'two hybrid',
    '80':   'partial DNA sequence identification by hybridization',
    '90':   'protein complementation assay',
    '97':   'reverse ras recruitment system',
    '111':  'dihydrofolate reductase reconstruction',
    '112':  'ubiquitin reconstruction',
    '228':  'cytoplasmic complementation assay',
    '230':  'membrane bound complementation assay',
    '231':  'mammalian protein protein interaction trap',
    '232':  'transcriptional complementation assay',
    '369':  'lex-a dimerization assay',
    '370':  'tox-r dimerization assay',
    '397':  'two hybrid array',
    '398':  'two hybrid pooling approach',
    '399':  'two hybrid fragment pooling approach',
    '432':  'one hybrid',
    '437':  'protein three hybrid',
    '438':  'rna three hybrid',
    '492':  'in vitro',
    '493':  'in vivo',
    '588':  'three hybrid',
    '655':  'lambda repressor two hybrid',
    '726':  'reverse two hybrid',
    '727':  'lexa b52 complementation',
    '728':  'gal4 vp16 complementation',
    '809':  'bimolecular fluorescence complementation',
    '895':  'protein kinase A complementation',
    '916':  'lexa vp16 complementation',
    '1037': 'Split renilla luciferase complementation',
    '1111': 'two hybrid bait or prey pooling approach',
    '1112': 'two hybrid prey pooling approach',
    '1113': 'two hybrid bait and prey pooling approach',
    '1203': 'split luciferase complementation',
    '1204': 'split firefly luciferase complementation',
    '1320': 'membrane yeast two hybrid',
    }

    main()