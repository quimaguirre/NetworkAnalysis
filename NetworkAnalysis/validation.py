import sys, os, re
import ConfigParser
import networkx as nx

import network_analysis as NA


class ConsensusPathDBParser(object):

    def __init__(self, ConsensusPathDB_file):

        self.ConsensusPathDB_file = ConsensusPathDB_file

        return

    def parse(self, output_network_file, output_network_format):
        """
        Parse the PPI network from ConsensusPathDB and return
        a network in raw/sif/multi-fields format.
        """
        G=nx.Graph()
        with open(self.ConsensusPathDB_file, 'r') as CPDB_fd:
            for line in CPDB_fd:
                if line[0] == '#':
                    continue

                srcs, pubmeds, participants, score = line.strip().split('\t')

                sources = srcs.lower().split(',')
                if pubmeds == '':
                    pubmeds = '-'
                else:
                    pubmeds = pubmeds.split(',')
                participants = participants.split(',')
                if score == 'NA':
                    continue
                else:
                    score = float(score)

                # Get the pairs from the participants
                pairs = set()
                if len(participants) == 1:
                    pairs.add(frozenset([participants[0],participants[0]]))
                elif len(participants) > 1:
                    for protein1 in participants:
                        for protein2 in participants:
                            if protein1 != protein2:
                                pairs.add(frozenset([ protein1, protein2 ]))

                # Add the pairs in the graph
                for pair in pairs:
                    node1, node2 = pair
                    if output_network_format == 'sif':
                        G.add_edge(node1,node2,score=score)
                    elif output_network_format == 'multi-fields':
                        G.add_edge(node1, node2, sources=sources, method_ids='-', method_names='-', pmids=pubmeds)
                    elif output_network_format == 'raw':
                        G.add_edge(node1,node2)

        # Write the network
        NA.write_network_file_from_networkx_graph(G, output_network_file, output_network_format, output_nodes_file=None, tissue_specific=False)

        return NA.Network(output_network_file, 'uniprotentry', output_network_format)



class I2DParser(object):

    def __init__(self, input_file):

        self.input_file = input_file

        return

    def parse(self, output_network_file, output_network_format):
        """
        Parse the PPI network from I2D and return
        a network in raw/sif/multi-fields format.
        """
        interaction_to_srcs = {}
        with open(self.input_file, 'r') as I2D_fd:
            first_line = I2D_fd.readline()
            for line in I2D_fd:

                src, node1, node2 = line.strip().split('\t')

                interaction = frozenset([node1, node2])
                interaction_to_srcs.setdefault(interaction, set())
                interaction_to_srcs[interaction].add(src.lower())

        G=nx.Graph()
        for interaction in interaction_to_srcs:
            if len(interaction) == 1:
                node1 = list(interaction)[0]
                node2 = list(interaction)[0]
            elif len(interaction) == 2:
                node1, node2 = interaction
            sources = interaction_to_srcs[interaction]
            if output_network_format == 'sif':
                G.add_edge(node1,node2,score=1.000000)
            elif output_network_format == 'multi-fields':
                G.add_edge(node1, node2, sources=sources, method_ids='-', method_names='-', pmids='-')
            elif output_network_format == 'raw':
                G.add_edge(node1,node2)

        # Write the network
        NA.write_network_file_from_networkx_graph(G, output_network_file, output_network_format, output_nodes_file=None, tissue_specific=False)

        return NA.Network(output_network_file, 'uniprotaccession', output_network_format)



class HippieParser(object):

    def __init__(self, input_file):

        self.hippie_file = input_file

        self.interactors = set()
        self.interactor2id2altid = {}
        self.interactor2taxid = {}

        self.interactions = set()
        self.interaction2interactor = {}
        self.interaction2methods = {}
        self.interaction2pubmeds = {}
        self.interaction2sources = {}
        self.interaction2score = {}

        self.formats = ['sif', 'raw', 'multi-fields']

        return



    def parse(self):

        print("\n.....PARSING THE HIPPIE INTERACTOME.....\n")

        hippie_file_fd = open(self.hippie_file,'r')

        first_line = hippie_file_fd.readline()

        # Obtain a dictionary: "field_name" => "position"
        fields_dict = self.obtain_header_fields(first_line)
        #ID Interactor A    ID Interactor B Alt IDs Interactor A    Alt IDs Interactor B    Aliases Interactor A    Aliases Interactor B    Interaction Detection Methods   Publication 1st Author  Publication Identifiers Taxid Interactor A  Taxid Interactor B  Interaction Types   Source Databases    Interaction Identifiers Confidence Value    Presence In Other Species


        for line in hippie_file_fd:

            # Split the line in fields
            fields = line.strip().split("\t")

            # entrez gene:216 entrez gene:216 uniprotkb:AL1A1_HUMAN   uniprotkb:AL1A1_HUMAN   -   -   MI:0493(in vivo)|MI:0018(two hybrid)    Rodriguez-Zavala JS (2002)|Rual JF (2005)|Rolland T (2014)  pubmed:12081471|pubmed:16189514|pubmed:25416956 taxid:9606(Homo sapiens)    taxid:9606(Homo sapiens)    -   MI:0468(hprd)|biogrid|MI:0469(intact)|MI:0471(mint)|i2d|rual05  -   0.76    

            #### Obtain IDENTIFIERS ####

            # Obtain the fields of interest
            identifier1 = fields[ fields_dict['ID Interactor A'] ]
            identifier2 = fields[ fields_dict['ID Interactor B'] ]
            alt_ident1 = fields[ fields_dict['Alt IDs Interactor A'] ]
            alt_ident2 = fields[ fields_dict['Alt IDs Interactor B'] ]

            # If we do not have identifiers, we stop the parsing
            if identifier1 == '-' and alt_ident1 == '-':
                print('Identifier unknown')
                print('Identifier: {} Alt identifier: {}'.format(identifier1, alt_ident1))
                sys.exit(10)
            if identifier2 == '-' and alt_ident2 == '-':
                print('Identifier unknown')
                print('Identifier: {} Alt identifier: {}'.format(identifier2, alt_ident2))
                sys.exit(10)

            # If we do not have main identifier, we use the alternative as main
            if identifier1 == '-' and alt_ident1 != '-': # For interactor 1
                id1 = alt_ident1
                alt_id1 = identifier1
            else:
                id1 = identifier1
                alt_id1 = alt_ident1
            if identifier2 == '-' and alt_ident2 != '-': # For interactor 2
                id2 = alt_ident2
                alt_id2 = identifier2
            else:
                id2 = identifier2
                alt_id2 = alt_ident2

            #### Obtain TAXID ####

            taxid1 = fields[ fields_dict['Taxid Interactor A'] ]
            taxid2 = fields[ fields_dict['Taxid Interactor B'] ]
            self.check_id(taxid1, 'taxid')
            self.check_id(taxid2, 'taxid')
            taxid1 = taxid1[len('taxid:'):].split('(')[0] # Get only the taxonomy ID: taxid:9606(Homo sapiens) --> 9606
            taxid2 = taxid2[len('taxid:'):].split('(')[0] # Get only the taxonomy ID: taxid:9606(Homo sapiens) --> 9606

            #### Obtain METHOD ID ####

            if fields[ fields_dict['Interaction Detection Methods'] ] == '-' or fields[ fields_dict['Interaction Detection Methods'] ] == '':
#### ---------> IF THERE IS NO METHOD, WE SKIP THE INTERACTION!!! ######
                continue
            else:
                all_methods = fields[ fields_dict['Interaction Detection Methods'] ].split('|') # MI:0493(in vivo)|MI:0018(two hybrid)
                methods = []
                for method in all_methods:
                    if method.startswith('MI:'):
                        methods.append(method[len('MI:'):].split('(')[0]) # Get only the method ID
                    else:
                        #print('Method does not start with MI: {}'.format(method))
                        #sys.exit(10)
                        pass

            #### Obtain PUBMED ID ####

            if fields[ fields_dict['Publication Identifiers'] ] == '-' or fields[ fields_dict['Publication Identifiers'] ] == '':
                pubmeds = '-' 
            else:
                pubmeds = fields[ fields_dict['Publication Identifiers'] ].split('|') # pubmed:12081471|pubmed:16189514|pubmed:25416956
                [ self.check_id(x, 'pubmed') for x in pubmeds ]
                pubmeds = [ x[len('pubmed:'):] for x in pubmeds ] # Get only the pubmed ID

            #### Obtain SOURCE DATABASES ID ####

            if fields[ fields_dict['Source Databases'] ] == '-' or fields[ fields_dict['Source Databases'] ] == '':
                sources = '-'
                print('No source for ids {} {}'.format(id1, id2))
            else:
                sources = fields[ fields_dict['Source Databases'] ].split('|') # MI:0468(hprd)|biogrid|MI:0469(intact)|MI:0471(mint)|i2d|rual05

            #### Obtain CONFIDENCE VALUE ####

            if fields[ fields_dict['Confidence Value'] ] == '-':
                score = '-'
                print('No score')
                sys.exit(10)
            else:
                score = float(fields[ fields_dict['Confidence Value'] ])
                if score < 0.5:
#### -------------> IF THE SCORE IS BELOW 0.5, WE SKIP THE INTERACTION ######
                    continue

            #### CREATE INTERACTION ID ####

            # Create an interaction id for the protein-protein interaction
            # ---> interaction id = identifier1 + '---' + geneid2
            interaction_1 = id1 + '---' + id2
            interaction_2 = id2 + '---' + id1
            interaction_3 = alt_id1 + '---' + alt_id2
            interaction_4 = alt_id2 + '---' + alt_id1

            #### CREATE INTERACTOR ID ####

            # We create an interactor id composed by the id + the alt_id
            interactor_1 = id1 + '---' + alt_id1
            interactor_2 = id2 + '---' + alt_id2

            #### INSERT THE FIELDS INTO DICTIONARIES

            # Check if interaction was already reported
            if not interaction_1 in self.interactions:
                interaction = interaction_1
            # If the identifier of the interaction is already used, we assign another one
            elif not interaction_2 in self.interactions:
                interaction = interaction_2
            elif not interaction_3 in self.interactions:
                interaction = interaction_3
            elif not interaction_4 in self.interactions:
                interaction = interaction_4
            # If all the possible identifiers are used, we stop... but this has not happened :)
            else:
                print('This interaction is already reported: {}'.format(interaction_1))
                print('The other way around was also reported: {}'.format(interaction_2))
                sys.exit(10)

            # Add the interaction
            self.interactions.add(interaction)

            # Add identifiers of the interaction
            self.interactors.add(interactor_1)
            self.interactors.add(interactor_2)

            self.interactor2id2altid.setdefault(interactor_1, {})
            self.interactor2id2altid[interactor_1]['id'] = id1
            self.interactor2id2altid[interactor_1]['altid'] = alt_id1

            self.interactor2id2altid.setdefault(interactor_2, {})
            self.interactor2id2altid[interactor_2]['id'] = id2
            self.interactor2id2altid[interactor_2]['altid'] = alt_id2

            self.interaction2interactor.setdefault(interaction, set())
            self.interaction2interactor[interaction].add(interactor_1)
            self.interaction2interactor[interaction].add(interactor_2)

            # Insert the taxIDs of the interactors
            self.interactor2taxid[interactor_1] = taxid1
            self.interactor2taxid[interactor_2] = taxid2

            # Insert methods of the interaction
            self.interaction2methods[interaction] = methods

            # Insert pubmeds of the interaction
            self.interaction2pubmeds[interaction] = pubmeds

            # Insert sources of the interaction
            self.interaction2sources[interaction] = sources

            # Insert score of the interaction
            self.interaction2score[interaction] = score

        hippie_file_fd.close()

        return


    def write_network_file(self, output_network_file, output_network_format, type_id):
        """
        Output the HIPPIE network in one of the following formats:
        raw/sif/multi-fields
        And it can be in the following types of ID: 
        uniprotentry or geneid
        """

        # Check the input parameters
        type_id = type_id.lower()
        if type_id != 'uniprotentry' and type_id != 'geneid':
            raise NA.IncorrectTypeID(type_id)
        output_network_format = output_network_format.lower()
        if output_network_format != 'sif' and output_network_format != 'raw' and output_network_format != 'multi-fields':
            raise NA.IncorrectNetworkFormat(output_network_format, self.formats)

        # Create the graph
        G=nx.Graph()
        for interaction in self.interactions:
            partners = self.interaction2interactor[interaction]

            if len(partners) == 1:
                node1 = list(partners)[0]
                node2 = list(partners)[0]
            elif len(partners) == 2:
                node1, node2 = partners
            else:
                print('Caution! Not one/two partners: {}'.format(partners))
                sys.exit(10)

            id1 = self.find_identifier(node1, type_id)
            id2 = self.find_identifier(node2, type_id)
            if not id1 or not id2:
                print('Skipped interaction: {} and {}'.format(node1, node2))
                continue

            methods = self.interaction2methods[interaction]
            pubmeds = self.interaction2pubmeds[interaction]
            sources = self.interaction2sources[interaction]
            score = float(self.interaction2score[interaction])

            source_regex = re.compile('MI:[0-9]{4}\(([a-zA-Z]+)\)')
            if sources != '-':
                src_names = set()
                for source in sources:
                    m = source_regex.search(source)
                    if m:
                        source = m.group(1)
                        src_names.add(source)
                    else:
                        src_names.add(source)
                sources = src_names

            if output_network_format == 'sif':
                G.add_edge(id1,id2,score=score)
            elif output_network_format == 'multi-fields':
                G.add_edge(id1, id2, sources=sources, method_ids=methods, method_names='-', pmids=pubmeds)
            elif output_network_format == 'raw':
                G.add_edge(id1, id2)

        # Write the network
        NA.write_network_file_from_networkx_graph(G, output_network_file, output_network_format, output_nodes_file=None, tissue_specific=False)

        return NA.Network(output_network_file, type_id, output_network_format)


    def find_identifier(self, interactor, input_type_id):
        hippie_type_id_to_correct_id = {'entrez gene' : 'geneid', 'uniprotkb' : 'uniprotentry'}
        interactor_id = self.interactor2id2altid[interactor]['id'] # Get the interactor id
        interactor_altid = self.interactor2id2altid[interactor]['altid'] # Get the interactor alternative id
        for identifier in (interactor_id, interactor_altid):
            if identifier == '-':
                continue
            hippie_type_id, value_id = self.analyze_identifier(identifier)
            if input_type_id.lower() == hippie_type_id_to_correct_id[hippie_type_id]:
                return value_id
        return None

    def obtain_header_fields(self, first_line):
        """ 
        Obtain a dictionary: "field_name" => "position" 
        """
        fields_dict = {}

        header_fields = first_line.strip().split("\t")
        for x in xrange(0, len(header_fields)):
            fields_dict[header_fields[x]] = x

        return fields_dict

    def check_id(self, identifier, type_id):
        """ 
        Checks if there is an id 
        """
        if identifier.startswith(type_id):
            identifier = identifier.split(type_id+':')[1]
        else:
            print('Not {} --> {}'.format(type_id,identifier))
            sys.exit(10)

        return

    def analyze_identifier(self, identifier):
        """ 
        From the identifier of an interactor, we return the type of identifier and the value
        Example ---> entrez gene:1029 ---> 'entrez gene', '1029'
        """
        type_identifier, value_identifier = identifier.split(':')

        return type_identifier, value_identifier



def translate_network_from_BIANA(network_instance, input_type_id, output_type_id, output_network_file):
    """
    Translate a network by querying BIANA.
    """
    import mysql.connector

    format_to_table = {
        'geneid' : 'externalEntityGeneID',
        'uniprotentry' : 'externalEntityUniprotEntry',
        'uniprotaccession' : 'externalEntityUniprotAccession',
        'genesymbol' : 'externalEntityGeneSymbol'
    }
    input_type_id = input_type_id.lower()
    if input_type_id not in format_to_table:
        raise NA.IncorrectNetworkFormat(input_type_id, format_to_table.keys())
    output_type_id = output_type_id.lower()
    if output_type_id not in format_to_table:
        raise NA.IncorrectNetworkFormat(output_type_id, format_to_table.keys())

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

    cursor = cnx.cursor()

    node_to_translation = {}
    for node in network_instance.get_nodes():

        query = (''' SELECT O.value
                     FROM {} I, {} O
                     WHERE I.externalEntityID = O.externalEntityID AND O.type = "unique" AND I.value = %s
                 '''.format(format_to_table[input_type_id], format_to_table[output_type_id]))

        cursor.execute(query, (node,))

        for items in cursor:

            translation = items[0]
            node_to_translation.setdefault(node, set())
            node_to_translation[node].add(translation)

    cursor.close()

    G=nx.Graph()
    if network_instance.network_format == 'sif':
        try:
            for u,v,d in network_instance.network.edges_iter():
                if 'score' in d:
                    score = d['score']
                else:
                    score = 1.000000
                if u in node_to_translation and v in node_to_translation:
                    trans1 = node_to_translation[u]
                    trans2 = node_to_translation[v]
                    for t1 in trans1:
                        for t2 in trans2:
                            G.add_edge(t1,t2,score=score)

        except:
            score = 1.000000
            for u,v in network_instance.network.edges_iter():
                if u in node_to_translation and v in node_to_translation:
                    trans1 = node_to_translation[u]
                    trans2 = node_to_translation[v]
                    for t1 in trans1:
                        for t2 in trans2:
                            G.add_edge(t1,t2,score=score)

    elif network_instance.network_format == 'multi-fields':
        for u,v,d in network_instance.network.edges_iter(data=True):
            if u in node_to_translation and v in node_to_translation:
                trans1 = node_to_translation[u]
                trans2 = node_to_translation[v]
                for t1 in trans1:
                    for t2 in trans2:
                        G.add_edge(t1,t2,d)

    elif network_instance.network_format == 'raw':
        for u,v in network_instance.network.edges_iter():
            if u in node_to_translation and v in node_to_translation:
                trans1 = node_to_translation[u]
                trans2 = node_to_translation[v]
                for t1 in trans1:
                    for t2 in trans2:
                        G.add_edge(t1,t2)

    # Write the network
    NA.write_network_file_from_networkx_graph(G, output_network_file, network_instance.network_format, output_nodes_file=None, tissue_specific=False)

    return NA.Network(output_network_file, output_type_id, network_instance.network_format)
