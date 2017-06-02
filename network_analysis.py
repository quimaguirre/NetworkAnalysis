import cPickle
import networkx as nx
import os, sys

import NetworkAnalysis.network_translation as NT
import NetworkAnalysis.tissue_specificity as TS

"""
    NetworkAnalysis
    Copyright (C) 2017 Joaquim Aguirre-Plans 
    Structural Bioinformatics Laboratory
"""

class Network(object):
    """ 
    Class defining a network object 
    """

    def __init__(self, network_file, node_file, type_id, network_format):
        """ 
        @param:    network_file
        @pdef:     Path to the file containing the edges of the network
        @ptype:    {String}

        @param:    node_file
        @pdef:     Path to the file containing the nodes of the network
        @ptype:    {String}

        @param:    type_id
        @pdef:     Type of IDs in the network
        @ptype:    {String}

        @param:    network_format
        @pdef:     Format of the network
        @ptype:    {String} {'sif' or 'multi-fields'}

        @raises: {IncorrectNetworkFormat} if the network_format is not in
        self.formats.
        @raises: {IncorrectTypeID} if the method translate_network is used with
        a network of type_id different from 'biana'
        """

        self.network_file = network_file
        self.node_file = node_file
        self.type_id = type_id
        self.network_format = network_format
        self.formats = ['sif', 'multi-fields']
        self.executable_path = '/home/quim/project/tissue_specificity/scripts'

        self._tissue_specific = False

        self.network = self.parse_network(tissue_specific=self._tissue_specific)

    ###########
    # METHODS #
    ###########

    def get_edges(self):
        return self.network.edges()

    def get_nodes(self):
        return self.network.nodes()

    def parse_network(self, tissue_specific=False):
        """
        Parse the network file using the Python module NetworkX.
        It is possible to parse the network in two formats:
            - 'sif' : <node1>\t<score>\t<node2>
            - 'multi-fields' : <node1>\t<node2>\t<sources>\t<method_ids>\t<method_names>\t<pmids>
            - 'multi-fields' + tissue_specific=True : <node1>\t<node2>\t<sources>\t<method_ids>\t<method_names>\t<pmids>\t<tissue_db>\t<tissue_additional>
        """

        print('Parsing network...\n')

        G=nx.Graph()

        network_fd = open(self.network_file, 'r')

        for line in network_fd:

            if self.network_format == 'sif':
                (node1, score, node2) = line.strip().split('\t')
                G.add_edge(node1,node2,score=score)

            elif self.network_format == 'multi-fields':
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
                    (node1, node2, sources, method_ids, method_names, pubmeds) = line.strip().split('\t')
                    sources = sources.split(';')
                    method_ids = method_ids.split(';')
                    method_names = method_names.split(';')
                    pubmeds = pubmeds.split(';')
                    G.add_edge(node1, node2, sources=sources, method_ids=method_ids, method_names=method_names, pmids=pubmeds)
            else:
                raise IncorrectNetworkFormat(self.network_format, self.formats)

        network_fd.close()

        print('Parsing network... finished!\n')

        return G


    def translate_network(self, translation_file, translation_id, translation_format, translated_network, translated_nodes):
        """
        Uses a BIANA translation file to return a new Network object with the 
        desired type of IDs and format

        @param:    translation_file
        @pdef:     File containing the translations from biana codes to the translation id
        @ptype:    String

        @param:    translation_id
        @pdef:     Type of ID in which the network will be translated
        @ptype:    String

        @param:    translation_format
        @pdef:     Format of the final network
        @ptype:    String

        @param:    translated_network
        @pdef:     File where the translated network will be written
        @ptype:    String

        @param:    translated_nodes
        @pdef:     File where the nodes file will be written
        @ptype:    String
        """

        if self.type_id != 'biana':
            raise IncorrectTypeID(self.type_id)

        NT.translate(self.network_file, self.node_file, translation_file, self.network_format, translation_format, translated_network, translated_nodes)

        return Network(translated_network, translated_nodes, translation_id, translation_format)


    def filter_network_by_tissue(self, filtered_network_file, filtered_nodes_file, tissue_object, permission):
        """
        Filter the network using only interactions where the proteins are 
        present in a Tissue object

        @param:    filtered_network_file
        @pdef:     File where the tissue-specific network will be written
        @ptype:    String

        @param:    filtered_nodes_file
        @pdef:     File where the nodes file will be written
        @ptype:    String

        @param:    tissue_object
        @pdef:     Tissue class object used to filter the network
        @ptype:    Tissue class object

        @param:    permission
        @pdef:     Level of permission to create the network
        @ptype:    Integer {0,1,2}
        """

        print('Filtering network by tissue...\n')
        TS.filter_network_tissue_specific(self.network_file, tissue_object, permission, filtered_network_file, filtered_nodes_file)
        print('Filtering network by tissue... finished!\n')

        return TissueSpecificNetwork(filtered_network_file, filtered_nodes_file, self.type_id, self.network_format, tissue_object, permission)


    def filter_network_by_method(self, methods_excluded=None, method_ids_excluded=None, methods_included=None, method_ids_included=None, output_network_file=None, output_nodes_file=None):
        """
        Filter the network: 
        - By excluding interactions only reported by methods in 'methods_excluded'
          and 'method_ids_excluded' list
        - By including interactions at least reported by one ofthe methods in 
          'methods_included' list and 'method_ids_included' list

        @param:    methods_excluded
        @pdef:     Method names list which will exclude interactions if they are
                   constituted by methods in this list
        @ptype:    List

        @param:    method_ids_excluded
        @pdef:     Same as methods_excluded for psi-mi IDs list
        @ptype:    List

        @param:    methods_included
        @pdef:     Method names list which will only include interactions that
                   at least contain one of the methods in this list
        @ptype:    List

        @param:    method_ids_included
        @pdef:     Same as methods_included for psi-mi IDs list
        @ptype:    List

        @param:    output_network_file
        @pdef:     File where the network will be written
        @ptype:    String

        @param:    output_nodes_file
        @pdef:     File where the nodes file will be written
        @ptype:    String

        @return:   Network object with the filtered network
        """
        if self.network_format != 'multi-fields':
            print('It is only possible to use this method with a multi-fields network\n')
            sys.exit(10)

        print('Filtering network by method...\n')

        fnet=nx.Graph()

        for u,v,d in self.network.edges_iter(data=True):
            skip_inc = False
            skip_exc = False

            # Check if at least one of the imprescindible methods is included
            if methods_included != None:
                skip_inc = True
                for method in d['method_names']:
                    if method in methods_included:
                        skip_inc = False
            if method_ids_included != None:
                skip_inc = True
                for method in d['method_ids']:
                    if method in methods_included:
                        skip_inc = False

            # Check if the interaction has at least another method apart from
            # the ones in methods_excluded/method_ids_excluded
            if methods_excluded != None:
                skip_exc = True
                for method in d['method_names']:
                    if method not in methods_excluded:
                        skip_exc = False

            if method_ids_excluded != None:
                method_ids_excluded = [str(x) for x in method_ids_excluded]
                skip_exc = True
                for method in d['method_ids']:
                    if method not in method_ids_excluded:
                        skip_exc = False

            if skip_inc == False and skip_exc == False:
                fnet.add_edge(u,v,d)

        self.write_network_file_from_networkx_graph(fnet, output_network_file, output_nodes_file, self.network_format, tissue_specific=self._tissue_specific)

        print('Filtering network by method... finished!\n')

        return Network(output_network_file, output_nodes_file, self.type_id, self.network_format)

    def filter_network_by_number_pubmeds(self, min_num_pubmeds, output_network_file, output_nodes_file):
        """
        Filter the network by minimum number of pubmed ids

        @param:    min_num_pubmeds
        @pdef:     Minimum number of pubmeds which an interaction has to have
        @ptype:    Integer

        @param:    output_network_file
        @pdef:     File where the network will be written
        @ptype:    String

        @param:    output_nodes_file
        @pdef:     File where the nodes file will be written
        @ptype:    String

        @return:   Network object with the filtered network
        """
        if self.network_format != 'multi-fields':
            print('It is only possible to use this method with a multi-fields network\n')
            sys.exit(10)

        print('Filtering network by number of pubmeds...\n')

        fnet=nx.Graph()

        for u,v,d in self.network.edges_iter(data=True):
            skip_inc = False
            skip_exc = False

            # Check if it has at least the minimum number of pubmeds required
            number_pubmeds = len(d['pmids'])
            if number_pubmeds >= min_num_pubmeds:
                fnet.add_edge(u,v,d)

        self.write_network_file_from_networkx_graph(fnet, output_network_file, output_nodes_file, self.network_format, tissue_specific=self._tissue_specific)

        print('Filtering network by number of pubmeds... finished!\n')

        return Network(output_network_file, output_nodes_file, self.type_id, self.network_format)

    def filter_network_by_database(self, databases_included, output_network_file, output_nodes_file):
        """
        Filter the network by interactions included only in certain databases.
        If the interaction is from one of the databases, it is included

        @param:    databases_included
        @pdef:     Databases from which the interaction must come from
        @ptype:    List

        @param:    output_network_file
        @pdef:     File where the network will be written
        @ptype:    String

        @param:    output_nodes_file
        @pdef:     File where the nodes file will be written
        @ptype:    String

        @return:   Network object with the filtered network
        """
        if self.network_format != 'multi-fields':
            print('It is only possible to use this method with a multi-fields network\n')
            sys.exit(10)

        print('Filtering network by databases...\n')

        fnet=nx.Graph()

        for u,v,d in self.network.edges_iter(data=True):
            skip_inc = False
            skip_exc = False

            # Check if at least one of the databases is included in the provided
            # list
            for database in d['sources']:
                if database in databases_included:
                    fnet.add_edge(u,v,d)
                    break

        self.write_network_file_from_networkx_graph(fnet, output_network_file, output_nodes_file, self.network_format, tissue_specific=self._tissue_specific)

        print('Filtering network by databases... finished!\n')

        return Network(output_network_file, output_nodes_file, self.type_id, self.network_format)

    def write_network_file_from_networkx_graph(self, input_network, output_network_file, output_nodes_file, output_network_format, tissue_specific=False):
        """
        Write a network file and a nodes file from a networkx graph object
        (Currently only available for multi-fields networks)
        """

        output_network_fd = open(output_network_file, 'w')

        for u,v,d in input_network.edges_iter(data=True):
            sources = ';'.join(d['sources'])
            method_ids = ';'.join(d['method_ids'])
            method_names = ';'.join(d['method_names'])
            pmids = ';'.join(d['pmids'])
            if output_network_format == 'multi-fields':
                if tissue_specific:
                    tissue_db = ';'.join(d['tissue_db'])
                    tissue_additional = ';'.join(d['tissue_additional'])
                    output_network_fd.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format( u,v,sources,method_ids,method_names,pmids,tissue_db,tissue_additional ))
                else:
                    output_network_fd.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format( u,v,sources,method_ids,method_names,pmids ))
            else:
                print('Not available format\n')
                sys.exit(10)

        output_network_fd.close()
        output_nodes_fd = open(output_nodes_file, 'w')

        for node in self.network.nodes():
            output_nodes_fd.write('{}\n'.format(node))

        output_nodes_fd.close()


class TissueSpecificNetwork(Network):
    """ 
    Child class of Network defining a tissue-specific network object 
    """

    def __init__(self, network_file, node_file, type_id, network_format, tissue_object, permission):
        """ 
        @param:    network_file
        @pdef:     Path to the file containing the edges of the network
        @ptype:    {String}

        @param:    node_file
        @pdef:     Path to the file containing the nodes of the network
        @ptype:    {String}

        @param:    type_id
        @pdef:     Type of IDs in the network
        @ptype:    {String}

        @param:    network_format
        @pdef:     Format of the network
        @ptype:    {String} {'sif' or 'multi-fields'}

        @param:    tissue_object
        @pdef:     Instance from Tissue Class associated to the Network
        @ptype:    {String}

        @param:    permission
        @pdef:     Level of permissivity of the network filtering
        @ptype:    {String}

        @raises: {IncorrectNetworkFormat} if the network_format is not in
        self.formats.
        @raises: {IncorrectTypeID} if the method translate_network is used with
        a network of type_id different from 'biana'
        """

        #super(TissueSpecificNetwork, self).__init__(network_file, node_file, type_id, network_format)

        self.network_file = network_file
        self.node_file = node_file
        self.type_id = type_id
        self.network_format = network_format
        self.formats = ['sif', 'multi-fields']
        self.executable_path = '/home/quim/project/tissue_specificity/scripts'

        self.tissue_object = tissue_object
        self.permission = permission

        self._tissue_specific = True 

        self.network = self.parse_network(tissue_specific=self._tissue_specific)

        self.hpa_edges = self.get_hpa_edges()
        self.jensen_edges = self.get_jensen_edges()

    ###########
    # METHODS #
    ###########

    def translate_network(self, translation_file, translation_id, translation_format, translated_network, translated_nodes):
        """
        Uses a BIANA translation file to return a new Network object with the 
        desired type of IDs and format
        """
        print('Translating network to {}...\n'.format(translation_id))

        if self.type_id != 'biana':
            raise IncorrectTypeID(self.type_id)

        NT.translate(self.network_file, self.node_file, translation_file, self.network_format, translation_format, translated_network, translated_nodes)

        print('Translating network to {}... finished!\n'.format(translation_id))

        return TissueSpecificNetwork(translated_network, translated_nodes, translation_id, translation_format, self.tissue_object, self.permission)

    def get_hpa_edges(self):
        """
        Obtain the edges in the tissue-specific network according to Human
        Protein Atlas 
        """
        hpa=nx.Graph()
        for u,v,d in self.network.edges_iter(data=True):
            if 'hpa' in d['tissue_db']:
                hpa.add_edge(u,v,d)
        return hpa.edges()

    def get_jensen_edges(self):
        """
        Obtain the edges in the tissue-specific network according to Tissues
        (Jensen Lab)
        """
        jensen=nx.Graph()
        for u,v,d in self.network.edges_iter(data=True):
            if 'jensen' in d['tissue_db']:
                jensen.add_edge(u,v,d)
        return jensen.edges()

    def get_union(self):
        """
        Obtain the common tissue-specific interactions between Human Protein
        Atlas and Tissues (Jensen Lab)
        """
        union = set(self.hpa_edges) & set(self.jensen_edges)

        return union

    def get_intersection(self):
        """
        Obtain the tissue-specific interactions in Human Protein Atlas, Tissues
        (Jensen Lab) or both
        """
        intersection = set(self.hpa_edges) | set(self.jensen_edges)

        return intersection



class Tissue(object):
    """ Class defining a tissue object """

    def __init__(self, tissue_terms_hpa, tissue_terms_jensen, jensen_conf=3, hpa_level='medium', hpa_rel='approved'):
        """ 
        @param:    tissue_terms_hpa
        @pdef:     Tissue terms of interest in Human Protein Atlas
        @ptype:    {List or String if only one term}

        @param:    tissue_terms_jensen
        @pdef:     Brenda Tissue Ontology names of interest in Tissues (Jensen Lab)
        @ptype:    {List or String if only one term}

        @param:    jensen_conf
        @pdef:     Level of confidence cut-off in Tissues (Jensen Lab) proteins
        @pdefault: 3
        @ptype:    {Integer or Float} {from 0 to 4}

        @param:    hpa_level
        @pdef:     Expression level cut-off in Human Protein Atlas proteins
        @pdefault: 'medium'
        @ptype:    {String} {'not detected','low','medium','high'}

        @param:    hpa_rel
        @pdef:     Reliability cut-off in Human Protein Atlas proteins
        @pdefault: 'approved'
        @ptype:    {String} {'uncertain','approved','supported'}
        """

        self.tissue_terms_hpa = self.check_tissue_terms(tissue_terms_hpa)
        self.tissue_terms_jensen = self.check_tissue_terms(tissue_terms_jensen)
        self.jensen_conf = jensen_conf
        self.hpa_level = hpa_level
        self.hpa_rel = hpa_rel

        # The scales of Human Protein Atlas levels and reliability:
        self.hpa_scale_level = {
        'not detected' : 0,
        'low' : 1,
        'medium' : 2,
        'high' : 3
        }
        self.hpa_scale_rel = {
        'uncertain' : 0,
        'approved' : 1,
        'supported' : 2
        }

        # Mapping files from tissue terms to biana user entities of tissues

        print('Generating tissue...\n')

        source_path = '/home/quim/project/tissue_specificity/scripts'
        BTOname_file = os.path.join(source_path,'BTOname2uE.pcl')
        HPA_tissue_file = os.path.join(source_path,'tissue2uEs.pcl')
        self.BTOname2uE = cPickle.load(open(BTOname_file))
        self.tissue2uEs = cPickle.load(open(HPA_tissue_file))

        self.user_entities_hpa = self.get_user_entities_hpa()
        self.user_entities_jensen = self.get_user_entities_jensen()
        self.all_tissue_user_entities = self.get_all_tissue_user_entities()

        # Mapping files from user entities of proteins to user entities of
        # tissues
        prot2tissues_file = os.path.join(source_path,'UEprot2UETissues.pcl')
        prot2HPA_file = os.path.join(source_path,'UEprot2UEHPA.pcl')
        self.UEprot2UETissues = cPickle.load(open(prot2tissues_file))
        self.UEprot2UEHPA = cPickle.load(open(prot2HPA_file))

        # Get all the tissue-specific proteins
        self.proteins_hpa = self.get_tissue_specific_proteins_hpa()
        self.proteins_jensen = self.get_tissue_specific_proteins_jensen()
        self.all_tissue_proteins = self.proteins_hpa | self.proteins_jensen

        print('Generating tissue... finished!\n')

    ###########
    # METHODS #
    ###########

    def check_tissue_terms(self, tissue_terms):
        """ If a string is introduced, it is added in a list object """
        if isinstance(tissue_terms, list):
            return tissue_terms
        elif isinstance(tissue_terms, str):
            termlist = []
            termlist.append(tissue_terms)
            return termlist
        else:
            print('Introduce a list in the tissue terms parameters!\n')
            raise ValueError

    def get_user_entities_hpa(self): 
        """ Returns user entities from Human Protein Atlas """
        self.user_entities_hpa = set([])
        for term in self.tissue_terms_hpa:
            if term in self.tissue2uEs:
                for uE in self.tissue2uEs[term]:
                    self.user_entities_hpa.add(uE)
        return self.user_entities_hpa

    def get_user_entities_jensen(self): 
        """ Returns user entities from Tissues (Jensen lab) """
        self.user_entities_jensen = set([])
        for term in self.tissue_terms_jensen:
            if term in self.BTOname2uE:
                self.user_entities_jensen.add(self.BTOname2uE[term])
        return self.user_entities_jensen

    def get_all_tissue_user_entities(self): 
        """ Returns user entities from Tissues (Jensen lab) """
        return self.user_entities_hpa | self.user_entities_jensen

    def get_tissue_specific_proteins_hpa(self): 
        """ 
        Returns all the user entities from Human Protein Atlas expressed in the
        tissues of interest above a certain level of expression and reliability
        """
        proteins_hpa = set()
        prot_hpa_to_values = {}
        for uE_protein in self.UEprot2UEHPA:
            for uE_tissue in self.UEprot2UEHPA[uE_protein]:
                if uE_tissue in self.user_entities_hpa:
                    level = self.UEprot2UEHPA[uE_protein][uE_tissue]['level']
                    reliability = self.UEprot2UEHPA[uE_protein][uE_tissue]['reliability']

                    # If our reliability is bigger or equal than the cut-off,
                    # we continue
                    if self.hpa_scale_rel[reliability] >= self.hpa_scale_rel[self.hpa_rel]:

                        # If our level of expression is bigger or equal than
                        # the cut-off, we continue
                        if self.hpa_scale_level[level] >= self.hpa_scale_level[self.hpa_level]:
                            proteins_hpa.add(uE_protein)
                            prot_hpa_to_values.setdefault(uE_protein, {})
                            prot_hpa_to_values[uE_protein].setdefault(uE_tissue, {})
                            prot_hpa_to_values[uE_protein][uE_tissue]['level'] = level
                            prot_hpa_to_values[uE_protein][uE_tissue]['reliability'] = reliability

        return proteins_hpa

    def get_tissue_specific_proteins_jensen(self): 
        """ 
        Returns all the user entities from Tissues (Jensen lab) expressed in the
        tissues of interest above a certain level of confidence
        """
        proteins_jensen = set()
        for uE_protein in self.UEprot2UETissues:
            for uE_tissue in self.UEprot2UETissues[uE_protein]:
                if uE_tissue in self.user_entities_jensen:
                    conf = self.UEprot2UETissues[uE_protein][uE_tissue]['confidence']
                    src = self.UEprot2UETissues[uE_protein][uE_tissue]['source']
                    evidence = self.UEprot2UETissues[uE_protein][uE_tissue]['evidence']

                    # If our confidence is bigger or equal than the cut-off, we
                    # continue
                    if float(conf) >= float(self.jensen_conf):
                        proteins_jensen.add(uE_protein)

        return proteins_jensen


class IncorrectNetworkFormat(NameError):
    """
    Subclass exception of the NameError which raises when a network format is
    not provided correctly
    """
    def __init__(self, network_format, formats):
        self.network_format = network_format
        self.formats = formats

    def __str__(self):
        return 'The network format {} is not valid. It must be one of the following formats: {}'.format(self.network_format, self.formats)

class IncorrectTypeID(Exception):
    """
    Exception that raises when a translation is done without having biana codes
    as initial type of IDs
    """
    def __init__(self, type_id):
        self.type_id = type_id

    def __str__(self):
        return 'The initial type of IDs of the network is not biana, it is {}. It is only possible to translate from BIANA codes'.format(self.type_id)


