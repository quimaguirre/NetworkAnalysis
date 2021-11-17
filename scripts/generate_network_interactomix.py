import argparse
import ConfigParser
import sys, os, re

import biana
try: from biana import *
except: sys.exit(10)

import methods_dictionaries as methods_dicts


def main():

    options = parse_user_arguments()
    generate_network(options)

def parse_user_arguments(*args, **kwds):

    parser = argparse.ArgumentParser(
        description = "Generate a protein-protein interaction network (implemented for Interactomix platform)",
        epilog      = "@oliva's lab 2019")
    parser.add_argument('-iseed','--seeds_input_file',dest='seed',action = 'store',
                        help = 'Seeds Input file (default is input_seed)')
    parser.add_argument('-radius','--radius_of_subnetwork_around_seeds',dest='radius',default=0,action = 'store',type=int,
                        help = '''Network is built in a radius of connections around the seed proteins.
                                  If 0, it creates the complete interactome''')
    parser.add_argument('-taxid','--TaxID',dest='taxid',action = 'store',default='9606',
                        help = 'Tax ID (i.e. human=9606 is default if TaxID=0 there is no restriction)')
    parser.add_argument('-stype','--seed_type',dest='stype',action = 'store',default='geneid',
                        help = 'Type of identifier for seeds (default is geneid)')
    parser.add_argument('-ttype','--translation_type',dest='ttype',action = 'store',default='accessionnumber',
                        help = '''Type of identifier for the output translation of codes (default is accessionnumber)
                        Using "proteinsequence" provides with the longest sequence of all codes''')
    parser.add_argument('-trans','--translation_of_nodes_file',dest='translation_file',action = 'store',default='translation_nodes.txt',
                        help = 'File with the translation of codes from BIANA to the selected type for all nodes')
    parser.add_argument('-strans','--translation_of_seeds_file',dest='translation_seeds_file',action = 'store',default='translation_seeds_to_BIANA_codes.txt',
                        help = 'File with the translation of codes from the introduced type of code to BIANA codes')
    parser.add_argument('-edge','--edge_file',dest='edge',action = 'store', default='biana_edges',
                        help = 'Output file with edges(default is biana_edges)')
    parser.add_argument('-node','--node_file',dest='node',action = 'store', default='biana_nodes',
                        help = 'Output file with nodes(default is biana_nodes)')
    parser.add_argument('-format','--output_format',dest='format',action = 'store',default='sif',
                        help = '''Format file of the edge file:\tsif (default), netscore, raw, multi-fields:\n
                                    'sif': <node1>\tscore\t<node2>\n
                                    'netscore': <node1>\t<node2>\t<score>\n
                                    'raw': <node1>\t<node2>\n
                                    'multi-fields' : <node1>\t<node2>\t<sources>\t<method_ids>\t<method_names>\t<pmids>\n''')
    parser.add_argument('-rAFF','--restricted_to_TAP',dest='restricted_to_TAP',action = 'store_true',
                        help = 'Flag to use interactions at least described by affinity methods (i.e. Tandem Affinity Purification)')
    parser.add_argument('-rY2H','--restricted_to_Y2H',dest='restricted_to_Y2H',action = 'store_true',
                        help = 'Flag to use interactions at least described by yeast two hybrid methods (Y2H)')
    parser.add_argument('-rUSER','--restricted_to_user',dest='restricted_to_user',action = 'store',default='restricted_methods',
                        help = 'File to use interactions described by the user selected methods')
    parser.add_argument('-eAFF','--except_TAP',dest='except_TAP',action = 'store_true',
                        help = 'Flag to use all interactions except those described by affinity methods (i.e. Tandem Affinity Purification)')
    parser.add_argument('-eY2H','--except_Y2H',dest='except_Y2H',action = 'store_true',
                        help = 'Flag to use all interactions except those described by yeast two hybrid methods (Y2H)')
    parser.add_argument('-eUSER','--except_user',dest='except_user',action = 'store',default='restricted_methods',
                        help = 'File to reject interactions described by the user selected methods')
    parser.add_argument('-v','--verbose',dest='verbose',action = 'store_true',
                        help = 'Flag to use verbose mode')
    options=parser.parse_args()

    """
    Example:
    python generate_network_interactomix.py -iseed example/sample1.txt -radius 1 -taxid 9606 -stype uniprotentry -ttype proteinsequence -trans example/output/example.proteinsequence.trans -strans example/output/example.seeds.trans -edge example/output/example.edges -node example/output/example.nodes -format raw -rY2H
    python /home/quim/PHD/Projects/BIANA/scripts/generate_network_interactomix.py -radius 0 -taxid 9606 -edge /home/quim/PHD/Projects/BIANA/outputs/BIANA_2020_geneID_seqtax_drugtarget/human_network_biana_2020.txt -node /home/quim/PHD/Projects/BIANA/outputs/BIANA_2020_geneID_seqtax_drugtarget/human_network_biana_2020_nodes.txt -trans /home/quim/PHD/Projects/BIANA/outputs/BIANA_2020_geneID_seqtax_drugtarget/human_network_biana_2020_translation.txt -ttype geneid -format multi-fields &> /home/quim/PHD/Projects/BIANA/outputs/BIANA_2020_geneID_seqtax_drugtarget/human_network_biana_2020.log
    """

    return options


def generate_network(options):
    """
    Generates a protein-protein interaction network extracting information from BIANA.
    """

    #----------------------#
    #   FIXED PARAMETERS   #
    #----------------------#

    # Parameters that I have decided to fix
    restricted_to_seeds = False
    minimum_number_of_methods = 1
    minimum_number_of_db = 1
    seed_score = 0.1


    #--------------------------------------#
    #   GET INFORMATION FROM CONFIG FILE   #
    #--------------------------------------#

    # Get the program path
    main_path = os.path.abspath(os.path.dirname(__file__))

    # Read the config file
    config_file = os.path.join(main_path, 'config.ini')
    config = ConfigParser.ConfigParser()
    config.read(config_file)


    #--------------------------------------#
    #   LOAD THE DICTIONARIES OF METHODS   #
    #--------------------------------------#

    # Get the affinity dictionary
    affinity_dict = methods_dicts.affinity_dict
    affinity=set(affinity_dict.keys())

    # Get the complementation dictionary
    complementation_dict = methods_dicts.complementation_dict
    complementation=set(complementation_dict.keys())


    #---------------------------------------#
    #   GET METHODS THAT WILL BE FILTERED   #
    #---------------------------------------#

    # Check if the user has introduced a file with methods that must be included
    if not fileExist(options.restricted_to_user):
        print "No restriction on methods selected by the user"
        user_selection=False
    else:
        use_methods=[]
        with open(options.restricted_to_user) as input_method_fd:
            for line in input_method_fd:
                fields = line.strip().split("\t")
                use_methods.append(fields[0])
        user_selection=True
        print "Input to use only Methods:",repr(use_methods)

    # Check if the user has introduced a file with methods that have to be excluded
    if not fileExist(options.except_user):
        print "No rejection of methods selected by the user"
        user_rejection=False
    else:
        no_methods=[]
        with open(options.except_user) as input_method_fd:
            for line in input_method_fd:
                fields = line.strip().split("\t")
                no_methods.append(fields[0])
        user_rejection=True
        print "Input of rejected Methods:",repr(no_methods)


    #---------------------------#
    #   START A BIANA SESSION   #
    #---------------------------#

    print "Open session"

    session = create_new_session( sessionID="biana_session",
                                  dbname=config.get('BIANA', 'database'),
                                  dbhost=config.get('BIANA', 'host'),
                                  dbuser=config.get('BIANA', 'user'),
                                  dbpassword=config.get('BIANA', 'password'),
                                  unification_protocol=config.get('BIANA', 'unification_protocol') )
    print "Continue"


    #------------------------------#
    #   DEFINE A USER ENTITY SET   #
    #------------------------------#

    # Create network network of expansion if the radius is larger than 0
    if restricted_to_seeds or options.radius>0:

        # Check if the seeds file exists
        if not fileExist(options.seed):
            print "File with seeds is missing or not found"
            sys.exit(10)
        else:
            level=options.radius
            seed_list = get_seeds_from_file(options.seed)

            # If we have Taxonomy restriction, we add it
            if options.taxid != "0":
                print("Check Proteome %s"%(repr(options.taxid)))
                proteome = session.create_new_user_entity_set(  identifier_description_list =seed_list,
                                                                attribute_restriction_list=[("taxid",options.taxid)],
                                                                id_type=options.stype,new_user_entity_set_id="proteome",
                                                                negative_attribute_restriction_list=[] )
            else:
                print('Proteome without Taxonomy restriction')
                proteome = session.create_new_user_entity_set(  identifier_description_list =seed_list,
                                                                id_type=options.stype,new_user_entity_set_id="proteome",
                                                                negative_attribute_restriction_list=[] )

    else:
        level=0
        proteome = session.create_new_user_entity_set( identifier_description_list = [("taxid",options.taxid)],
                                                      attribute_restriction_list=[], id_type="embedded",
                                                      new_user_entity_set_id="proteome",
                                                      negative_attribute_restriction_list=[] )


    #----------------------------------------------------#
    #   SELECT THE INTERACTIONS OF THE USER ENTITY SET   #
    #----------------------------------------------------#

    print ("Selecting interactions")

    # Select interactions that have been detected at least by affinity technology
    if options.restricted_to_TAP:
        print ('Using interactions at least described by affinity methods (i.e. Tandem Affinity Purification)')
        session.create_network( user_entity_set_id = "proteome" , level = level, relation_type_list=["interaction"] ,
                               relation_attribute_restriction_list = [("Method_id",400)],
                                #relation_attribute_restriction_list  = [("psimi_name","affinity technology")],
                               include_relations_last_level = (not restricted_to_seeds) , use_self_relations = False)

    # Select interactions that have been detected at least by yeast two hybrid
    elif options.restricted_to_Y2H:
        print ('Using interactions at least described by yeast-two-hybrid methods (Y2H)')
        session.create_network( user_entity_set_id = "proteome" , level = level, relation_type_list=["interaction"] ,
                               relation_attribute_restriction_list = [("Method_id",18)],
                               #relation_attribute_restriction_list = [("psimi_name","y2h2")],
                               include_relations_last_level = (not restricted_to_seeds) , use_self_relations = False)

    # Select all interactions
    else:
        session.create_network( user_entity_set_id = "proteome" , level = level, relation_type_list=["interaction"] ,
                               include_relations_last_level = (not restricted_to_seeds) , use_self_relations = False)



    # Summary of interactions
    out_network = open(options.edge,'w')
    all_interactions = proteome.getRelations()
    print "Num interactions:", len(all_interactions)


    #--------------------------------------#
    #   FILTER THE SELECTED INTERACTIONS   #
    #--------------------------------------#

    nodes=set()

    # Get all the user entity ids from the user entity set 'proteome'
    all_uEs = proteome.get_user_entity_ids()
    # Obtain a dictionary user entity ID => type
    uEId_to_type = session.dbAccess.get_user_entity_type(config.get('BIANA', 'unification_protocol'), all_uEs)

    skip_interactions=0
    for (uE_id1, uE_id2) in all_interactions:

	#self.dbAccess.get_external_entities_dict( externalEntityIdsList = [external_entity_relation_id] )

        # Get TYPE of user entity
        uE1_type = uEId_to_type[uE_id1]
        uE2_type = uEId_to_type[uE_id2]
        # If type is not protein, we skip the interaction
        if uE1_type != 'protein' or uE2_type != 'protein':
            if options.verbose:
                print('Skipping interaction because the type of one of the user entities is not protein!')
                print('Node 1: {}\tType: {}'.format(uE_id1, uE1_type))
                print('Node 2: {}\tType: {}'.format(uE_id2, uE2_type))
            skip_interactions=skip_interactions+1
            continue

        eErIDs_list = proteome.get_external_entity_relation_ids(uE_id1, uE_id2)

        method_names = set()
        method_ids = set()
        source_databases = set()
        use_method_ids=set()
        pubmed_ids = set()
        unused_method_names = set()


        relationObj_dict = session.dbAccess.get_external_entities_dict(
                                externalEntityIdsList = eErIDs_list, attribute_list = [],
                                relation_attribute_list = ["method_id","psimi_name","pubmed"], participant_attribute_list = [] )

        num_methods=0
        for current_eErID in eErIDs_list:
            relationObj = relationObj_dict[current_eErID]
            if options.verbose:
                print "Interaction: (",uE_id1,",",uE_id2,")"
                print relationObj

            #if relationObj.get_attribute(attribute_identifier="psimi_name") is not None:
            #    print "\t".join([ x.value for x in relationObj.get_attribute(attribute_identifier="psimi_name") ])
            #if relationObj.get_attribute(attribute_identifier="method_id") is not None:
            #print "\t".join([ x.value for x in relationObj.get_attribute(attribute_identifier="method_id") ])
            #print relationObj.get_attributes_dict()
            #print [ x.value for x in relationObj.get_attributes_dict()["psimi_name"] ]
            #print [ x.value for x in relationObj.get_attributes_dict()["method_id"] ]

            if "psimi_name" in relationObj.get_attributes_dict():
                method_names.update([ str(x.value) for x in relationObj.get_attributes_dict()["psimi_name"] ])
            if "method_id" in relationObj.get_attributes_dict():
                method_ids.update([ x.value for x in relationObj.get_attributes_dict()["method_id"]])
            if "pubmed" in relationObj.get_attributes_dict():
                pubmed_ids.update([ x.value for x in relationObj.get_attributes_dict()["pubmed"]])
            source_databases.add(str(session.dbAccess.get_external_database(
                                    database_id = relationObj.get_source_database()) ))
            if options.except_TAP:
                for m in method_ids:
                    if m not in affinity:
                        use_method_ids.add(m)
                        #print "Add", m
                    else:
                        unused_method_names.add(affinity_dict[m])
            elif options.except_Y2H:
                #print "check Y2H"
                for m in method_ids:
                    if m not in complementation:
                        use_method_ids.add(m)
                        #print "Add", m
                    else:
                        unused_method_names.add(complementation_dict[m])
            elif user_rejection:
                for m in method_ids:
                    if m not in no_methods:
                        use_method_ids.add(m)
            elif user_selection:
                for m in method_ids:
                    #print "Check",repr(use_methods)
                    if m in set(use_methods):
                        use_method_ids.add(m)
                    if options.verbose:
                        print "Not among selected methods ",m
            else:
                use_method_ids.update(method_ids)

        if len(source_databases) > 0:
            info_sources=";".join([str(x) for x in source_databases])
        else:
            if options.verbose:
                print('Skipping interaction it has no source database!')
                print('Node 1: {}\tNode 2: {}'.format(uE_id1, uE_id2))
            skip_interactions=skip_interactions+1
            continue
        if len(method_names) > 0:
            method_names = [x for x in method_names if x not in unused_method_names] # Remove method names that were excluded
            info_methods=";".join([str(x) for x in method_names])
        else:
            info_methods='-'
        if len(use_method_ids) > 0:
            info_methods_ids=";".join([str(x) for x in use_method_ids])
        else:
            if options.verbose:
                print('Skipping interaction it has no method!')
                print('Node 1: {}\tNode 2: {}'.format(uE_id1, uE_id2))
            skip_interactions=skip_interactions+1
            continue
        if len(pubmed_ids) > 0:
            info_pubmed_ids=";".join([str(x) for x in pubmed_ids])
        else:
            info_pubmed_ids='-'
        num_databases=len(source_databases)
        num_methods=len(use_method_ids)
        num_pubmeds = len(pubmed_ids)

        if options.verbose:
            print "Methods",num_methods,info_methods,"\tSelected:",info_methods_ids
            print "Databases",num_databases,info_sources
            print "Pubmeds",num_pubmeds,info_pubmed_ids

        # Check if the number of methods is higher than the minimum established
        if num_methods >= minimum_number_of_methods:
            use=True
        else:
            use=False

        # Check if the number of database is higher than the minimum established
        if use and num_databases >= minimum_number_of_db:
            use=True
        else:
            use=False

        if not use:
            skip_interactions=skip_interactions+1
        #print method_names, method_ids, source_databases


        #----------------------#
        #   OUTPUT EDGE FILE   #
        #----------------------#

        if use:
            #print uE_id1, uE_id/2
            nodes.add(uE_id1)
            nodes.add(uE_id2)
            #print "Attribute ",(uE_id1,uE_id2).get_attribute(

            if options.format == 'multi-fields' :
                out_network.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".
                              format(uE_id1,uE_id2,info_sources,info_methods_ids,info_methods,info_pubmed_ids))
            elif options.format == 'netscore':
                out_network.write('\t{}\t{}\t{:.2f}\n'.format(uE_id1,uE_id2,1.))
            elif options.format == 'raw':
                out_network.write("{}\t{}\n".format(uE_id1,uE_id2))
            else:
                # If the format is not multi-fields, netscore or raw, the output format is sif
                out_network.write("{}\t{:.2f}\t{}\n".format(uE_id1,1.,uE_id2))


    print "Num neglected interactions:", skip_interactions
    out_network.close()


    #---------------------------------------#
    #   OUTPUT NODE AND TRANSLATION FILES   #
    #---------------------------------------#

    # If we wanted the complete interactome, the translation will be done differently

    if options.radius <= 0:

        # Output node file

        out_proteins = open(options.node,'w')
        for protein in nodes:
            if options.format == 'multi-fields':
                out_proteins.write("{0}\t{1:.2f}\t{2:.2f}\t{3:.2f}\n".format(protein,1.,1.,0.1))
            elif options.format == 'netscore':
                out_proteins.write("{0}\t{1:.2f}\t{2:.2f}\t{3:.2f}\n".format(protein,1.,1.,0.1))
            else:
                out_proteins.write("{0}\t{1:.2f}\n".format(protein,0.1))
        out_proteins.close()


        ################################# TRANSLATION ####################################
        out_translation = open(options.translation_file,'w')

        # TRANSLATION TO 'stype'
        trans_stype=False
        if options.stype != 'proteinsequence' and options.stype != options.ttype:
            trans_stype = True
            out_trans_stype = open(options.translation_file+'.'+options.stype+'.trans','w')

        for protein in nodes:
            uE = session.get_user_entity(protein)
            translate=set()
            translate_stype=set()
            if options.ttype == "proteinsequence":
                maxlen=0;
                for current_id in uE.get_attribute(attribute_identifier=options.ttype):
                    if maxlen < len(current_id.value.get_sequence().upper()):
                        maxlen=len(current_id.value.get_sequence().upper())
                translation=",".join([str(current_id.value.get_sequence().upper()) for current_id in uE.get_attribute(attribute_identifier=options.ttype) if len(str(current_id.value.get_sequence().upper())) == maxlen ] )
                #print "Translation",protein,translation
                #print("{0}\t'{1}'\n".format(protein,translation))
            else:
                ##### TRANSLATION TO 'ttype'
                for current_id in uE.get_attribute(attribute_identifier=options.ttype):
                    translate.add(current_id.value.upper())
                translation="','".join(["{0}".format(x) for x in translate])
            out_translation.write("{0}\t'{1}'\n".format(protein,translation))
            ##### TRANSLATION TO STYPE
            if trans_stype:
                for current_id in uE.get_attribute(attribute_identifier=options.stype):
                    translate_stype.add(current_id.value.upper())
                translation_stype="','".join(["{0}".format(x) for x in translate_stype])
                out_trans_stype.write("{0}\t'{1}'\n".format(protein,translation_stype))
        out_translation.close()
        if trans_stype:
            out_trans_stype.close()
        ####################################################################################


    # If we wanted a network of expansion, the translation will be done differently
    elif options.radius > 0:

        # Read the seeds

        seeds=set()
        input_seed = open(options.seed,'r')
        for line in input_seed:
            fields = line.strip().split("\t")
            seeds.add(fields[0].lower())
        input_seed.close()

        # Output node file

        out_proteins = open(options.node,'w')
        translate={}
        for protein in nodes:
            score=seed_score
            uE = session.get_user_entity(protein)
            for current_id in uE.get_attribute(attribute_identifier=options.stype):
                if current_id.value.lower() in seeds:
                    translate.setdefault(current_id.value.lower(),[])
                    translate[current_id.value.lower()].append(protein)
                    score=1.0
            if options.format == 'multi-fields':
                out_proteins.write("{0}\t{1:.2f}\t{2:.2f}\t{3:.2f}\n".format(protein,1.,1.,score))
            elif options.format == 'netscore':
                out_proteins.write("{0}\t{1:.2f}\t{2:.2f}\t{3:.2f}\n".format(protein,1.,1.,score))
            else:
                out_proteins.write("{0}\t{1:.2f}\n".format(protein,score))
        out_proteins.close()


        # Get the IDS of single nodes that were not previously found in the network

        single=set()
        for uE_id in proteome.get_unconnected_nodes():
            single.add(uE_id)
        for protein in single:
            uE = session.get_user_entity(protein)
            for current_id in uE.get_attribute(attribute_identifier=options.stype):
                if current_id.value.lower() in seeds:
                    translate.setdefault(current_id.value.lower(),[])
                    translate[current_id.value.lower()].append(protein)

        
        # Get all IDS of SEEDS, defined as "proteome", and check missing codes to be
        # added for translation

        allseed=set()
        for uE_id in proteome.get_user_entity_ids():
            allseed.add(uE_id)
        for protein in allseed:
            if protein not in single and protein not in nodes:
                uE = session.get_user_entity(protein)
                for current_id in uE.get_attribute(attribute_identifier=options.stype):
                    if current_id.value.lower() in seeds:
                        translate.setdefault(current_id.value.lower(),[])
                        translate[current_id.value.lower()].append(protein)


        ################################# TRANSLATION ####################################
        out_translation = open(options.translation_seeds_file,'w')
        for s in seeds:
            if s == '': continue
            if s in translate:
                codes=set(translate[s])
                translation="','".join([str(x) for x in codes])
                #out_translation.write("%s\t'%s'\n" % (s.upper(),translation))
                out_translation.write("{0}\t'{1}'\n".format(s.upper(),translation))
            else:
                out_translation.write("{0}\t'Unknown'\n".format(s.upper()))
        out_translation.close()

        # Output translation file

        # TRANSLATION TO 'ttype'
        out_translation = open(options.translation_file,'w')

        # TRANSLATION TO 'stype'
        trans_stype=False
        if options.stype != 'proteinsequence' and options.stype != options.ttype:
            trans_stype = True
            out_trans_stype = open(options.translation_file+'.'+options.stype+'.trans','w')

        for protein in nodes:
            uE = session.get_user_entity(protein)
            translate=set()
            translate_stype=set()
            if options.ttype == "proteinsequence":
                maxlen=0;
                for current_id in uE.get_attribute(attribute_identifier=options.ttype):
                    if maxlen < len(current_id.value.get_sequence().upper()):
                        maxlen=len(current_id.value.get_sequence().upper())
                translation=",".join([str(current_id.value.get_sequence().upper()) for current_id in uE.get_attribute(attribute_identifier=options.ttype) if len(str(current_id.value.get_sequence().upper())) == maxlen ] )
                #print "Translation",protein,translation
                #print("{0}\t'{1}'\n".format(protein,translation))
            else:
                for current_id in uE.get_attribute(attribute_identifier=options.ttype):
                    translate.add(current_id.value.upper())
                translation="','".join(["{0}".format(x) for x in translate])
            out_translation.write("{0}\t'{1}'\n".format(protein,translation))

            ##### TRANSLATION TO STYPE
            if trans_stype:
                for current_id in uE.get_attribute(attribute_identifier=options.stype):
                    translate_stype.add(current_id.value.upper())
                translation_stype="','".join(["{0}".format(x) for x in translate_stype])
                out_trans_stype.write("{0}\t'{1}'\n".format(protein,translation_stype))            

        out_translation.close()
        if trans_stype:
            out_trans_stype.close()
        ####################################################################################


    print('Generation of the network done!')

    return



def fileExist(file):
    """
    Checks if a file exists AND is a file
    """
    return os.path.exists(file) and os.path.isfile(file)


def get_seeds_from_file(seed_file):
    """
    Obtain the seeds from a file and introduce them to a Python list.
    The seeds must be separated by new lines!
    """
    seed_set = set()
    with open(seed_file, 'r') as seed_file_fd:
        for line in seed_file_fd:
            fields = line.strip().split('\t')
            seed_set.add(fields[0])
    return list(seed_set)



if  __name__ == "__main__":
    main()



