from biana import *
import sets

from biana.utilities import identifier_utilities



taxid = 9606
#taxid = 10116
level = 0
restricted_to_seeds = False
output_prefix = 'human_reactome_may18'
#output_prefix = 'rat_reactome_may18'
output_path = '/home/quim/data/networks/human_reactome_may18'
#output_path = '/home/quim/data/networks/rat_reactome_may18'

# START BIANA SESSION

biana_session_object = create_new_session(sessionID='session_ID', 
                dbname='test_BIANA_MAY_2018',
                dbhost='localhost',
                dbuser='quim',
                dbpassword='', 
                unification_protocol='geneID_seqtax_drugtarget')


# CREATE A LIST WITH ALL THE SEED IDENTIFIERS
#list_input_identifiers = [("uniprotentry","ACE_YEAST"),("uniprotentry","PGH2_HUMAN"),("uniprotentry","RND3_HUMAN"),("uniprotentry","ACE_YEAST"),("uniprotentry","PGH2_HUMAN"),("uniprotentry","RND3_HUMAN")]
list_input_identifiers = [("taxid",taxid)]
list_input_restriction_identifiers = []
list_input_negative_restriction_identifiers = []


# CREATE THE SET

user_entity_set_object = biana_session_object.create_new_user_entity_set( identifier_description_list = list_input_identifiers, 
                            attribute_restriction_list = list_input_restriction_identifiers,
                            negative_attribute_restriction_list = list_input_negative_restriction_identifiers,                  
                            id_type='embedded', 
                            new_user_entity_set_id="my_user_entity_set")


# CREATE NETWORK COMMANDS

biana_session_object.create_network(user_entity_set_id = 'my_user_entity_set', 
                    level = level, 
                    include_relations_last_level = (not restricted_to_seeds), 
                    relation_type_list = ["reaction"], 
                    relation_attribute_restriction_list = [],
                    use_self_relations = False,
                    expansion_attribute_list = [],
                    expansion_relation_type_list = [], 
                    attribute_network_attribute_list = [], 
                    group_relation_type_list = [])


# OUTPUT COMMANDS

biana_session_object.output_user_entity_set_network_in_sif_format(user_entity_set_id='my_user_entity_set', 
                    output_path = output_path,
                    output_prefix = output_prefix, 
                    node_attributes = ["uniprotaccession", "genesymbol"], 
                    participant_attributes = ["role"], 
                    relation_attributes = ['psimi_name', 'Pubmed'], 
                    include_tags = True,
                    output_1_value_per_attribute = True,
                    only_selected = False)
