import argparse
import ConfigParser
import time
import sys, os, re

from biana.biana_commands import available_sessions, create_new_session
from biana.utilities import biana_output_converter
from context import NetworkAnalysis




def main():

    create_network_in_GUILDify_style()


#################
#################
# MAIN FUNCTION #
#################
#################

def create_network_in_GUILDify_style():
    """
    Generates the profiles of the input drug
    """

    # Start marker for time measure
    start = time.time()

    # Get the program path
    main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

    # Read the config file
    config_file = os.path.join(main_path, 'config.ini')
    config = ConfigParser.ConfigParser()
    config.read(config_file)

    tax_id = '9606'
    output_path = os.path.abspath(os.path.join(os.path.join(os.path.dirname(__file__), '..'), 'workspace/GUILDify_MAY_2018'))
    create_directory(output_path)
    network_file_prefix = 'network'
    node_info_file = os.path.join(output_path, "node_info.txt")
    filtered_network_file = os.path.join(output_path, "network_filtered.sif")

    print "Open session"
    biana_session = create_new_session( sessionID="biana_session",
                                        dbname=config.get('BIANA', 'database'),
                                        dbhost=config.get('BIANA', 'host'),
                                        dbuser=config.get('BIANA', 'user'),
                                        dbpassword=config.get('BIANA', 'password'),
                                        unification_protocol=config.get('BIANA', 'unification_protocol') )
    print "Continue"

    unification_protocol = config.get('BIANA', 'unification_protocol')
    NODE_ATTRIBUTES = [ "uniprotaccession", "genesymbol", "name", "geneid" ]
    if not fileExist(filtered_network_file):
        user_entity_ids_all = create_interactome_using_biana(biana_session, tax_id, output_path, network_file_prefix)
    else:
        # Get user entities from node info file
        user_entity_ids_all = set()
        with open(node_info_file) as node_fd:
            first = node_fd.readline()
            for line in node_fd:
                fields = line.strip().split('\t')
                user_entity_ids_all.add(fields[0])
    ueids = output_user_entities_info(biana_session, unification_protocol, user_entity_ids_all, NODE_ATTRIBUTES, node_info_file)
    #create_filtered_network_file(os.path.join(output_path, network_file_prefix), filtered_network_file, ueids)

    # End marker for time
    end = time.time()
    print('\nTIME OF EXECUTION: {:.3f} seconds or {:.3f} minutes.\n'.format(end - start, (end - start) / 60))

    return

#######################
#######################
# SECONDARY FUNCTIONS #
#######################
#######################

def remove_isoforms(entry_ids):
    """
    Remove the isoforms from a list of uniprot accessions
    """
    # Check how many non-isoforms are there
    non_isoforms = [ entry for entry in entry_ids.split('; ') if len(entry.split('-')) == 1 ]

    # Remove the isoforms
    old_entry_ids = entry_ids.split('; ')
    new_entry_ids = set()
    for entry in old_entry_ids:
        if len(non_isoforms) > 0:
            if len(entry.split('-')) == 1:
                new_entry_ids.add(entry)
        else:
            new_entry_ids.add(entry)
    entry_ids = '; '.join(list(new_entry_ids))
    return entry_ids


def create_interactome_using_biana(biana_session, tax_id, output_path, network_file_prefix):
    """
        Creates ppi network files using BIANA
        network_file_prefix: network_file_prefix.sif (network in sif format) also output network_file_attribute.eda (Cytoscape edge attribute files for each attribute)
    """
    identifier_description_list = [("taxid", tax_id)]
    relation_type_list = []
    relation_attribute_restriction_list = []
    relation_attributes=["pubmed", "method_id"] #, "stringscore"]
    relation_type_list = ["interaction"] #, "functional_association"] # ["interaction","reaction","complex","pathway"] 
    # Create a new User Entity Set (group of biomolecules) in this session
    user_entity_set = biana_session.create_new_user_entity_set(identifier_description_list=identifier_description_list, attribute_restriction_list=[], id_type="embedded", new_user_entity_set_id="User_Entity_Set_1")
    # Fetch relations of biomolecules in the set 
    #biana_session.create_network(user_entity_set_id = "User_Entity_Set_1" , level = 0, relation_type_list=relation_type_list, relation_attribute_restriction_list=relation_attribute_restriction_list, include_relations_last_level=True, use_self_relations=False)
    # Export all the information of this set and its network to a file in a tab separated format
    #biana_session.output_user_entity_set_network_in_sif_format(user_entity_set_id = "User_Entity_Set_1", output_path = output_path, output_prefix = network_file_prefix, node_attributes = [], participant_attributes = [], relation_attributes = relation_attributes, output_1_value_per_attribute = False, include_tags = False)
    print('Network finished')
    return map(lambda x: str(int(x)), user_entity_set.get_user_entity_ids())

def output_user_entities_info(biana_session, unification_protocol, user_entity_ids, node_attributes, file_name):
    ueids = set()
    f = open(file_name, 'w')
    f.write("BIANA ID\tUniProt ID\t%s\tEquivalent Entries\n" % "\t".join(node_attributes[1:]))
    for user_entity_id in user_entity_ids:
        attribute_values_uniprot_all = None
        attribute_values_kegg_all = None
        skip_flag = False
        inner_values = [ user_entity_id ]
        for attribute in node_attributes:
            if attribute == "uniprotaccession":
                attribute_values = biana_session.dbAccess.get_user_entity_attributes( unification_protocol_name = unification_protocol, listUserEntityID = [user_entity_id], attribute_identifier = attribute, only_uniques = True).values()
                if len(attribute_values) > 0: attribute_values = set(attribute_values[0])
                else: attribute_values = set()
                attribute_values_uniprot_all = biana_session.dbAccess.get_user_entity_attributes( unification_protocol_name = unification_protocol, listUserEntityID = [user_entity_id], attribute_identifier = attribute, only_uniques = False).values()
                if len(attribute_values_uniprot_all) > 0: attribute_values_uniprot_all = set(attribute_values_uniprot_all[0])
                else: attribute_values_uniprot_all = set()
                attribute_values_uniprot_all -= attribute_values

                temp_str = "; ".join([ x for x in attribute_values])  
                temp_str = remove_isoforms(temp_str)

                # attribute_values = biana_session.dbAccess.get_user_entity_attributes( unification_protocol_name = unification_protocol, listUserEntityID = [user_entity_id], attribute_identifier = "kegggene", only_uniques = True).values()
                # if len(attribute_values) > 0: attribute_values = set(attribute_values[0])
                # else: attribute_values = set()
                # attribute_values_kegg_all = biana_session.dbAccess.get_user_entity_attributes( unification_protocol_name = unification_protocol, listUserEntityID = [user_entity_id], attribute_identifier = "kegggene", only_uniques = False).values()
                # if len(attribute_values_kegg_all) > 0: attribute_values_kegg_all = set(attribute_values_kegg_all[0])
                # else: attribute_values_kegg_all = set()
                # attribute_values_kegg_all_new = set()
                # for attribute_value in attribute_values_kegg_all:
                #     index = attribute_value.find(":")
                #     if index != -1:
                #         attribute_values_kegg_all_new.add(attribute_value[index:])
                # attribute_values_kegg_all = attribute_values_kegg_all_new
                # attribute_values_kegg_all -= attribute_values
                # temp_str2 = "; ".join([ "KEGG:%s" % x for x in attribute_values])  
                # if temp_str2 != "": 
                #     if temp_str != "":
                #         temp_str += "; " 
                #     temp_str += temp_str2

                if temp_str == "": 
                    temp_str = "-"
                    skip_flag = True
                inner_values.append( temp_str )
            else:
                DESCRIPTION_EXP = re.compile("(.*?);")
                if attribute == "description":
                    attribute_values = biana_session.dbAccess.get_user_entity_attributes( unification_protocol_name = unification_protocol, listUserEntityID = [user_entity_id], attribute_identifier = attribute, only_uniques = False).values()
                else:
                    attribute_values = biana_session.dbAccess.get_user_entity_attributes( unification_protocol_name = unification_protocol, listUserEntityID = [user_entity_id], attribute_identifier = attribute, only_uniques = True).values()
                if len(attribute_values) > 0: attribute_values = attribute_values[0]
                if attribute == "name":
                    if len(attribute_values) > 1:
                        attribute_values = [max(attribute_values, key=len)] # If multiple names, get the longest one
                if attribute == "description":
                    attribute_values_new = []
                    for value in attribute_values:
                        frases = value.split("RecName: Full=")
                        for frase in frases:
                            if frase == "":
                                continue
                            m = re.search(DESCRIPTION_EXP, frase)
                            if m:
                                attribute_values_new.append(m.group(0).rstrip(";"))
                    attribute_values = attribute_values_new

                temp_str = "; ".join([ x.replace("\n"," ") for x in attribute_values])  
                if temp_str == "": temp_str = "-"

                # Skip the entries without GeneID
                if attribute.lower() == "geneid" and temp_str == "-":
                    skip_flag = True

                inner_values.append( temp_str )

        temp_str = "; ".join([ x for x in attribute_values_uniprot_all]) 

        # temp_str2 = "; ".join([ "KEGG:%s" % x for x in attribute_values_kegg_all])
        # if temp_str2 != "":
        #     if temp_str != "":
        #         temp_str += "; " 
        #     temp_str += temp_str2

        if temp_str == "": temp_str = "-"
        inner_values.append(temp_str)
        if skip_flag:
            continue
        f.write("%s\n" % "\t".join(inner_values))
        ueids.add(user_entity_id)
    f.close()
    return ueids

def create_filtered_network_file(network_file_prefix, filtered_network_file, ueids):
    """
        Creates a new sif file from the given network file where interactions coming from non-tap experiments are filtered. Then filters interactions that are not in given user_entity_id set.
    """
    network_file_method_attribute = network_file_prefix + "_method_id.eda"
    network_file_source_attribute = network_file_prefix + "_source.eda"
    #biana_output_converter.filter_network_by_interaction_type(network_attribute_file_name = network_file_method_attribute, network_out_file_name = network_file_prefix + "_y2h.sif", interaction_type="y2h")
    #biana_output_converter.filter_network_by_interaction_type(network_attribute_file_name = network_file_method_attribute, network_out_file_name = network_file_prefix + "_tap.sif", interaction_type="tap")
    #biana_output_converter.filter_network_by_interaction_type(network_attribute_file_name = network_file_method_attribute, network_out_file_name = network_file_prefix + "_no_tap.sif", interaction_type="tap", reverse_selection=True)
    #biana_output_converter.filter_network_by_interaction_type(network_attribute_file_name = network_file_method_attribute, network_out_file_name = filtered_network_file + ".no_tap", interaction_type="tap", reverse_selection=True)
    valid_ids = set([0,4,96,676,729,19,6,7,858,59,109]) # TAP
    biana_output_converter.filter_network_by_interaction_attribute_value(network_attribute_file_name = network_file_method_attribute, network_out_file_name = filtered_network_file + ".no_tap", accept_attribute_value = lambda x: int(x) not in valid_ids)

    #interaction_to_sources = get_interaction_sources(network_file_source_attribute)
    with open(filtered_network_file, 'w') as f:
        for line in open(filtered_network_file + ".no_tap"):
            id1, dummy, id2 = line.split()
            # Filter self interactions
            if id1 == id2:
                continue
            # Remove singleton interacions (that has evidence only from one database)
            #id_pair = sorted([id1, id2])
            #if is_singleton(interaction_to_sources[(id_pair[0], id_pair[1])]):
            #   continue
            # Do not include ambigous user entities
            if id1 in ueids and id2 in ueids:
                f.write(line)
    return

def is_singleton(sources):
    if len(sources) > 1:
        return False
    # Multiple evidence from the same db counts
    source = sources[0]
    idx = source.find("(")
    if int(source[idx+1:-1]) > 1:
        return False
    return True
def get_interaction_sources(network_file_source_attribute):
    with open(network_file_source_attribute) as f:
        interaction_to_sources = {}
        for line in f:
            words = line.strip().split()
            if len(words) == 1:
                continue
            if words[3] != "=":
                raise Exception("sif attribute file format error! %s" % line)
            id_pair = sorted([words[0], words[2]])
            interaction_to_sources.setdefault((id_pair[0],id_pair[1]), []).append(words[4])
    return interaction_to_sources


def fileExist(file):
    """
    Checks if a file exists AND is a file
    """
    return os.path.exists(file) and os.path.isfile(file)

def create_directory(directory):
    """
    Checks if a directory exists and if not, creates it
    """
    try:
        os.stat(directory)
    except:
        os.mkdir(directory)
    return


if  __name__ == "__main__":
    main()
