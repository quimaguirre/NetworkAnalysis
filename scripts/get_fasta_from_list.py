import argparse
import ConfigParser
import pandas as pd
import time
import os, sys, re, uuid

import biana
try: from biana import *
except: sys.exit(10)
from biana.BianaObjects import UserEntitySet


def main():
    """
    Requires:
    - BIANA package installed
    - BIANA database available

    Usage:
    python create_network.py -i <input_file> -t <type_of_protein_identifier> -taxid <taxonomy restriction (by default None = all species)) -o <output_file>

    Examples of command:
    python /home/quim/PHD/Projects/BIANA/scripts/get_fasta_from_list.py -i /home/quim/PHD/Projects/BIANA/data/human.txt -t uniprotaccession -o /home/quim/PHD/Projects/BIANA/outputs/human.fasta -v
    python /home/quim/PHD/Projects/BIANA/scripts/get_fasta_from_list.py -i /home/quim/PHD/Projects/BIANA/data/yeast.txt -t uniprotaccession -o /home/quim/PHD/Projects/BIANA/outputs/yeast.fasta -v
    python /home/quim/PHD/Projects/BIANA/scripts/get_fasta_from_list.py -i /home/quim/PHD/Projects/BIANA/data/marta_data/genes_subset.txt -t genesymbol -taxid 9606 -o /home/quim/PHD/Projects/BIANA/outputs/marta_outputs/genes_subset.fasta -v
    python /home/quim/PHD/Projects/BIANA/scripts/get_fasta_from_list.py -i /home/quim/PHD/Projects/BIANA/data/marta_data/ensembl_subset.txt -t ensembl -taxid 9606 -o /home/quim/PHD/Projects/BIANA/outputs/marta_outputs/ensembl_subset.fasta -v
    python /home/quim/PHD/Projects/BIANA/scripts/get_fasta_from_list.py -i /home/quim/PHD/Projects/BIANA/data/marta_data/genes.txt -t genesymbol -taxid 9606 -o /home/quim/PHD/Projects/BIANA/outputs/marta_outputs/genes.fasta -v
    python /home/quim/PHD/Projects/BIANA/scripts/get_fasta_from_list.py -i /home/quim/PHD/Projects/BIANA/data/marta_data/sequences.txt -t ensembl -taxid 9606 -o /home/quim/PHD/Projects/BIANA/outputs/marta_outputs/ensembl.fasta -v
    """

    options = parse_user_arguments()
    get_fasta_from_list(options)

def parse_user_arguments(*args, **kwds):

    parser = argparse.ArgumentParser(
        description = "Generate a protein-protein interaction network",
        epilog      = "@oliva's lab 2017")
    parser.add_argument('-i','--input_file',dest='input_file',action = 'store',
                        help = 'Input file')
    parser.add_argument('-t','--id_type',dest='id_type',action = 'store',default='uniprotaccession',
                        help = '''Type of identifier of the proteins 
                        e.g. uniprotaccession, uniprotentry, geneid, genesymbol 
                        (default is uniprotaccession)''')
    parser.add_argument('-taxid','--taxid',dest='taxid',action = 'store', default=None,
                        help = 'Input taxonomy identifier to restrict the results to a concrete species. default = None (all species)')
    parser.add_argument('-o','--output_file',dest='output_file',action = 'store',
                        help = 'Output file')
    parser.add_argument('-v','--verbose',dest='verbose',action = 'store_true',
                        help = 'Flag to use verbose mode')
    options=parser.parse_args()

    return options

def get_fasta_from_list(options):
    """
    Gets the protein sequences in FASTA format from a list of protein identifiers (separated by new lines).
    """

    # Start marker for time measure
    start = time.time()

    # Load config file
    scripts_path = os.path.abspath(os.path.dirname(__file__))
    config_file = os.path.join(scripts_path, 'config.ini')
    config = ConfigParser.ConfigParser()
    config.read(config_file)

    fetcher = BIANAInfoFetcher(config)
    fetcher.get_fasta(options.input_file, options.output_file, options.id_type, options.taxid, options.verbose)

    # End marker for time
    end = time.time()
    print('\n  TIME OF EXECUTION: {:.3f} seconds or {:.3f} minutes.\n'.format(end - start, (end - start) / 60))


class BIANAInfoFetcher():
    """
    Class to fetch information from BIANA
    """

    def __init__(self, config):
        """
        BIANA info fetcher
        """

        # Create BIANA session
        self.biana_session_id = str(uuid.uuid4())
        session = create_new_session( sessionID=self.biana_session_id,
                                    dbname=config.get('BIANA', 'database'),
                                    dbhost=config.get('BIANA', 'host'),
                                    dbuser=config.get('BIANA', 'user'),
                                    dbpassword=config.get('BIANA', 'password'),
                                    unification_protocol=config.get('BIANA', 'unification_protocol') )
        if session is None:
            # Mysql related error
            raise ValueError("Error in creating a new BIANA session")

        self.biana_session = available_sessions[self.biana_session_id]
        self.output_result = ''
        self.protein_to_taxid_to_sequences = {}

        return


    def get_fasta(self, input_file, output_file, id_type, taxid, verbose=False):
        """
        Reads an input file with protein identifiers and fetches their corresponding FASTA sequences.
        """

        id_type = id_type.lower()
        check_taxid = False
        if taxid and taxid.lower() != 'none':
            check_taxid = True

        # Read input file
        ids = self.read_input_file(input_file)
        #print(ids)


        # Obtain the user entity IDs associated to the input proteins
        missing_proteins = set()
        multiple_sequence_proteins = set()
        list_input_restriction_identifiers = []
        for protein_id in ids:

            # CREATE A LIST WITH ALL THE SEED IDENTIFIERS
            # Example: list_input_identifiers = [("uniprotentry","ACE_YEAST"),("uniprotentry","PGH2_HUMAN"),("uniprotentry","RND3_HUMAN")]
            #list_input_identifiers = [ (id_type, identifier) for identifier in ids ]
            list_input_identifiers = [ (id_type, protein_id) ]
            #print(list_input_identifiers)
            list_input_restriction_identifiers = [] # TAXID as a restriction attribute is too time-consuming, so better not to put it!
            list_input_negative_restriction_identifiers = []


            # CREATE THE SET
            proteome = self.biana_session.create_new_user_entity_set(   identifier_description_list =list_input_identifiers,
                                                                        attribute_restriction_list=list_input_restriction_identifiers,
                                                                        negative_attribute_restriction_list = list_input_negative_restriction_identifiers,
                                                                        id_type='embedded',
                                                                        only_uniques=True,
                                                                        new_user_entity_set_id='proteome'
                                                                    )
            if proteome:
                user_entity_ids = proteome.get_user_entity_ids()
            else:
                print('Protein {} not found in BIANA'.format(protein_id))
                missing_proteins.add(protein_id)
                continue

            #print('User entity set created.')
            #print('User entities selected: {}'.format(', '.join([str(user_entity) for user_entity in user_entity_ids])))


            # GET EXTERNAL ENTITIES ASSOCIATED WITH USER ENTITIES
            output_list = []
            node_attributes = [id_type, "taxid", "proteinsequence"]
            self.biana_session.output_user_entity_details(  user_entity_set = proteome, 
                                                            user_entity_id_list = user_entity_ids, 
                                                            out_method = self.output_method_for_biana, 
                                                            attributes = node_attributes, 
                                                            include_level_info = False, 
                                                            include_degree_info=False, 
                                                            include_tags_info=False, 
                                                            include_tags_linkage_degree_info=[], 
                                                            substitute_node_attribute_if_not_exists=False, 
                                                            output_1_value_per_attribute=True, 
                                                            output_format="tabulated", 
                                                            include_command_in_rows=False, 
                                                            output_only_unique_values=True
                                                        )

            #print('Output created.')
            #print(self.output_result)

            # Process output results
            result_attributes=['bianaid']+ node_attributes
            if(self.output_result != '' and self.output_result != None):
                n_line = 1
                for line in self.output_result.split('\n'):
                    fields = line.split('\t')
                    if n_line == 1:
                        header = line
                    elif len(fields) == len(result_attributes):
                        # Get ID
                        id_index = result_attributes.index(id_type)
                        id_result = fields[id_index]
                        # Get protein sequence
                        ps_index = result_attributes.index('proteinsequence')
                        ps_result = fields[ps_index]
                        # Check if taxID restriction
                        taxid_index = result_attributes.index('taxid')
                        taxid_result = fields[taxid_index]
                        if check_taxid:
                            if taxid_result != taxid:
                                continue
                        # Append results
                        self.protein_to_taxid_to_sequences.setdefault(id_result, {})
                        self.protein_to_taxid_to_sequences[id_result].setdefault(taxid_result, set()).add(ps_result)
                    n_line += 1

            # Clean the output result for this protein
            self.output_result = ''


        # Output FASTA
        with open(output_file, 'w') as out_fd:
            n=80 # Maximum 80 residues by line
            for protein_id in self.protein_to_taxid_to_sequences:
                for protein_taxid in self.protein_to_taxid_to_sequences[protein_id]:
                    seq = self.protein_to_taxid_to_sequences[protein_id][protein_taxid]
                    # If multiple sequences per tax ID, get the one with maximum length
                    if len(seq) > 1:
                        seq = max(seq)
                        multiple_sequence_proteins.add(protein_id)
                    else:
                        seq = list(seq)[0]
                    chunks = [seq[i:i+n] for i in range(0, len(seq), n)]
                    #print(chunks)
                    chunks = '\n'.join(chunks)
                    out_fd.write('>{} (species={})\n{}\n'.format(protein_id, protein_taxid, chunks))
        if verbose:
            print('Output file created.')

        print('Proteins with multiple sequences for the same species: {}'.format(', '.join(list(multiple_sequence_proteins))))
        print('Proteins not found in BIANA: {}'.format(', '.join(list(missing_proteins))))

        return


    def output_method_for_biana(self, output):
        """
        Append an output string to the string of output_result
        """
        self.output_result += output
        return


    def read_input_file(self, input_file):
        """
        Reads the input file.
        """
        ids = set()
        with open(input_file, 'r') as input_f:
            for line in input_f:
                fields = line.strip().split('\t')
                if fields[0] == '':
                    continue
                else:
                    ids.add(fields[0])
        return list(ids)


def fileExist(file):
    """
    Checks if a file exists AND is a file.
    """
    return os.path.exists(file) and os.path.isfile(file)


def create_directory(directory):
    """
    Checks if a directory exists and if not, creates it.
    """
    try:
        os.stat(directory)
    except:
        os.mkdir(directory)
    return


if  __name__ == "__main__":
    main()



