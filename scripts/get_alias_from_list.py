import argparse
import ConfigParser
import sys, os, re

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
    python get_alias_from_list.py -i <input_file> -ii <input_proteins_identifier> -io <output_proteins_identifier> -o <output_file>

    Example of command:
    python /home/quim/PHD/Projects/BIANA/scripts/get_alias_from_list.py -i /home/quim/PHD/Projects/BIANA/data/yeast_symbols.txt -ii genesymbol -io sgd -o /home/quim/PHD/Projects/BIANA/outputs/yeast_sgd.fasta -v
    """

    options = parse_user_arguments()
    get_fasta(options)

def parse_user_arguments(*args, **kwds):

    parser = argparse.ArgumentParser(
        description = "Generate a protein-protein interaction network",
        epilog      = "@oliva's lab 2017")
    parser.add_argument('-i','--input_file',dest='input_file',action = 'store',
                        help = 'Input file')
    parser.add_argument('-ii','--id_input',dest='id_input',action = 'store',default='uniprotaccession',
                        help = '''Type of identifier of the input proteins 
                        e.g. uniprotaccession, uniprotentry, geneid, genesymbol 
                        (default is uniprotaccession)''')
    parser.add_argument('-io','--id_output',dest='id_output',action = 'store',default='uniprotaccession',
                        help = '''Type of identifier of the output proteins 
                        e.g. uniprotaccession, uniprotentry, geneid, genesymbol 
                        (default is uniprotaccession)''')
    parser.add_argument('-o','--output_file',dest='output_file',action = 'store',
                        help = 'Output file')
    parser.add_argument('-v','--verbose',dest='verbose',action = 'store_true',
                        help = 'Flag to use verbose mode')
    options=parser.parse_args()

    return options


def get_fasta(options):
    """
    Generates a protein-protein interaction network extracting information from BIANA.
    """
    # Load config file
    scripts_path = os.path.abspath(os.path.dirname(__file__))
    config_file = os.path.join(scripts_path, 'config.ini')
    config = ConfigParser.ConfigParser()
    config.read(config_file)


    # Read input file
    ids = read_input_file(options.input_file)
    #print(ids)


    # START BIANA SESSION
    session = create_new_session( sessionID="biana_session",
                                  dbname=config.get('BIANA', 'database'),
                                  dbhost=config.get('BIANA', 'host'),
                                  dbuser=config.get('BIANA', 'user'),
                                  dbpassword=config.get('BIANA', 'password'),
                                  unification_protocol=config.get('BIANA', 'unification_protocol') )


    # Obtain the user entity IDs associated to the input proteins
    protein_to_sequence = {}
    for protein_id in ids:
        proteome = session.create_new_user_entity_set(
                    identifier_description_list =[protein_id],
                    id_type=options.id_input,
                    new_user_entity_set_id="proteome",
                    only_uniques=True)
        if not proteome:
            print('Protein {} not found in BIANA'.format(protein_id))
            continue
        matched_user_entities = proteome.get_user_entity_ids()
        attribute_values = session.dbAccess.get_user_entity_attributes( unification_protocol_name = config.get('BIANA', 'unification_protocol'), listUserEntityID = matched_user_entities, attribute_identifier = options.id_output, only_uniques = True).values()
        if len(attribute_values) > 0: attribute_values = attribute_values[0] # Normally we obtain the list of attributes inside a list, so we can get directly the first element
        #if len(attribute_values) > 0: attribute_values = max(attribute_values) # If there are more than one sequence, we get the longest one
        protein_to_sequence[protein_id] = attribute_values
        if options.verbose:
            print('PROTEIN: {}\nSEQUENCE: {}'.format(protein_id, attribute_values))

    if len(protein_to_sequence) == 0:
        print('The proteins introduced have not been found in BIANA!')
        sys.exit(10)

    sys.exit(0)
    # Output FASTA
    with open(options.output_file, 'w') as out_fd:
        n=80 # Maximum 80 residues by line
        for protein_id in protein_to_sequence:
            seq = protein_to_sequence[protein_id]
            chunks = [seq[i:i+n] for i in range(0, len(seq), n)]
            print(chunks)
            chunks = '\n'.join(chunks)
            out_fd.write('>{}\n{}\n'.format(protein_id, chunks))
    if options.verbose:
        print('Output file created.')

    return


def read_input_file(input_file):
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



