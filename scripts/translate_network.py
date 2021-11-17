import argparse
import time
import sys, os, re

from context import NetworkAnalysis
import NetworkAnalysis.network_analysis as NA



def main():

    options = parse_user_arguments()
    translate_network(options)

def parse_user_arguments(*args, **kwds):
    """
    Parses the arguments of the program
    """

    parser = argparse.ArgumentParser(
        description = "Comparison of PPI networks",
        epilog      = "@oliva's lab 2017")
    parser.add_argument('-in','--input_network_file',dest='input_network_file',action = 'store',
                        help = """ Input file with a protein-protein interaction network in SIF format. """)
    parser.add_argument('-if','--input_network_format',dest='input_network_format',action = 'store',
                        help = """ Input network format. """)
    parser.add_argument('-on','--output_network_file',dest='output_network_file',action = 'store',
                        help = """ Input file with a protein-protein interaction network in SIF format. """)
    parser.add_argument('-of','--output_network_format',dest='output_network_format',action = 'store',
                        help = """ Output network format. """)
    parser.add_argument('-trans','--translation_file',dest='translation_file',action = 'store',
                        help = """ Input file with the translation file of biana codes to geneID """)

    options=parser.parse_args()

    return options

#################
#################
# MAIN FUNCTION #
#################
#################

def translate_network(options):
    """
    Translation of a network from BIANA IDs using a translation file.
    """

    # Start marker for time measure
    start = time.time()

    # Define the input network
    input_network_file = options.input_network_file
    type_id = 'biana'
    input_network_format = options.input_network_format
    network_instance = NA.Network(input_network_file, type_id, input_network_format)

    # Translate the main network
    output_network_file = options.output_network_file
    translation_file = options.translation_file
    output_network_format = options.output_network_format
    translation_id = None
    network_instance.translate_network(translation_file, translation_id, output_network_format, output_network_file)

    # End marker for time
    end = time.time()
    print('\nTIME OF EXECUTION: {:.3f} seconds or {:.3f} minutes.\n'.format(end - start, (end - start) / 60))

    return

if  __name__ == "__main__":
    main()
