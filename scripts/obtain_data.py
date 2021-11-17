import argparse
import ConfigParser
import cPickle
import mysql.connector
import networkx as nx
import sys, os, re

from context import NetworkAnalysis
import NetworkAnalysis.drug as NA_drug

"""
    NetworkAnalysis
    2017 Joaquim Aguirre-Plans 
    Structural Bioinformatics Laboratory
    Universitat Pompeu Fabra
"""


def main():

    options = parse_user_arguments()
    create_tissue_specific_network(options)

def parse_user_arguments(*args, **kwds):
    """
    Parses the arguments of the program
    """

    parser = argparse.ArgumentParser(
        description = "Obtain tissue-specificity data from BIANA and save it in Pickle files",
        epilog      = "@oliva's lab 2017")
    parser.add_argument('-p','--pickles_path',dest='pickles_path',action = 'store',default=os.path.join(os.path.join(os.path.dirname(__file__), '..'), 'NetworkAnalysis/pickles'),
                        help = """Define the directory where the data will be stored. """)

    options=parser.parse_args()

    return options


#################
#################
# MAIN FUNCTION #
#################
#################

def create_tissue_specific_network(options):
    """
    Generates the profiles of the input drug
    """

    #--------------------------------------#
    #   GET INFORMATION FROM CONFIG FILE   #
    #--------------------------------------#

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


    #----------------------#
    #   OBTAIN BTO FILES   #
    #----------------------#

    BTOterm_file = os.path.join(options.pickles_path, 'BTOterm2uE.pcl')
    BTOname_file = os.path.join(options.pickles_path, 'BTOname2uE.pcl')
    if not fileExist(BTOterm_file) or not fileExist(BTOname_file):
        BTOterm2uE, BTOname2uE = obtain_uEs_from_Tissues(cnx, config.get('BIANA', 'unification_protocol'))
        cPickle.dump(BTOterm2uE, open(BTOterm_file, 'w')) 
        cPickle.dump(BTOname2uE, open(BTOname_file, 'w')) 

    #----------------------------------#
    #   OBTAIN HPA AND TISSUES FILES   #
    #----------------------------------#

    HPA_tissue_file = os.path.join(options.pickles_path, 'tissue2uEs.pcl')
    HPA_complete_file = os.path.join(options.pickles_path, 'tissue2cell2uE.pcl')
    if not fileExist(HPA_tissue_file) or not fileExist(HPA_complete_file):
        tissue2uEs, tissue2cell2uE = obtain_uEs_from_HPA(cnx, config.get('BIANA', 'unification_protocol'))
        cPickle.dump(tissue2uEs, open(HPA_tissue_file, 'w'))
        cPickle.dump(tissue2cell2uE, open(HPA_complete_file, 'w'))

    prot2tissues_file = os.path.join(options.pickles_path, 'UEprot2UETissues.pcl')
    if not fileExist(prot2tissues_file):
        UEprot2UETissues = obtain_uE_prot_2_uE_Tissues(cnx, config.get('BIANA', 'unification_protocol'))
        cPickle.dump(UEprot2UETissues, open(prot2tissues_file, 'w')) 

    prot2HPAmic_file = os.path.join(options.pickles_path, 'UEprot2UEHPAmic.pcl')
    if not fileExist(prot2HPAmic_file):
        UEprot2UEHPAmic = obtain_uE_prot_2_uE_HPAmic(cnx, config.get('BIANA', 'unification_protocol'))
        cPickle.dump(UEprot2UEHPAmic, open(prot2HPAmic_file, 'w')) 

    prot2HPARNAseq_file = os.path.join(options.pickles_path, 'UEprot2UEHPARNAseq.pcl')
    if not fileExist(prot2HPARNAseq_file):
        UEprot2UEHPARNAseq = obtain_uE_prot_2_uE_HPARNAseq(cnx, config.get('BIANA', 'unification_protocol'))
        cPickle.dump(UEprot2UEHPARNAseq, open(prot2HPARNAseq_file, 'w')) 

    #-----------------------------#
    #   OBTAIN CODES OF METHODS   #
    #-----------------------------#

    psimi2method_file = os.path.join(options.pickles_path, 'psimi2method.pcl')
    if not fileExist(psimi2method_file):
        key_attribute_table = NA_drug.return_key_attribute_table(cnx, ontology_name='psimiobo')
        psimi2method = obtain_psimi_to_method(cnx, key_attribute_table)
        cPickle.dump(psimi2method, open(psimi2method_file, 'w')) 

    #-----------------------------------#
    #   PARSE WANG LIVER INTERACTIONS   #
    #-----------------------------------#

    wang_liver_file = os.path.join(options.pickles_path, 'wang_liver_network.pcl')
    if not fileExist(wang_liver_file):
        wang_liver_network = obtain_wang_liver_interactions(cnx, config.get('BIANA', 'unification_protocol'))
        cPickle.dump(wang_liver_network, open(wang_liver_file, 'w'))

    #-------------------------#
    #   PARSE HIPPIE SCORES   #
    #-------------------------#

    hippie_scores_file = os.path.join(main_path, 'NetworkAnalysis/hippie_scores/experimental_scores.tsv')
    psimi2score_file = os.path.join(options.pickles_path, 'psimi2score.pcl')
    if not fileExist(psimi2score_file):
        psimi2score = parse_hippie_scores(hippie_scores_file)
        cPickle.dump(psimi2score, open(psimi2score_file, 'w'))

    #------------------------------#
    #   OBTAIN HOUSEKEEPING DATA   #
    #------------------------------#

    # Obtain housekeeping genes data
    HPA_housekeeping_file = os.path.join(main_path, 'NetworkAnalysis/housekeeping/tissue_specificity_rna_any_expressed.tab')
    elieis_housekeeping_file = os.path.join(main_path, '/home/quim/project/tissue_specificity/housekeeping/HK_genes.txt')
    hpa_geneid_dump = os.path.join(options.pickles_path, 'hpa_hk_geneIDs.pcl')
    hpa_uE_dump = os.path.join(options.pickles_path, 'hpa_hk_uEs.pcl')
    elieis_geneid_dump = os.path.join(options.pickles_path, 'eisenberg_hk_geneIDs.pcl')
    elieis_uE_dump = os.path.join(options.pickles_path, 'eisenberg_hk_uEs.pcl')
    try:
        hpa_hk_geneIDs = cPickle.load(open(hpa_geneid_dump))
        hpa_hk_uEs = cPickle.load(open(hpa_uE_dump))
        eisenberg_hk_geneIDs = cPickle.load(open(elieis_geneid_dump))
        eisenberg_hk_uEs = cPickle.load(open(elieis_uE_dump))
    except:
        hpa_hk_geneIDs, hpa_hk_uEs, eisenberg_hk_geneIDs, eisenberg_hk_uEs = obtain_housekeeping_genes(cnx, config.get('BIANA', 'unification_protocol'), HPA_housekeeping_file, elieis_housekeeping_file, hpa_geneid_dump, hpa_uE_dump, elieis_geneid_dump, elieis_uE_dump)


    return


#######################
#######################
# SECONDARY FUNCTIONS #
#######################
#######################


def fileExist(file):
    """
    Checks if a file exists AND is a file
    """
    return os.path.exists(file) and os.path.isfile(file)

def obtain_uEs_from_Tissues(cnx, unification_protocol):
    """
    Obtain dictionary uE : {'BTO_term' : ... , 'BTO_name' : ...}
    """

    up_table = NA_drug.return_unification_protocol_table(cnx, unification_protocol)
    cursor = cnx.cursor()

    # Get the user entity of the TISSUE, not the BTOelement!!
    # Because BTO and TISSUES are not properly unified, as they have different entities (tissue and BTOelement).
    # So, using this command, I am able to get the user entity of the tissue, and not the BTOelement.
    query = (''' SELECT U.userEntityID, B1.value, N.value 
                 FROM externalEntityBTO_name N, externalEntityBTO B1, externalEntityBTO B2, {} U 
                 WHERE N.externalEntityID = B1.externalEntityID AND B1.value = B2.value AND B1.externalEntityID != B2.externalEntityID AND B2.externalEntityID = U.externalEntityID
             '''.format(up_table))

    cursor.execute(query)

    BTOterm2uE = {}
    BTOname2uE = {}

    for items in cursor:

        uE, BTO_term, BTO_name = items

        BTO_name = BTO_name.lower()

        if BTO_term not in BTOterm2uE:
            BTOterm2uE[BTO_term] = uE
        else:
            if uE != BTOterm2uE[BTO_term]:
                print('BTO_term {} has multiple uEs'.format(BTO_term))
                sys.exit(10)

        # There can be more than one BTO term with the same BTO name
        # This is why we add the user entities in a set
        if BTO_name not in BTOname2uE:
            BTOname2uE.setdefault(BTO_name, set())
            BTOname2uE[BTO_name].add(uE)

    cursor.close()

    return BTOterm2uE, BTOname2uE


def obtain_uEs_from_HPA(cnx, unification_protocol):
    """
    Obtain dictionary uE : {'tissue_name' : ... , 'cell_type' : ...}
    """

    up_table = NA_drug.return_unification_protocol_table(cnx, unification_protocol)
    cursor = cnx.cursor()

    query = (''' SELECT uE.userEntityID, HPAT.value, HPAC.value 
                 FROM {} uE, externalEntityHumanProteinAtlas_tissue HPAT, externalEntityHumanProteinAtlas_celltype HPAC 
                 WHERE uE.externalEntityID = HPAT.externalEntityID AND uE.externalEntityID = HPAC.externalEntityID
             '''.format(up_table))

    cursor.execute(query)

    tissue2uEs = {}
    tissue2cell2uE = {}

    for items in cursor:
        uE, tissue_name, cell_type = items

        tissue_name = tissue_name.lower()
        cell_type = cell_type.lower()

        tissue2uEs.setdefault(tissue_name, set())
        tissue2uEs[tissue_name].add(uE)

        tissue2cell2uE.setdefault(tissue_name, {})
        if cell_type not in tissue2cell2uE[tissue_name]:
            tissue2cell2uE[tissue_name][cell_type] = uE
        else:
            print('Tissue {} and cell_type {} have multiple uEs'.format(tissue_name, cell_type))
            sys.exit(10)

    cursor.close()

    return tissue2uEs, tissue2cell2uE


def obtain_uE_prot_2_uE_Tissues(cnx, unification_protocol):
    """
    Obtain dictionary uE_prot : {'uE_tissue' : {'confidence' : ..., 'source' : ..., 'evidence' : ...} }
    """     

    print('\n.....Obtaining dictionary of user entity proteins to user entity TISSUES.....\n')

    up_table = NA_drug.return_unification_protocol_table(cnx, unification_protocol)
    cursor = cnx.cursor()

    query = (''' SELECT UP.userEntityID, UT.userEntityID, TC.value, TS.value, TE.value
                 FROM {} UP, externalEntityRelationParticipant RP, externalEntityRelationParticipant RT, externalEntityTissuesConfidence TC, externalEntityTissuesSource TS, externalEntityTissuesEvidence TE, {} UT, externalEntity ET 
                 WHERE UP.externalEntityID = RP.externalEntityID AND RP.externalEntityID != RT.externalEntityID AND RP.externalEntityRelationID = RT.externalEntityRelationID AND RT.externalEntityRelationID = TC.externalEntityID AND RT.externalEntityRelationID = TS.externalEntityID AND RT.externalEntityRelationID = TE.externalEntityID AND RT.externalEntityID = UT.externalEntityID AND RT.externalEntityID = ET.externalEntityID AND ET.type = 'tissue'
             '''.format(up_table, up_table))

    cursor.execute(query)

    UEprot2UETissues = {}

    for items in cursor:

        uEprot, uEtissues, confidence, source, evidence = items

        source = source.lower()

        UEprot2UETissues.setdefault(uEprot, {})
        UEprot2UETissues[uEprot].setdefault(uEtissues, {})
        UEprot2UETissues[uEprot][uEtissues]['confidence'] = confidence
        UEprot2UETissues[uEprot][uEtissues]['source'] = source
        UEprot2UETissues[uEprot][uEtissues]['evidence'] = source

    cursor.close()

    print('\nProtein 2 TISSUES dictionary obtained!\n')

    return UEprot2UETissues


def obtain_uE_prot_2_uE_HPAmic(cnx, unification_protocol):
    """
    Obtain dictionary uE_prot : {'uE_tissue' : {'level' : ..., 'reliability' : ...} }
    """     

    print('\n.....Obtaining dictionary of user entity proteins to user entity HUMAN PROTEIN ATLAS (microarray).....\n')

    up_table = NA_drug.return_unification_protocol_table(cnx, unification_protocol)
    cursor = cnx.cursor()

    query = (''' SELECT UP.userEntityID, UT.userEntityID, TL.value, TR.value 
                 FROM {} UP, externalEntityRelationParticipant RP, externalEntityRelationParticipant RT, externalEntityHumanProteinAtlas_level TL, externalEntityHumanProteinAtlas_reliability TR, {} UT, externalEntity ET 
                 WHERE UP.externalEntityID = RP.externalEntityID AND RP.externalEntityID != RT.externalEntityID AND RP.externalEntityRelationID = RT.externalEntityRelationID AND RT.externalEntityRelationID = TL.externalEntityID AND RT.externalEntityRelationID = TR.externalEntityID AND RT.externalEntityID = UT.externalEntityID AND RT.externalEntityID = ET.externalEntityID AND ET.type = 'tissue'
             '''.format(up_table, up_table))

    cursor.execute(query)

    UEprot2UEHPA = {}

    for items in cursor:

        uEprot, uEtissues, level, reliability = items

        level = level.lower()
        reliability = reliability.lower()

        UEprot2UEHPA.setdefault(uEprot, {})
        UEprot2UEHPA[uEprot].setdefault(uEtissues, {})
        UEprot2UEHPA[uEprot][uEtissues]['level'] = level
        UEprot2UEHPA[uEprot][uEtissues]['reliability'] = reliability

    cursor.close()

    print('\nProtein 2 HUMAN PROTEIN ATLAS (microarray) dictionary obtained!\n')

    return UEprot2UEHPA


def obtain_uE_prot_2_uE_HPARNAseq(cnx, unification_protocol):
    """
    Obtain dictionary uE_prot : {'uE_tissue' : {'value' : ..., 'unit' : ...} }
    """     

    print('\n.....Obtaining dictionary of user entity proteins to user entity HUMAN PROTEIN ATLAS (RNAseq).....\n')

    up_table = NA_drug.return_unification_protocol_table(cnx, unification_protocol)
    cursor = cnx.cursor()

    query = (''' SELECT UP.userEntityID, UT.userEntityID, RNA.value, RNA.unit 
                 FROM {} UP, externalEntityRelationParticipant RP, externalEntityRelationParticipant RT, externalEntityHumanProteinAtlas_RNAseq_value RNA, {} UT, externalEntity ET 
                 WHERE UP.externalEntityID = RP.externalEntityID AND RP.externalEntityID != RT.externalEntityID AND RP.externalEntityRelationID = RT.externalEntityRelationID AND RT.externalEntityRelationID = RNA.externalEntityID AND RT.externalEntityID = UT.externalEntityID AND RT.externalEntityID = ET.externalEntityID AND ET.type = 'tissue'
             '''.format(up_table, up_table))

    cursor.execute(query)

    UEprot2UEHPARNAseq = {}

    for items in cursor:

        uEprot, uEtissues, value, unit = items

        value = float(value)
        unit = unit.lower()

        if unit != 'tpm':
            print('Incorrect RNAseq unit for uE protein {} and uE tissue {}: {}'.format(uEprot, uEtissues, unit))
            sys.exit(10)

        UEprot2UEHPARNAseq.setdefault(uEprot, {})
        UEprot2UEHPARNAseq[uEprot].setdefault(uEtissues, {})
        UEprot2UEHPARNAseq[uEprot][uEtissues]['value'] = value
        UEprot2UEHPARNAseq[uEprot][uEtissues]['unit'] = unit

    cursor.close()

    print('\nProtein 2 HUMAN PROTEIN ATLAS (RNAseq) dictionary obtained!\n')

    return UEprot2UEHPARNAseq


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


def obtain_wang_liver_interactions(cnx, unification_protocol):
    """
    Obtain the liver-specific interactions in Wang et al.
    """     

    print('\n.....Obtaining liver interactions from Wang.....\n')

    up_table = NA_drug.return_unification_protocol_table(cnx, unification_protocol)
    cursor = cnx.cursor()

    query = (''' SELECT G1.value, G2.value 
                 FROM {} U1, {} U2, externalEntityRelationParticipant R1, externalEntityRelationParticipant R2, externalEntityPubmed P, externalEntityGeneID G1, externalEntityGeneID G2 
                 WHERE U1.userEntityID != U2.userEntityID AND U1.externalEntityID = R1.externalEntityID AND U2.externalEntityID = R2.externalEntityID AND R1.externalEntityID != R2.externalEntityID AND R1.externalEntityRelationID = R2.externalEntityRelationID AND R1.externalEntityRelationID = P.externalEntityID AND P.value = 21988832 AND U1.externalEntityID = G1.externalEntityID AND U2.externalEntityID = G2.externalEntityID
             '''.format(up_table, up_table))

    cursor.execute(query)

    G=nx.Graph()

    for items in cursor:

        node1 = items[0]
        node2 = items[1]
        G.add_edge(node1,node2)

    cursor.close()

    print('\nWang liver interactions obtained!\n')

    return G

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

def obtain_housekeeping_genes(cnx, unification_protocol, HPA_housekeeping_file, elieis_housekeeping_file, HPA_geneid_output, HPA_uE_output, eisenberg_geneid_output, eisenberg_uE_output):
    """
    Parses the housekeeping files and obtains a dictionary with them.
    """

    print('\n.....The housekeeping set of genes has not been parsed. Parsing the files.....\n')

    # Obtain Human Protein Atlas HOUSKEEPING geneIDs --> http://www.proteinatlas.org/humanproteome/housekeeping

    hpa_hk_geneIDs = set()
    hpa_hk_uEs = set()

    hpa_housekeeping_fd = open(HPA_housekeeping_file, 'r')

    first_line = hpa_housekeeping_fd.readline()

    fields_dict = obtain_header_fields(first_line)
    # Gene  Gene synonym    Ensembl Gene description    Chromosome  Position    Protein class   Evidence    HPA evidence    UniProt evidence    MS evidence Antibody    Reliability (IH)    Reliability (Mouse Brain)   Reliability (IF)    Subcellular location    RNA tissue category RNA TS  RNA TS TPM  TPM max in non-specific

    for line in hpa_housekeeping_fd:
        fields = line.strip().split("\t")
        gene = fields[ fields_dict['gene'] ]
        ensembl = fields[ fields_dict['ensembl'] ]
        reliability = fields[ fields_dict['reliability (ih)'] ].lower()
        if reliability != '' and reliability != '-' and reliability != 'uncertain':
            print(reliability)
            pass
        else:
            print('Skipping {} for low reliability'.format(gene))
            continue
        if ensembl != '' and ensembl != '-':
            uEs, geneids = obtain_uE_and_geneid_from_ensembl(cnx, unification_protocol, ensembl)
            for geneid in geneids:
                hpa_hk_geneIDs.add(int(geneid))
            for uE in uEs:
                hpa_hk_uEs.add(int(uE))
        else:
            print('Missing ensembl in HPA housekeeping for gene: {}'.format(gene))
            sys.exit(10)

    hpa_housekeeping_fd.close()


    # Obtain Eisenberg HOUSKEEPING geneIDs --> https://www.tau.ac.il/~elieis/HKG/
    
    eisenberg_hk_geneIDs = set()
    eisenberg_hk_uEs = set()

    elieis_housekeeping_fd = open(elieis_housekeeping_file, 'r')

    for line in elieis_housekeeping_fd:
        fields = line.strip().split("\t")
        gene = fields[0]
        refseq = fields[1]
        if refseq != '' or refseq != '-':
            uEs, geneids = obtain_uE_and_geneid_from_refseq(cnx, unification_protocol, refseq)
            for geneid in geneids:
                eisenberg_hk_geneIDs.add(int(geneid))
            for uE in uEs:
                eisenberg_hk_uEs.add(int(uE))
        else:
            print('Missing RefSeq in HPA housekeeping for gene: {}'.format(gene))
            sys.exit(10)

    elieis_housekeeping_fd.close()

    cPickle.dump(hpa_hk_geneIDs, open(HPA_geneid_output, 'w')) 
    cPickle.dump(hpa_hk_uEs, open(HPA_uE_output, 'w')) 
    cPickle.dump(eisenberg_hk_geneIDs, open(eisenberg_geneid_output, 'w')) 
    cPickle.dump(eisenberg_hk_uEs, open(eisenberg_uE_output, 'w')) 

    return hpa_hk_geneIDs, hpa_hk_uEs, eisenberg_hk_geneIDs, eisenberg_hk_uEs

def obtain_uE_and_geneid_from_ensembl(cnx, unification_protocol, ensembl):
    """
    Obtain geneIDs and user entities sets from their corresponding Ensembl.
    """     

    up_table = NA_drug.return_unification_protocol_table(cnx, unification_protocol)
    cursor = cnx.cursor()

    query = (''' SELECT U1.userEntityID, G.value 
                 FROM externalEntityEnsembl EN, {} U1, {} U2, externalEntityGeneID G 
                 WHERE EN.externalEntityID = U1.externalEntityID AND U1.userEntityID = U2.userEntityID AND U2.externalEntityID = G.externalEntityID AND EN.value = %s
             '''.format(up_table, up_table))

    cursor.execute(query, (ensembl,))

    uEs = set()
    geneids = set()

    for items in cursor:

        uE, geneid = items

        uEs.add(uE)
        geneids.add(geneid)

    cursor.close()

    return uEs, geneids

def obtain_uE_and_geneid_from_refseq(cnx, unification_protocol, refseq):
    """
    Obtain geneIDs and user entities sets from their corresponding RefSeq
    """     

    up_table = NA_drug.return_unification_protocol_table(cnx, unification_protocol)
    cursor = cnx.cursor()

    query = (''' SELECT U1.userEntityID, G.value 
                 FROM externalEntityRefSeq R, {} U1, {} U2, externalEntityGeneID G 
                 WHERE R.externalEntityID = U1.externalEntityID AND U1.userEntityID = U2.userEntityID AND U2.externalEntityID = G.externalEntityID AND R.value = %s
             '''.format(up_table, up_table))

    cursor.execute(query, (refseq,))

    uEs = set()
    geneids = set()

    for items in cursor:

        uE, geneid = items

        uEs.add(uE)
        geneids.add(geneid)

    cursor.close()

    return uEs, geneids

def obtain_header_fields(first_line, sep='\t'):
    """ 
    Obtain a dictionary: "field_name" => "position" 
    """
    fields_dict = {}

    header_fields = first_line.strip().split(sep)
    for x in xrange(0, len(header_fields)):
        fields_dict[header_fields[x].lower()] = x

    return fields_dict



if  __name__ == "__main__":
    main()