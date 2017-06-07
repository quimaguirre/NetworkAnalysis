import cPickle
import mysql.connector
import networkx as nx
import sys, re, os

"""
    NetworkAnalysis
    2017 Joaquim Aguirre-Plans 
    Structural Bioinformatics Laboratory
    Universitat Pompeu Fabra
"""


def main():

    cnx = mysql.connector.connect(user='quim', password='',
                                  host='localhost',
                                  database='BIANA_APR_2017')

    hippie_scores_file = '/home/quim/project/tissue_specificity/hippie_scores/experimental_scores.tsv'

    BTOterm2uE, BTOname2uE = obtain_uEs_from_Tissues(cnx)
    tissue2uEs, tissue2cell2uE = obtain_uEs_from_HPA(cnx)
    UEprot2UETissues = obtain_uE_prot_2_uE_Tissues(cnx)
    UEprot2UEHPA = obtain_uE_prot_2_uE_HPA(cnx)
    psimi2method = obtain_psimi_to_method(cnx)
    wang_liver_network = obtain_wang_liver_interactions(cnx)
    psimi2score = parse_hippie_scores(hippie_scores_file)


    BTOterm_file = 'scripts/'+'BTOterm2uE.pcl'
    cPickle.dump(BTOterm2uE, open(BTOterm_file, 'w')) 

    BTOname_file = 'scripts/'+'BTOname2uE.pcl'
    cPickle.dump(BTOname2uE, open(BTOname_file, 'w')) 

    HPA_tissue_file = 'scripts/'+'tissue2uEs.pcl'
    cPickle.dump(tissue2uEs, open(HPA_tissue_file, 'w')) 

    HPA_complete_file = 'scripts/'+'tissue2cell2uE.pcl'
    cPickle.dump(tissue2cell2uE, open(HPA_complete_file, 'w')) 

    prot2tissues_file = 'scripts/'+'UEprot2UETissues.pcl'
    cPickle.dump(UEprot2UETissues, open(prot2tissues_file, 'w')) 

    prot2HPA_file = 'scripts/'+'UEprot2UEHPA.pcl'
    cPickle.dump(UEprot2UEHPA, open(prot2HPA_file, 'w')) 

    psimi2method_file = 'scripts/'+'psimi2method.pcl'
    cPickle.dump(psimi2method, open(psimi2method_file, 'w')) 

    wang_liver_file = 'scripts/'+'wang_liver_network.pcl'
    cPickle.dump(wang_liver_network, open(wang_liver_file, 'w'))

    psimi2score_file = 'scripts/'+'psimi2score.pcl'
    cPickle.dump(psimi2score, open(psimi2score_file, 'w'))

    return


def obtain_uEs_from_Tissues(cnx):
    """
    Obtain dictionary uE : {'BTO_term' : ... , 'BTO_name' : ...}
    """

    cursor = cnx.cursor()

    query = (''' SELECT uE.userEntityID, BT.value, BN.value 
                 FROM userEntityUnification_protocol_3 uE, externalEntityBTO_name BN, externalEntityBTO_term BT 
                 WHERE uE.externalEntityID = BN.externalEntityID AND uE.externalEntityID = BT.externalEntityID
             ''')

    cursor.execute(query)

    BTOterm2uE = {}
    BTOname2uE = {}

    for items in cursor:

        uE, BTO_term, BTO_name = items

        BTO_name = BTO_name.lower()

        if BTO_term not in BTOterm2uE:
            BTOterm2uE[BTO_term] = uE
        else:
            print('BTO_term {} has multiple uEs'.format(BTO_term))
            sys.exit(10)

        if BTO_name not in BTOname2uE:
            BTOname2uE[BTO_name] = uE
        else:
            print('BTO_name {} has multiple uEs'.format(BTO_name))
            sys.exit(10)

    cursor.close()

    return BTOterm2uE, BTOname2uE


def obtain_uEs_from_HPA(cnx):
    """
    Obtain dictionary uE : {'tissue_name' : ... , 'cell_type' : ...}
    """

    cursor = cnx.cursor()

    query = (''' SELECT uE.userEntityID, HPAT.value, HPAC.value 
                 FROM userEntityUnification_protocol_3 uE, externalEntityHumanProteinAtlas_tissue HPAT, externalEntityHumanProteinAtlas_celltype HPAC 
                 WHERE uE.externalEntityID = HPAT.externalEntityID AND uE.externalEntityID = HPAC.externalEntityID
             ''')

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


def obtain_uE_prot_2_uE_Tissues(cnx):
    """
    Obtain dictionary uE_prot : {'uE_tissue' : {'confidence' : ..., 'source' : ..., 'evidence' : ...} }
    """     

    print('\n.....Obtaining dictionary of user entity proteins to user entity TISSUES.....\n')

    cursor = cnx.cursor()

    query = (''' SELECT UP.userEntityID, UT.userEntityID, TC.value, TS.value, TE.value
                 FROM userEntityUnification_protocol_3 UP, externalEntityRelationParticipant RP, externalEntityRelationParticipant RT, externalEntityTissuesConfidence TC, externalEntityTissuesSource TS, externalEntityTissuesEvidence TE, userEntityUnification_protocol_3 UT, externalEntity ET 
                 WHERE UP.externalEntityID = RP.externalEntityID AND RP.externalEntityID != RT.externalEntityID AND RP.externalEntityRelationID = RT.externalEntityRelationID AND RT.externalEntityRelationID = TC.externalEntityID AND RT.externalEntityRelationID = TS.externalEntityID AND RT.externalEntityRelationID = TE.externalEntityID AND RT.externalEntityID = UT.externalEntityID AND RT.externalEntityID = ET.externalEntityID AND ET.type = 'tissue'
             ''')

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


def obtain_uE_prot_2_uE_HPA(cnx):
    """
    Obtain dictionary uE_prot : {'uE_tissue' : {'level' : ..., 'reliability' : ...} }
    """     

    print('\n.....Obtaining dictionary of user entity proteins to user entity HUMAN PROTEIN ATLAS.....\n')

    cursor = cnx.cursor()

    query = (''' SELECT UP.userEntityID, UT.userEntityID, TL.value, TR.value 
                 FROM userEntityUnification_protocol_3 UP, externalEntityRelationParticipant RP, externalEntityRelationParticipant RT, externalEntityHumanProteinAtlas_level TL, externalEntityHumanProteinAtlas_reliability TR, userEntityUnification_protocol_3 UT, externalEntity ET 
                 WHERE UP.externalEntityID = RP.externalEntityID AND RP.externalEntityID != RT.externalEntityID AND RP.externalEntityRelationID = RT.externalEntityRelationID AND RT.externalEntityRelationID = TL.externalEntityID AND RT.externalEntityRelationID = TR.externalEntityID AND RT.externalEntityID = UT.externalEntityID AND RT.externalEntityID = ET.externalEntityID AND ET.type = 'tissue'
             ''')

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

    print('\nProtein 2 HUMAN PROTEIN ATLAS dictionary obtained!\n')

    return UEprot2UEHPA


def obtain_psimi_to_method(cnx):
    """
    Obtain dictionary uE_prot : {'psi-mi code' : 'method_name' }
    """     

    print('\n.....Obtaining dictionary of PSI-MI codes to Method names.....\n')

    cursor = cnx.cursor()

    query = (''' SELECT K.value, P.value FROM externalEntitypsimi_name P, key_attribute_2 K where P.externalEntityID = K.externalEntityID
             ''')

    cursor.execute(query)

    psimi2method = {}

    for items in cursor:

        psimi = items[0]
        method = items[1]
        psimi2method[psimi] = method

    cursor.close()

    print('\nPSI-MI 2 METHOD dictionary obtained!\n')

    return psimi2method

def obtain_wang_liver_interactions(cnx):
    """
    Obtain the liver-specific interactions in Wang et al.
    """     

    print('\n.....Obtaining liver interactions from Wang.....\n')

    cursor = cnx.cursor()

    query = (''' SELECT G1.value, G2.value 
                 FROM userEntityUnification_protocol_3 U1, userEntityUnification_protocol_3 U2, externalEntityRelationParticipant R1, externalEntityRelationParticipant R2, externalEntityPubmed P, externalEntityGeneID G1, externalEntityGeneID G2 
                 WHERE U1.userEntityID != U2.userEntityID AND U1.externalEntityID = R1.externalEntityID AND U2.externalEntityID = R2.externalEntityID AND R1.externalEntityID != R2.externalEntityID AND R1.externalEntityRelationID = R2.externalEntityRelationID AND R1.externalEntityRelationID = P.externalEntityID AND P.value = 21988832 AND U1.externalEntityID = G1.externalEntityID AND U2.externalEntityID = G2.externalEntityID
             ''')

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



if  __name__ == "__main__":
    main()