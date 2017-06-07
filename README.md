
# NetworkAnalysis

2017 Joaquim Aguirre-Plans 
Structural Bioinformatics Laboratory
Universitat Pompeu Fabra

The NetworkAnalysis package permits to manage networks retrieved from BIANA 
(Biological Interactions and Network Analysis), obtain tissue-specific networks
and validate them.

## Getting started

### 1. Create the network using the following script (requires BIANA):

```
python scripts/generate_netscore_files_vapr2017.py -radius 0 -taxid 9606 -node main/human_nodes.eAFF.biana.010617.txt -edge main/human_edges.eAFF.biana.010617.txt -trans main/human_network.eAFF.010617 -ttype geneid -format multi-fields -eAFF -v
```

* radius	: Expansion of the network. If 0, it plots the complete interactome
* taxid		: Taxonomy ID used to filter the network by species-specific interactions
* node 		: Name of the output file that will contain the nodes of the network (in BIANA codes)
* edge 		: Name of the output file that will contain the edges of the network (in BIANA codes)
* trans 	: Name of the output file that will contain the translation from BIANA codes to <ttype>
* ttype 	: Type of code in which you want to translate your network
* format 	: Format of the network. It can be multi-fields or guild
* eAFF 		: Filter the network by the methods included in the affinity_dict dictionary inside the script
* v 		: Verbose

(There is an example of network in the package)



### 2. Load the network using the Network class of our package

```
network_file = 'human_edges.eAFF.biana.010617.txt'
node_file = 'human_nodes.eAFF.biana.010617.txt'
type_id = 'biana'
network_format = 'multi-fields'

main_network = NA.Network(network_file, node_file, type_id, network_format)
```


### 3. Do you want to play with the nodes or edges? Easy:

```
main_network.get_edges()
main_network.get_nodes()
```


### 4. Translate the network to <ttype> with this fancy method!

```
translation_file = 'human_network.eAFF.010617.geneid.trans')
translation_id = 'geneid'
translated_network = 'human_edges.eAFF.geneid.txt'
translated_nodes = 'human_nodes.eAFF.geneid.txt'
geneid_net = network.translate_network(translation_file, translation_id, network_format, translated_network, translated_nodes)
```


### 5. Create more specific networks by filtering the main network:

#### 5.1. Filter by psimi method IDs:

* Excluding the interactions only containing psimi methods in the exclusion list:

```
	method_ids_excluded = [512, 1024, 13, 1038]
	method_network_file = 'human_edges.eAFF.biana.method.txt'
	method_nodes_file = 'human_nodes.eAFF.biana.method.txt'
	method_network = main_network.filter_network_by_method(methods_excluded=None, method_ids_excluded=method_ids_excluded, methods_included=None, method_ids_included=None, output_network_file=method_network_file, output_nodes_file=method_nodes_file)
```

* Including interactions that are at least reported by psimi methods in the inclusion list:

```
	method_ids_included = [18, 1112]
	method_network_file = 'human_edges.eAFF.biana.method.txt'
	method_nodes_file = 'human_nodes.eAFF.biana.method.txt'
	method_network = main_network.filter_network_by_method(methods_excluded=None, method_ids_excluded=None, methods_included=None, method_ids_included=method_ids_included, output_network_file=method_network_file, output_nodes_file=method_nodes_file)
```

#### 5.2. Filter by number of pubmeds

Including only interactions with equal or larger number of pubmed IDs than in min_num_pubmeds:

```
pmid_network_file = 'human_edges.eAFF.biana.pmid.txt'
pmid_nodes_file = 'human_nodes.eAFF.biana.pmid.txt'
min_num_pubmeds = 2
pmid_network = main_network.filter_network_by_number_pubmeds(min_num_pubmeds, output_network_file=pmid_network_file, output_nodes_file=pmid_nodes_file)
```

#### 5.3. Filter by database:

Including only interactions from the sources included in the databases_included list:

```
databases = {
    'intact' : 'intact [release 2017_04 of 05-apr-2017]',
    'biogrid': 'biogrid [release 3.4.147 (31-mar-2017)]',
    'irefindex' : 'irefindex [14.0 (last edited in 2016-07-23)]',
    'hippie': 'hippie [v2.0 (06/24/2016)]',
    'hprd' : 'hprd [release 2010_04 of 13-apr-2010]',
    'dip': 'dip [release 2017_02 of 05-feb-2017]',
    'gpcr' : 'gpcr [1-apr-2017]'
}
intgrid_network_file = 'human_edges.eAFF.biana.intgrid.txt'
intgrid_nodes_file = 'human_nodes.eAFF.biana.intgrid.txt'
databases_included = [databases['intact'], databases['biogrid']]
intgrid_network = main_network.filter_network_by_database(databases_included, intgrid_network_file, intgrid_nodes_file)
```


### 6. Tissue specificity!

Create tissue-specific networks using information from Human Protein Atlas and/or Tissues (from Jensen Lab)

#### 6.1. First, define a tissue:

```
tissue_terms_hpa = 'kidney'
tissue_terms_jensen = 'kidney'
kidney = NA.Tissue(tissue_terms_hpa, tissue_terms_jensen, jensen_conf=3, hpa_level='low', hpa_rel='approved')
```

#### 6.2. Second, filter the network by the proteins in this tissue:

```
tissue_network_file = 'kidney_edges.biana.biana.txt'
tissue_nodes_file = 'kidney_nodes.biana.biana.txt'
permission = 0
tissue_network = main_network.filter_network_by_tissue(tissue_network_file, tissue_nodes_file, kidney, permission)
```

#### 6.3. Validate the network using housekeeping datasets of Human Protein Atlas and Einseberg dataset

* Obtain the nodes/edges from Human Protein Atlas, from Einsemberg, the intersection and the union:

```
    hpa_edges=tissue_network.hpa_edges
    jensen_edges=tissue_network.jensen_edges
    union=tissue_network.get_union()
    intersection=tissue_network.get_intersection()
```

* Calculate coverage of housekeeping genes in the tissue-specific network:

```
    num_all_hk_tissue = float(tissue_network.get_number_of_housekeeping_genes(all_hk_genes))
    num_hpa_hk_tissue = float(tissue_network.get_number_of_housekeeping_genes(hpa_hk_genes))
    num_eis_hk_tissue = float(tissue_network.get_number_of_housekeeping_genes(eisenberg_hk_genes))
    num_all_nonhk_tissue = float(len(set(tissue_network.get_nodes()))) - num_all_hk_tissue
    per_all_hk_tissue = num_all_hk_tissue / num_all_hk_network * 100
    per_hpa_hk_tissue = num_hpa_hk_tissue / num_hpa_hk_network * 100
    per_eis_hk_tissue = num_eis_hk_tissue / num_eis_hk_network * 100
```

* Calculate a contingency table:

```
    contingency_table = np.array([[num_all_hk_tissue, num_all_nonhk_tissue], [num_all_hk_network, num_all_nonhk_network]])
    chi2, pval, dof, expected = NA.calculate_contingency_table(contingency_table)
```

* Print the results:

```
    print('{} {}-specific network'.format(abbr.upper(), abbr_tis.upper()))
    print('Number of edges: {}'.format(len(tissue_network.get_edges())))
    print('Number of nodes: {}'.format(len(tissue_network.get_nodes())))
    print('Interactions using only Tissues (Jensen lab): {}'.format(len(jensen_edges)))
    print('Interactions using only Human Protein Atlas: {}'.format(len(hpa_edges)))
    print('Intersection: {}'.format(len(intersection)))
    print('Union: {}\n'.format(len(union)))
    print('Percentage of (all) HK genes in the tissue with respect to main network: {:.2f}%\t{:.0f} of {:.0f} genes'.format(per_all_hk_tissue, num_all_hk_tissue, num_all_hk_network))
    print('Percentage of (HPA) HK genes in the tissue with respect to main network: {:.2f}%\t{:.0f} of {:.0f} genes'.format(per_hpa_hk_tissue, num_hpa_hk_tissue, num_hpa_hk_network))
    print('Percentage of (Eisenberg) HK genes in the tissue with respect to main network: {:.2f}%\t{:.0f} of {:.0f} genes'.format(per_eis_hk_tissue, num_eis_hk_tissue, num_eis_hk_network))
    print('Contingency table: {}'.format(contingency_table))
    print('Contingency table result:\tchi2: {}\tp-value: {}\tdegrees of freedom: {}\n'.format(chi2, pval, dof))
```


### 7. Example:

Example of a complete analysis by running the following script:

```
python create_tissue_specific_networks_v1.py
```


