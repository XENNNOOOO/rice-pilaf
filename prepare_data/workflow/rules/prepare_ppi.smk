from os import path
# Resolve Network Dependency
for key, path in config['networks'].items():
    config["networks"][key] = path.format(network_dir=config["network_dir"])


# removes species id from raw data from STRING. During development,
# this is the version 12.0
rule strip_string_species_id:
    input:
        path.join(config["network_dir"], "ppi_raw/{network}.txt")
    output:
        path.join(config["network_dir"], "{network}.txt")
    script:
        "python strip_string_species_id.py "\
        "{input} {output}"

# no clue where Nb_gene_descriptions.csv is sources from
rule prepare_uniprot_to_gene:
    input:
        "{0}/Nb/Nb_gene_descriptions.csv".format(config["gene_desc_dir"])
    output:
        "{0}/msu_mapping/uniprot_to_msu.pickle".format(config["gene_id_mapping_dir"])
    shell:
        "python scripts/ppi_util/prepare_uniprot_to_gene.py " \
        "{{input}} {0}/msu_mapping uniprot_to_msu".format(config["gene_id_mapping_dir"])

rule get_proteins_from_network:
    input:
        path.join(config["network_dir"], "{network}.txt")
    output:
        "{0}/all_proteins/{{network}}/uniprot/all-proteins.txt".format(config["ppi_dir"])
    shell:
        "python scripts/network_util/get-nodes-from-network.py " \
        "{{input}} {0}/all_proteins/{{wildcards.network}}/uniprot " \
        "--name all-proteins".format(config["ppi_dir"])

rule convert_all_proteins_to_genes:
    input:
        proteins_file="{0}/all_proteins/{{network}}/uniprot/all-proteins.txt".format(config["ppi_dir"]),
        protein_to_gene_mapping="{0}/msu_mapping/uniprot_to_msu.pickle".format(config["gene_id_mapping_dir"])
    output:
        "{0}/all_genes/{{network}}/MSU/all-genes.txt".format(config["raw_enrich_dir"])
    shell:
        "python scripts/ppi_util/convert_all_prot_to_gene.py " \
        "{{input.proteins_file}} {{input.protein_to_gene_mapping}} " \
        "{0}/all_genes/{{wildcards.network}}/MSU".format(config["raw_enrich_dir"])
