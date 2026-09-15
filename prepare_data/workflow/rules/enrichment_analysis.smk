from os import path
import re

#####################
#   Preprocessing   #
#####################

# Gene extraction

COEXPRESSION_NETWORKS_PATTERN = "|".join(config["coexpression_networks"])

rule get_genes_from_network:
    input:
        path.join(config["network_dir"], "{network}.txt")
    output:
        "{0}/all_genes/{{network}}/MSU/all-genes.txt".format(config["raw_enrich_dir"])
    wildcard_constraints:
        network = COEXPRESSION_NETWORKS_PATTERN
    shell:
        "python scripts/network_util/get-nodes-from-network.py " \
        "{{input}} {0}/all_genes/{{wildcards.network}}/MSU " \
        "--name all-genes".format(config["raw_enrich_dir"])

# ID conversion files

rule ricegeneid_msu_to_transcript_id:
    input:
        "{0}/all_genes/{{network}}/MSU/all-genes.txt".format(config["raw_enrich_dir"])
    output:
        "{0}/temp/{{network}}/all-transcript-id.txt".format(config["raw_enrich_dir"]),
        "{0}/temp/{{network}}/all-na-transcript-id.txt".format(config["raw_enrich_dir"])
    shell:
        "Rscript --vanilla scripts/enrichment_analysis/util/ricegeneid-msu-to-transcript-id.r " \
        "-g {{input}} " \
        "-o {0}/temp/{{wildcards.network}}".format(config["raw_enrich_dir"])

rule msu_to_transcript_id:
    input:
        all_transcript="{0}/temp/{{network}}/all-transcript-id.txt".format(config["raw_enrich_dir"]),
        all_na_transcript="{0}/temp/{{network}}/all-na-transcript-id.txt".format(config["raw_enrich_dir"]),
        rap_to_msu="{0}/rap_db/{1}".format(config["raw_enrich_dir"], config["rap_to_msu_file"]),
        rap_to_transcript="{0}/rap_db/{1}".format(config["raw_enrich_dir"], config["rap_to_transcript_file"])
    output:
        "{0}/mapping/{{network}}/msu-to-transcript-id.pickle".format(config["raw_enrich_dir"])
    shell:
        "python scripts/enrichment_analysis/util/msu-to-transcript-id.py " \
        "{{input.all_transcript}} {{input.all_na_transcript}} {{input.rap_to_msu}} {{input.rap_to_transcript}} " \
        "{0}/mapping/{{wildcards.network}}".format(config["raw_enrich_dir"])

rule transcript_to_msu_id:
    input:
        "{0}/mapping/{{network}}/msu-to-transcript-id.pickle".format(config["raw_enrich_dir"])
    output:
        "{0}/msu_mapping/{{network}}/transcript-to-msu-id.pickle".format(config["gene_id_mapping_dir"])
    shell:
        "python scripts/enrichment_analysis/util/transcript-to-msu-id.py " \
        "{{input}} {0}/msu_mapping/{{wildcards.network}}".format(config["gene_id_mapping_dir"])

rule convert_all_genes_msu_to_transcript:
    input:
        all_genes = "{0}/all_genes/{{network}}/MSU/all-genes.txt".format(config["raw_enrich_dir"]),
        mapping_file="{0}/mapping/{{network}}/msu-to-transcript-id.pickle".format(config["raw_enrich_dir"])
    output:
        "{0}/all_genes/{{network}}/transcript/all-genes.tsv".format(config["raw_enrich_dir"])
    shell:
        "python scripts/enrichment_analysis/util/file-convert-msu.py " \
        "{{input.all_genes}} {{input.mapping_file}} " \
        "{0}/all_genes/{{wildcards.network}} transcript".format(config["raw_enrich_dir"])

# rule convert_all_genes_msu_to_rap:
#     input:
#         all_genes = "{0}/all_genes/{{network}}/MSU/all-genes.txt".format(config["raw_enrich_dir"]),
#         mapping_file="{0}/msu_mapping/msu_to_rap.pickle".format(config["gene_id_mapping_dir"])
#     output:
#         "{0}/all_genes/{{network}}/rap/all-genes.tsv".format(config["raw_enrich_dir"])
#     shell:
#         "python scripts/enrichment_analysis/util/file-convert-msu.py " \
#         "{{input.all_genes}} {{input.mapping_file}} " \
#         "{0}/all_genes/{{wildcards.network}} rap".format(config["raw_enrich_dir"])

# Module Conversion

# rule convert_modules_msu_to_rap:
#     input:
#         module_file = "{0}/{{network}}/MSU/{{algo}}/{{value}}/{{algo}}-module-list.tsv".format(config["network_mod_dir"]),
#         mapping_file="{0}/msu_mapping/msu_to_rap.pickle".format(config["gene_id_mapping_dir"])
#     output:
#         "{0}/{{network}}/{{algo}}/{{value}}/rap/{{algo}}-module-list.tsv".format(config["network_mod_dir"])
#     shell:
#         "python scripts/enrichment_analysis/util/file-convert-msu.py " \
#         "{{input.module_file}} {{input.mapping_file}} " \
#         "{0}/{{wildcards.network}}/{{wildcards.algo}}/{{wildcards.value}} rap ".format(config["network_mod_dir"])

rule convert_modules_msu_to_transcript:
    input:
        module_file = path.join(config["network_mod_dir"], "{network}/MSU/{algo}/{value}/{algo}-module-list.tsv"),
        mapping_file = path.join(config["raw_enrich_dir"], "mapping/{network}/msu-to-transcript-id.pickle")
    output:
        path.join(config["network_mod_dir"], "{network}/transcript/{algo}/{value}/{algo}-module-list.tsv")
    shell:
        "python scripts/enrichment_analysis/util/module-convert-msu.py " \
        "{input.module_file} {input.mapping_file} {output}"

# Ontology Preparation

rule prepare_go_annotations:
    input:
        agrigo_file= path.join((config["raw_enrich_dir"]), "go/agrigo.tsv"),
        oryzabase_file= path.join((config["raw_enrich_dir"]), "go/OryzabaseGeneListAll_20230322010000.txt"),
        rap_db_file= path.join((config["raw_enrich_dir"]), "rap_db/IRGSP-1.0_representative_annotation_2023-03-15.tsv"),
        all_genes_file= path.join((config["raw_enrich_dir"]), "all_genes/{network}/transcript/all-genes.tsv"),
        msu_to_transcript= path.join((config["raw_enrich_dir"]), "mapping/{network}/msu-to-transcript-id.pickle"),
    output:
         path.join((config["raw_enrich_dir"]), "go/{network}/go-annotations.tsv")
    shell:
        "python scripts/enrichment_analysis/util/aggregate-go-annotations.py " \
        "{input.agrigo_file} {input.oryzabase_file} {input.rap_db_file} " \
        "{input.all_genes_file} {input.msu_to_transcript} " \
         f"{path.join((config['raw_enrich_dir']), 'go/{wildcards.network}')}"

rule prepare_to_annotations:
    input:
        oryzabase_file = path.join(config["raw_enrich_dir"], "go/OryzabaseGeneListAll_20230322010000.txt")
    output:
        path.join(config["raw_enrich_dir"], "to/{network}/to-annotations.tsv"),
        path.join(config["raw_enrich_dir"], "to/{network}/to-id-to-name.tsv")
    shell:
        "python scripts/enrichment_analysis/util/aggregate-to-annotations.py " \
        "{input.oryzabase_file} " \
        f"{path.join(config['raw_enrich_dir'], 'to/{wildcards.network}')}"

rule prepare_po_annotations:
    input:
        oryzabase_file = path.join(config["raw_enrich_dir"], "go/OryzabaseGeneListAll_20230322010000.txt")
    output:
        path.join(config["raw_enrich_dir"], "po/{network}/po-annotations.tsv"),
        path.join(config["raw_enrich_dir"], "po/{network}/po-id-to-name.tsv")
    shell:
        "python scripts/enrichment_analysis/util/aggregate-po-annotations.py " \
        "{input.oryzabase_file} " \
        f"{path.join(config['raw_enrich_dir'], 'po/{wildcards.network}')}"

##################################
#   Enrichment Analysis Proper   #
##################################

checkpoint count_modules:
    # this is a checkpoint rather than a rule so that
    # get_enriched_module_list can request for the file.
    input:
        mod_list="{0}/{{network}}/MSU/{{algo}}/{{param}}/{{algo}}-module-list.tsv".format(config["network_mod_dir"])
    output:
        count_file="{0}/temp/{{network}}/{{algo}}/{{param}}/module_count.txt".format(config["raw_enrich_dir"])
    run:
        with open(input.mod_list) as f:
            line_count = sum(1 for line in f)
        with open(output.count_file, "w") as f:
            f.write(str(line_count))

def get_enriched_module_list(wildcards):
    count_file = checkpoints.count_modules.get(
        network=wildcards.network,
        algo=wildcards.algo,
        param=wildcards.param
    ).output[0]
    
    with count_file.open() as cf:
        count = int(cf.read().strip())
        
    return expand(
        path.join(config['app_enrich_dir'], 
        "{network}/output/{algo}/{param}/{enrichment_type}/{enrichment_subtype}/results/{enrichment_subtype}-df-{index}.tsv"),
        index=range(1, count + 1),
        **wildcards
    )
    
rule finish_enrichment:
    input:
        get_enriched_module_list
    output:
        temp(touch(path.join(config['app_enrich_dir'],
        "{network}/output/{algo}/{param}/{enrichment_type}/{enrichment_subtype}/results/done")))

def ceo(file_path):
    # clean enrichment output
    # unnecessary, only for legibility and style
    return re.sub(r"\/results\/\w+-df-\d+.tsv$","",file_path)

rule execute_gene_ontology_enrichment_analysis:
    input:
        mod_list = path.join(config['network_mod_dir'], '{network}/MSU/{algo}/{param}/{algo}-module-list.tsv'),
        all_genes = path.join(config['raw_enrich_dir'], 'all_genes/{network}/MSU/all-genes.txt'),
        go_annotations = path.join(config['raw_enrich_dir'], 'go/{network}/go-annotations.tsv')
    output:
        path.join(config['app_enrich_dir'], "{network}/output/{algo}/{param}/ontology_enrichment/go/results/go-df-{index}.tsv")
    params:
        output_dir = lambda wildcards, output: ceo(output[0])
    shell:
        "Rscript --vanilla scripts/enrichment_analysis/ontology_enrichment/go-enrichment.r " \
        "-g {input.mod_list} -i {wildcards.index} -b {input.all_genes} -m {input.go_annotations} " \
        "-o {params.output_dir}"
        
rule execute_plant_ontology_enrichment_analysis:
    input:
        mod_list = path.join(config['network_mod_dir'], '{network}/MSU/{algo}/{param}/{algo}-module-list.tsv'),
        all_genes = path.join(config['raw_enrich_dir'], 'all_genes/{network}/MSU/all-genes.txt'),
        po_annotations = path.join(config['raw_enrich_dir'], 'po/{network}/po-annotations.tsv'),
        po_id_to_name = path.join(config['raw_enrich_dir'], 'po/{network}/po-id-to-name.tsv'),
    output:
        path.join(config['app_enrich_dir'], "{network}/output/{algo}/{param}/ontology_enrichment/po/results/po-df-{index}.tsv")
    params:
        output_dir = lambda wildcards, output: ceo(output[0])
    shell:
        "Rscript --vanilla scripts/enrichment_analysis/ontology_enrichment/po-enrichment.r " \
        "-g {input.mod_list} -i {wildcards.index} -b {input.all_genes} -m {input.po_annotations} -t {input.po_id_to_name} " \
        "-o {params.output_dir}"

rule execute_trait_ontology_enrichment_analysis:
    input:    
        mod_list = path.join(config['network_mod_dir'], '{network}/MSU/{algo}/{param}/{algo}-module-list.tsv'),
        all_genes = path.join(config['raw_enrich_dir'], 'all_genes/{network}/MSU/all-genes.txt'),
        to_annotations = path.join(config['raw_enrich_dir'], 'to/{network}/to-annotations.tsv'),
        to_id_to_name = path.join(config['raw_enrich_dir'], 'to/{network}/to-id-to-name.tsv'),
    output:
        path.join(config['app_enrich_dir'], "{network}/output/{algo}/{param}/ontology_enrichment/to/results/to-df-{index}.tsv")
    params:
        output_dir = lambda wildcards, output: ceo(output[0])
    shell:
        "Rscript --vanilla scripts/enrichment_analysis/ontology_enrichment/to-enrichment.r " \
        "-g {input.mod_list} -i {wildcards.index} -b {input.all_genes} -m {input.to_annotations} -t {input.to_id_to_name} " \
        "-o {params.output_dir}"
 
rule execute_overrepresentation_pathway_enrichment_analysis_via_clusterprofiler:
    threads: 3 # https://www.kegg.jp/kegg/rest/ "limit to 3 requests per second"
    input:    
        mod_list = path.join(config['network_mod_dir'], '{network}/transcript/{algo}/{param}/{algo}-module-list.tsv'),
        all_genes = path.join(config['raw_enrich_dir'], 'all_genes/{network}/transcript/all-genes.tsv'),
    output:
        path.join(config['app_enrich_dir'], "{network}/output/{algo}/{param}/pathway_enrichment/ora/results/ora-df-{index}.tsv")
    params:
        output_dir = lambda wildcards, output: ceo(output[0])
    shell:
        "Rscript --vanilla scripts/enrichment_analysis/pathway_enrichment/ora-enrichment.r " \
        "-g {input.mod_list} -i {wildcards.index} -b {input.all_genes} " \
        "-o {params.output_dir}"

rule execute_topology_based_pathway_enrichment_analysis_via_pathway_express:
    threads: 3 # https://www.kegg.jp/kegg/rest/ "limit to 3 requests per second"
    input:    
        mod_list = path.join(config['network_mod_dir'], '{network}/transcript/{algo}/{param}/{algo}-module-list.tsv'),
        all_genes = path.join(config['raw_enrich_dir'], 'all_genes/{network}/transcript/all-genes.tsv'),
    output:
        path.join(config['app_enrich_dir'], "{network}/output/{algo}/{param}/pathway_enrichment/pe/results/pe-df-{index}.tsv")
    params:
        output_dir = lambda wildcards, output: ceo(output[0])
    shell:
        "Rscript --vanilla scripts/enrichment_analysis/pathway_enrichment/pe-enrichment.r " \
        "-g {input.mod_list} -i {wildcards.index} -b {input.all_genes} " \
        "-o {params.output_dir}"

rule execute_topology_based_pathway_enrichment_analysis_via_spia:
    input:    
        mod_list = path.join(config['network_mod_dir'], '{network}/transcript/{algo}/{param}/{algo}-module-list.tsv'),
        all_genes = path.join(config['raw_enrich_dir'], 'all_genes/{network}/transcript/all-genes.tsv'),
        spia_path = path.join(config['raw_enrich_dir'], 'kegg_dosa/SPIA'),
    output:
        path.join(config['app_enrich_dir'], "{network}/output/{algo}/{param}/pathway_enrichment/spia/results/spia-df-{index}.tsv")
    params:
        output_dir = lambda wildcards, output: ceo(output[0])
    shell:
        "Rscript --vanilla scripts/enrichment_analysis/pathway_enrichment/spia-enrichment.r " \
        "-g {input.mod_list} -i {wildcards.index} -b {input.all_genes} -s {input.spia_path} " \
        "-o {params.output_dir}"

