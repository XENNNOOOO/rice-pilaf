# Resolve Network Dependency
for key, path in config['networks'].items():
    config["networks"][key] = path.format(network_dir=config["network_dir"])
from os import path

rule module_detect_clusterone:
    input:
        expand(
            "{mod_detect_dir}/{network}/temp/clusterone/clusterone-results-{density}.csv",
            density = config['clusterone_min_density'].keys(),
            mod_detect_dir = config['mod_detect_dir'],
            network = config["networks"].keys()
            ),
        expand(
            "{0}/{{network}}/clusterone/{{density}}/{1}/clusterone-module-list.tsv".format(config['network_mod_dir'], "uniprot"),
            density = config['clusterone_min_density'].keys(),
            network = config['networks'].keys()
        )

rule execute_clusterone:
    input:
        lambda wildcards: config["networks"][wildcards.network]
    output:
        "{mod_detect_dir}/{network}/temp/clusterone/clusterone-results-{density}.csv"
    params:
        value = lambda wildcards: config['clusterone_min_density'][int(wildcards.density)],
        clusterone_jar_path = config['clusterone_jar_path']
    shell:
        "java -jar {params.clusterone_jar_path} " \
        "--output-format csv " \
        "--min-density {params.value} " \
        "{input} > {output}"

# Note: Each different format should have its own rule
rule get_mod_clusterone_uniprot:
    input:
        "{0}/{{network}}/temp/clusterone/clusterone-results-{{density}}.csv".format(config['mod_detect_dir'])
    output:
        "{0}/{{network}}/clusterone/{{density}}/{1}/clusterone-module-list.tsv".format(config['network_mod_dir'], "uniprot")
    shell:
        "python scripts/module_util/get-modules-from-clusterone-results.py " \
        "{{input}} {0}/{{wildcards.network}}/clusterone/{{wildcards.density}}/uniprot".format(config['network_mod_dir'])


rule module_detect_fox:
    input:
        expand(
            "{mod_detect_dir}/{network}/temp/fox/fox-int-module-list-{wcc}.txt",
            mod_detect_dir = config['mod_detect_dir'],
            network = config['networks'].keys(),
            wcc = config['wcc_threshold'].keys()
        ),
        expand(
            "{0}/{{network}}/fox/{{wcc}}/{1}/fox-module-list.tsv".format(config['network_mod_dir'], "uniprot"),
            network = config['networks'].keys(),
            wcc = config['wcc_threshold'].keys()
        )

rule generate_int_edge_list:
    input:
        lambda wildcards: config["networks"][wildcards.network]
    output:
        edge_list = "{mod_detect_dir}/{network}/mapping/int-edge-list.txt",
        node_mapping = "{mod_detect_dir}/{network}/mapping/int-edge-list-node-mapping.pickle"
    shell:
        "python scripts/network_util/convert-to-int-edge-list.py " \
        "{input} {wildcards.mod_detect_dir}/{wildcards.network}/mapping"

rule execute_fox:
    input:
        "{mod_detect_dir}/{network}/mapping/int-edge-list.txt"
    output:
        "{mod_detect_dir}/{network}/temp/fox/fox-int-module-list-{wcc}.txt"
    params:
        wcc_val = lambda wildcards: config['wcc_threshold'][int(wildcards.wcc)]
    shell:
        "scripts/module_detection/execute_lazyfox.sh {input} {params.wcc_val} {output}"

rule get_mod_fox_uniprot:
    input:
        fox_result = "{0}/{{network}}/temp/fox/fox-int-module-list-{{wcc}}.txt".format(config['mod_detect_dir']),
        node_mapping = "{0}/{{network}}/mapping/int-edge-list-node-mapping.pickle".format(config['mod_detect_dir'])
    output:
        "{0}/{{network}}/fox/{{wcc}}/{1}/fox-module-list.tsv".format(config['network_mod_dir'], "uniprot")
    shell:
        "python scripts/module_util/restore-node-labels-in-modules.py " \
        "{{input.fox_result}} {{input.node_mapping}} " \
        "{0}/{{wildcards.network}}/fox/{{wildcards.wcc}}/uniprot fox".format(config['network_mod_dir'])

rule execute_coach:
    input:
        path.join(config['mod_detect_dir'], "{network}/mapping/int-edge-list.txt")
    params:
        affinity_threshold_value = lambda wildcards: config['affinity_thresholds'][wildcards.affinity_threshold]
    output:
        path.join(config['mod_detect_dir'], "{network}/temp/coach/coach-int-module-list-{affinity_threshold}.csv")
    shell:
        "python scripts/module_detection/detect-modules-via-coach.py " \
        "--affinity_threshold {params.affinity_threshold_value} {input} " \
        f"{path.join(config['mod_detect_dir'], '{wildcards.network}/temp/coach/')}"

# rule relabel_nodes_in_coach_modules:
#     input:
#         coach_results = path.join(config['mod_detect_dir'], "{network}/temp/coach/coach-int-module-list-{affinity_threshold}.csv"),
#         node_mapping = path.join(config['mod_detect_dir'], "{network}/mapping/int-edge-list-node-mapping.pickle")
#     output:
#         path.join(config['mod_detect_dir'], "{network}/coach/{affinity_threshold}/uniprot/coach-module-list.tsv")
#     shell:
#         "python scripts/module_util/restore-node-labels-in-modules.py " \
#         "{input.coach_results} {input.node_mapping} " \
#         f"{path.join(config['mod_detect_dir'], '{wildcards.network}/coach/{wildcards.affinity_threshold}/uniprot')} " \
#         "coach"

rule execute_demon:
    input:
        path.join(config['mod_detect_dir'], "{network}/mapping/int-edge-list.txt")
    params:
        merging_threshold_value = lambda wildcards: config['merging_thresholds'][wildcards.merging_threshold]
    output:
        path.join(config['mod_detect_dir'], "{network}/temp/demon/demon-int-module-list-{merging_threshold}.csv")
    shell:
        "python scripts/module_detection/detect-modules-via-demon.py " \
        "--epsilon {params.merging_threshold_value} {input} " \
        f"{path.join(config['mod_detect_dir'], '{wildcards.network}/temp/demon/')}"

# rule relabel_nodes_in_demon_modules:
#     input:
#         demon_results = path.join(config['mod_detect_dir'], "{network}/temp/demon/demon-int-module-list-{merging_threshold}.csv"),
#         node_mapping = path.join(config['mod_detect_dir'], "{network}/mapping/int-edge-list-node-mapping.pickle")
#     output:
#         path.join(config['mod_detect_dir'], "{network}/demon/{merging_threshold}/uniprot/demon-module-list.tsv")
#     shell:
#         "python scripts/module_util/restore-node-labels-in-modules.py " \
#         "{input.demon_results} {input.node_mapping} " \
#         f"{path.join(config['mod_detect_dir'], '{wildcards.network}/demon/{wildcards.merging_threshold}/uniprot')} " \
#         "demon"

rule relabel_nodes_from_int_to_id_in_detected_modules:
    input:
        demon_results = path.join(config['mod_detect_dir'], "{network}/temp/{algo}/{algo}-int-module-list-{merging_threshold}.csv"),
        node_mapping = path.join(config['mod_detect_dir'], "{network}/mapping/int-edge-list-node-mapping.pickle")
    output:
        path.join(config['mod_detect_dir'], "{network}/{algo}/{merging_threshold}/uniprot/{algo}-module-list.tsv")
    shell:
        "python scripts/module_util/restore-node-labels-in-modules.py " \
        "{input.demon_results} {input.node_mapping} " \
        f"{path.join(config['mod_detect_dir'], '{wildcards.network}/{wildcards.algo}/{wildcards.merging_threshold}/uniprot')} " \
        "{wildcards.algo}"
