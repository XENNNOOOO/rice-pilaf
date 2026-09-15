from os import path

################################
# Install Detection Algorithms #
################################

rule install_clusterone_and_fox:
    output:
        "scripts/module_detection/cluster_one-1.0.jar",
        "scripts/module_detection/LazyFox"
    shell:
        "scripts/module_detection/install_clusterone_and_lazyfox.sh"

###########################
# Module Detection Proper #
###########################

# CLUSTERONE 

rule execute_clusterone:
    input:
        network = path.join(config['network_dir'], "{network}.txt"),
        clusterone_jar_path = "scripts/module_detection/cluster_one-1.0.jar"
    output:
        path.join(config['mod_detect_dir'], "{network}/temp/clusterone/clusterone-results-{param}.csv")
    params:
        min_density = lambda wildcards: config['module_detection_params']['clusterone'][wildcards.param],
    shell:
        "java -jar {input.clusterone_jar_path} " \
        "--output-format csv " \
        "--min-density {params.min_density} " \
        "{input.network} > {output}"

rule get_clusterone_modules:
    input:
        path.join(config['mod_detect_dir'], "{network}/temp/clusterone/clusterone-results-{param}.csv"),
    output:
        path.join(config['network_mod_dir'], "{network}/clusterone/{param}/{format}/clusterone-module-list.tsv")
    wildcard_constraints:
        format = r"MSU|uniprot"
    shell:
        "python scripts/module_util/get-modules-from-clusterone-results.py " \
        "{input} "\
        f"{path.join(config['network_mod_dir'], '{wildcards.network}/clusterone/{wildcards.param}/{wildcards.format}')}"


# UTILITY (For FOX, COACH, and DEMON)

rule generate_int_edge_list:
    input:
        path.join(config['network_dir'], "{network}.txt")
    output:
        edge_list = path.join(config['mod_detect_dir'], "{network}/mapping/int-edge-list.txt"),
        node_mapping = path.join(config['mod_detect_dir'], "{network}/mapping/int-edge-list-node-mapping.pickle")
    shell:
        "python scripts/network_util/convert-to-int-edge-list.py " \
        "{input} "\
        f"{path.join(config['mod_detect_dir'], '{wildcards.network}/mapping')}"

rule restore_node_labels_from_int_to_id:
    input:
        result = path.join(config['mod_detect_dir'], "{network}/temp/{algo}/{algo}-int-module-list-{param}.csv"),
        node_mapping = path.join(config['mod_detect_dir'], "{network}/mapping/int-edge-list-node-mapping.pickle"),
    output:
        path.join(config['network_mod_dir'], "{network}/{algo}/{param}/{format}/{algo}-module-list.tsv")
    wildcard_constraints:
        format = r"MSU|uniprot"
    shell:
        "python scripts/module_util/restore-node-labels-in-modules.py "
        "{input.result} {input.node_mapping} "
        f"{path.join(config['network_mod_dir'], '{wildcards.network}/{wildcards.algo}/{wildcards.param}/{wildcards.format}')} "
        "{wildcards.algo}"

# FOX

rule execute_fox:
    # the script creates a temp folder and removes it after. 
    # parallel calls to this script will cause errors.
    threads: 1
    input:
        path.join(config['mod_detect_dir'], "{network}/mapping/int-edge-list.txt")
    output:
        path.join(config['mod_detect_dir'], "{network}/temp/fox/fox-int-module-list-{param}.csv")
    params:
        wcc_threshold = lambda wildcards: config['module_detection_params']['fox'][wildcards.param]
    shell:
        "scripts/module_detection/execute_lazyfox.sh {input} {params.wcc_threshold} {output}"

# COACH

rule execute_coach:
    input:
        path.join(config['mod_detect_dir'], "{network}/mapping/int-edge-list.txt")
    params:
        affinity_threshold_value = lambda wildcards: config['module_detection_params']['coach'][wildcards.param]
    output:
        path.join(config['mod_detect_dir'], "{network}/temp/coach/coach-int-module-list-{param}.csv")
    shell:
        "python scripts/module_detection/detect-modules-via-coach.py " \
        "--affinity_threshold {params.affinity_threshold_value} {input} " \
        f"{path.join(config['mod_detect_dir'], '{wildcards.network}/temp/coach/')}"

# DEMON

rule execute_demon:
    input:
        path.join(config['mod_detect_dir'], "{network}/mapping/int-edge-list.txt")
    params:
        merging_threshold_value = lambda wildcards: config['module_detection_params']['demon'][wildcards.param]
    output:
        path.join(config['mod_detect_dir'], "{network}/temp/demon/demon-int-module-list-{param}.csv")
    shell:
        "python scripts/module_detection/detect-modules-via-demon.py " \
        "--epsilon {params.merging_threshold_value} {input} " \
        f"{path.join(config['mod_detect_dir'], '{wildcards.network}/temp/demon/')}"