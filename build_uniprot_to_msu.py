"""
Builds static/app_data/gene_id_mapping/msu_mapping/uniprot_to_msu.pickle

Pipeline:
  1. Collect the unique UniProt accessions actually used in the STRING
     network (static/app_data/networks/STRING-Physical.txt).
  2. Submit them to UniProt's ID mapping API (UniProtKB_AC-ID -> EnsemblPlants),
     which returns RAP-style gene IDs (e.g. Os01g0100100) for rice.
  3. Translate those RAP IDs to MSU (LOC_Os...) IDs using RAP-MSU.txt
     (downloaded from RAP-DB's "ID converter" section).
  4. Pickle a dict: { uniprot_accession: [msu_gene_id, ...], ... }

Run from the repo root:
    python3 build_uniprot_to_msu.py \
        static/app_data/networks/STRING-Physical.txt \
        /tmp/RAP-MSU.txt \
        static/app_data/gene_id_mapping/msu_mapping/uniprot_to_msu.pickle
"""

import os
import pickle
import sys
import time

import requests

UNIPROT_API = "https://rest.uniprot.org/idmapping"


def collect_accessions(network_file):
    accessions = set()
    with open(network_file) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            accessions.add(parts[0])
            accessions.add(parts[1])
    print(f"Found {len(accessions)} unique UniProt accessions in {network_file}")
    return sorted(accessions)


def load_rap_to_msu(rap_msu_file):
    """
    RAP-MSU.txt format: RAP_ID<TAB>MSU_ID.transcript,MSU_ID.transcript,...  (or 'None')
    Strips transcript suffixes (.1, .2, ...) and dedupes to gene-level MSU IDs.
    """
    rap_to_msu = {}
    with open(rap_msu_file) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 2:
                continue
            rap_id, msu_field = parts
            if msu_field == "None":
                continue
            genes = set()
            for msu_transcript in msu_field.split(","):
                gene = msu_transcript.split(".")[0]
                if gene:
                    genes.add(gene)
            if genes:
                rap_to_msu[rap_id] = sorted(genes)
    print(f"Loaded {len(rap_to_msu)} RAP->MSU entries from {rap_msu_file}")
    return rap_to_msu


def submit_id_mapping_job(accessions, batch_size=5000):
    """
    UniProt's ID mapping API is async: submit a job, poll until finished,
    then page through results. Batches large accession lists to stay well
    under UniProt's per-request limits.
    """
    all_results = []  # list of (from_accession, to_rap_id)

    for i in range(0, len(accessions), batch_size):
        batch = accessions[i : i + batch_size]
        print(f"Submitting batch {i // batch_size + 1} ({len(batch)} accessions)...")

        resp = requests.post(
            f"{UNIPROT_API}/run",
            data={
                "ids": ",".join(batch),
                "from": "UniProtKB_AC-ID",
                "to": "Ensembl_Genomes",
            },
        )
        if not resp.ok:
            print(f"UniProt API error {resp.status_code}: {resp.text}")
        resp.raise_for_status()
        job_id = resp.json()["jobId"]

        # Poll until the job finishes
        while True:
            status_resp = requests.get(f"{UNIPROT_API}/status/{job_id}")
            status_resp.raise_for_status()
            status_data = status_resp.json()
            if "results" in status_data or "failedIds" in status_data:
                break
            if status_data.get("jobStatus") == "FINISHED":
                break
            time.sleep(2)

        # Fetch results, following pagination if present
        url = f"{UNIPROT_API}/results/{job_id}?size=500"
        while url:
            results_resp = requests.get(url)
            results_resp.raise_for_status()
            data = results_resp.json()

            for entry in data.get("results", []):
                from_acc = entry["from"]
                to_val = entry["to"]
                # EnsemblPlants entries can come back as a dict
                # ({"geneId": ..., "transcriptId": ...}) or a plain string,
                # depending on API version -- handle both.
                if isinstance(to_val, dict):
                    rap_id = to_val.get("geneId") or to_val.get("transcriptId")
                else:
                    rap_id = str(to_val)
                if rap_id:
                    all_results.append((from_acc, rap_id))

            # Look for a "next" link in the Link header for pagination
            next_url = None
            link_header = results_resp.headers.get("Link", "")
            for part in link_header.split(","):
                if 'rel="next"' in part:
                    next_url = part[part.find("<") + 1 : part.find(">")]
                    break
            url = next_url

        print(f"  -> {len(all_results)} total mappings so far")

    return all_results


def main(network_file, rap_msu_file, output_pickle):
    accessions = collect_accessions(network_file)
    rap_to_msu = load_rap_to_msu(rap_msu_file)

    uniprot_to_rap = submit_id_mapping_job(accessions)

    uniprot_to_msu = {}
    unmapped = 0
    for uniprot_acc, rap_id in uniprot_to_rap:
        msu_genes = rap_to_msu.get(rap_id)
        if not msu_genes:
            unmapped += 1
            continue
        uniprot_to_msu.setdefault(uniprot_acc, [])
        for gene in msu_genes:
            if gene not in uniprot_to_msu[uniprot_acc]:
                uniprot_to_msu[uniprot_acc].append(gene)

    print(f"Mapped {len(uniprot_to_msu)} UniProt accessions to MSU gene IDs")
    print(f"({unmapped} RAP IDs had no MSU crosswalk entry)")

    os.makedirs(os.path.dirname(output_pickle), exist_ok=True)
    with open(output_pickle, "wb") as f:
        pickle.dump(uniprot_to_msu, f)
    print(f"Wrote {output_pickle}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(
            "Usage: python3 build_uniprot_to_msu.py "
            "<STRING-Physical.txt> <RAP-MSU.txt> <output_pickle_path>"
        )
        sys.exit(1)

    main(sys.argv[1], sys.argv[2], sys.argv[3])
