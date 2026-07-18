"""
Builds static/app_data/gene_id_mapping/msu_mapping/msu_to_rap.pickle
by inverting the RAP->MSU crosswalk (RAP-MSU.txt) we already downloaded.

Usage:
    python3 build_msu_to_rap.py /tmp/RAP-MSU.txt static/app_data/gene_id_mapping/msu_mapping/msu_to_rap.pickle
"""

import os
import pickle
import sys


def main(rap_msu_file, output_pickle):
    msu_to_rap = {}

    with open(rap_msu_file) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 2:
                continue
            rap_id, msu_field = parts
            if msu_field == "None":
                continue
            for msu_transcript in msu_field.split(","):
                gene = msu_transcript.split(".")[0]
                if not gene:
                    continue
                # Keep the first RAP ID seen for a given MSU gene
                # (a gene can rarely appear under more than one RAP locus)
                msu_to_rap.setdefault(gene, rap_id)

    print(f"Built {len(msu_to_rap)} MSU->RAP mappings")

    os.makedirs(os.path.dirname(output_pickle), exist_ok=True)
    with open(output_pickle, "wb") as out:
        pickle.dump(msu_to_rap, out)
    print(f"Wrote {output_pickle}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 build_msu_to_rap.py <RAP-MSU.txt> <output_pickle_path>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
