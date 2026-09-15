import csv
import os
import pickle


def convert_msu(msu_id_file, mapping_file, output_name, skip_no_matches):
    output_dir = os.path.dirname(output_name)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(msu_id_file) as msu_file, open(mapping_file, "rb") as mapping, open(output_name, "w") as output:
        mapping_dict = pickle.load(mapping)

        csv_reader = csv.reader(msu_file, delimiter="\t")
        for line in csv_reader:
            output_set = set()
            for msu_id in line:
                if len(mapping_dict[msu_id]) != 0:
                    output_set = output_set.union(mapping_dict[msu_id])

            output.write("\t".join(list(output_set)))

            if skip_no_matches and len(output_set) > 0:
                output.write("\n")
            elif not skip_no_matches:
                output.write("\n")

    print(f"Generated {output_name}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "msu_module_file",
        help="text file containing the list of MSU accessions to be converted",
    )
    parser.add_argument(
        "mapping_file",
        help="pickled dictionary mapping the MSU accessions to the target IDs",
    )
    parser.add_argument(
        "output_name",
        help="output filename (including directory) for the file containing the equivalent IDs after conversion",
    )
    parser.add_argument(
        "--skip_no_matches",
        action="store_true",
        help="accessions that cannot be converted will be skipped",
    )

    args = parser.parse_args()

    convert_msu(
        args.msu_module_file,
        args.mapping_file,
        args.output_name,
        args.skip_no_matches,
    )
