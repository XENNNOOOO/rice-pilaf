NODE_A = 0
NODE_B = 1
SCORE = 2
SEPARATOR = " "
SPECIES_PREFIX_SEPARATOR = "."
HAS_HEADER = True

def extract_interactions(STRING_file):
    all_interactions = []
    
    in_header = HAS_HEADER
    
    with open(STRING_file) as network:
        for row in network:
            row = row.rstrip()
            col = row.split(SEPARATOR)
            
            if in_header:
                in_header = False
                continue
            
            values = []
            values.append(remove_prefix(col[NODE_A], SPECIES_PREFIX_SEPARATOR))
            values.append(remove_prefix(col[NODE_B], SPECIES_PREFIX_SEPARATOR))
            values.append(col[SCORE])
            
            all_interactions.append(values)  
            
    return all_interactions

def remove_prefix(string, separator):
    return string.split(separator)[1]

def display(list):
    first_ten = list[:10]
    print(first_ten)

def remove_header(list):
    return list[1:]

def output_to_file(list, output_dir):
    with open(f"{output_dir}", "w") as f:
        for row in list:
            f.write("\t".join(row))
            f.write("\n")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "STRING_file",
        help="the downloaded PPI file from STRING version 12.0.",
    )
    """
    example of the file's contents:

    protein1 protein2 combined_score
    39947.A0A075DNI2 39947.B9FNL5 187
    39947.A0A075DNI2 39947.Q6H432 293
    39947.A0A075DNI2 39947.A0A0P0VQJ4 182
    39947.A0A075DNI2 39947.A0A0P0VJ93 223
    """

    parser.add_argument(
        "output_dir", help="output directory for the converted module list"
    )
    args = parser.parse_args()

    result = extract_interactions(args.STRING_file)
    output_to_file(result, args.output_dir)