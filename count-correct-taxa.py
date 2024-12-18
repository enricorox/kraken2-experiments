import argparse


def read_ground_truth_file(filename):
    """
    Read from file like
    S0R0	1292	Staphylococcus warneri
    or
    taxid_1042876.1	1042876

    and save tax id associated to read id
    :param filename: tsv file to read
    :return: a map between read id and tax id
    """
    truth = {}
    with open(filename, 'r') as truth_file:
        for line in truth_file:
            values = line.split('\t')
            truth[values[0]] = values[1]
    return truth


def read_results_file(filename):
    """
    Read from file like
    SRR1804065.1	816

    :param filename: tsv file to read
    :return: a map between read id and tax id
    """
    return read_ground_truth_file(filename)


def read_taxa_file(filename):
    taxa = set()
    with open(filename, 'r') as taxa_file:
        for line in taxa_file:
            taxa.add(line)
    return taxa


def compare(res, taxa, truth):
    match = 0
    selected_match = 0
    for key, value in res.items():
        if value == truth.get(key):
            match += 1
            if value in taxa:
                selected_match += 1
    return match, selected_match


def main():
    # define parser
    parser = argparse.ArgumentParser(description="Count the number of correct taxa")
    parser.add_argument("--res", type=str, help="Kraken2 result file", required=True)
    parser.add_argument("--taxa", type=str, help="file containing selected taxa to compare", required=True)
    parser.add_argument("--truth", type=str, help="ground truth file", required=True)

    # parse arguments
    args = parser.parse_args()

    # read input files
    truth = read_ground_truth_file(args.truth)
    res = read_results_file(args.res)
    taxa = read_taxa_file(args.taxa)

    total_queries = len(res)
    total_taxa = len(taxa)
    total_taxa_in_queries = len(taxa & set(res.values()))
    match, selected_math = compare(res, taxa, truth)

    print(f'Number of queries: {total_queries}')
    print(f'Number of matches: {match} ({match / total_queries * 100:.2f}%)')
    print()
    print(f'Input taxa: {total_taxa}')
    print(f'Selected taxa: {total_taxa_in_queries}')
    print(f'Selected taxa matches: {selected_math} ({selected_math / total_taxa_in_queries * 100:.2f}%)')


if __name__ == '__main__':
    main()
