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
            truth[values[0]] = int(values[1])
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
            taxa.add(int(line))
    return taxa


def count_matches(res, truth):
    match = 0
    for key, value in res.items():
        if value == truth.get(key):
            match += 1
    return match


def filter_results(res, taxa):
    new_res = {}
    for key, value in res.items():
        if value in taxa:
            new_res[key] = value
    return new_res


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

    total_found_reads = len(res)
    total_selected_taxa = len(taxa)
    total_true_reads = len(truth)
    filtered_res = filter_results(res, taxa)
    total_filtered_found_reads = len(filtered_res)
    matches = count_matches(res, truth)
    filtered_matches = count_matches(filtered_res, truth)

    print(f"Number of true reads: {total_true_reads}, unique true taxa: {len(set(truth.values()))}")
    print(f'Number of found reads: {total_found_reads}, unique found taxa: {len(set(res.values()))}')
    print(f'Number of matches: {matches} ({matches / total_found_reads * 100:.2f}%)')
    print()
    print(f'Input taxa: {total_selected_taxa}')
    print(f'Filtered found reads: {total_filtered_found_reads}, unique filtered found taxa: {len(set(filtered_res.values()))}')
    print(f'Filtered matches: {filtered_matches} ({filtered_matches / total_filtered_found_reads * 100:.2f}%)')


if __name__ == '__main__':
    main()
