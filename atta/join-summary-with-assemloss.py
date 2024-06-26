#! /usr/bin/env python
import pandas as pd
import sys
import argparse
import glob


def get_single_glob_match(pattern):
    matches = glob.glob(pattern)
    if not matches:
        print(f"No matches to glob pattern '{pattern}'; are files missing?!",
              file=sys.stderr)
        raise Exception
    elif len(matches) > 1:
        print(f"Too many matches to glob pattern '{pattern}'; what's up?",
              file=sys.stderr)
        print(matches)
        raise Exception
    else:
        assert len(matches) == 1
        return matches[0]


def main():
    p = argparse.ArgumentParser()
    p.add_argument('grist_summary_csv')
    p.add_argument('atta_directory')
    p.add_argument('-o', '--output-csv', required=True)
    args = p.parse_args()

    summary = pd.read_csv(args.grist_summary_csv)

    accessions = list(summary['genome_id'])
    sample_id = set(summary['sample_id'])
    assert len(sample_id) == 1
    sample_id = list(sample_id)[0]

    print(f'sample_id: {sample_id}')
    print(f'got {len(accessions)} genomes.')

    all_accs = []
    for acc in accessions:
        print(f'globbing on {acc}')
        base_pattern = f'{args.atta_directory}/{sample_id}.x.{acc}.*'
        trim_pattern = base_pattern + '.trim.csv'
        assem_pattern = base_pattern + '.assem.csv'
        mapassem_pattern = base_pattern + '.mapassem.csv'

        trim_csv = get_single_glob_match(trim_pattern)
        assem_csv = get_single_glob_match(assem_pattern)
        mapassem_csv = get_single_glob_match(mapassem_pattern)

        try:
            df1 = pd.read_csv(trim_csv)
            df2 = pd.read_csv(assem_csv)
            df3 = pd.read_csv(mapassem_csv)
        except pd.errors.EmptyDataError:
            print(f'SKIPPING {sample_id} {acc}')

        assert len(df1) == 1
        assert len(df2) == 1
        assert len(df3) == 1

        trim_cont = list(df1.similarity)[0]
        assem_cont = list(df2.similarity)[0]
        mapassem_cont = list(df3.similarity)[0]

        d = dict(genome_id=acc, sample_id=sample_id, trim_cont=trim_cont,
                 assem_cont=assem_cont, mapassem_cont=mapassem_cont)
        all_accs.append(d)

    # turn into dataframe
    df = pd.DataFrame.from_dict(all_accs)

    # join!
    join_df = summary.merge(df, on=['genome_id', 'sample_id'], how='left') 
                                                                
    # write!
    join_df.to_csv(args.output_csv)
    

if __name__ == '__main__':
    sys.exit(main())
