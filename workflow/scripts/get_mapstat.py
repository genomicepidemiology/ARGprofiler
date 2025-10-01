#!/bin/python3

import os
import sys
import argparse
import pandas as pd
import subprocess
import re
from collections import defaultdict
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        '-f', '--file',
        type=str,
        help='mapstat file',
        required=True,
        dest='file'
    )

    parser.add_argument(
        '-o', '--output',
        type=str,
        default='results/sql/',
        dest='output'
    )


    parser.add_argument(
        '--motus',
        action='store_true',
        help = '',
        dest='motus'
    )

    parser.add_argument(
        '--pad_taxa',
        action='store_true',
        help='If some taxonmic level is unknown, pad it with the Highest level + id',
        dest='pad_taxa'
    )

    parser.add_argument(
        '--run',
        action='store_true',
        help='Run the SQL command',
        dest='run'
    )
    return parser.parse_args()
    

def read_mapstat(f,min_identity=0):
    df = pd.read_csv(f, sep='\t', skiprows=6)

    # remove zero counts
    df = df.loc[df['fragmentCountAln'] > 0,]

    if df.shape[0] == 0: 
        print("No counts are above zero for file:", f)
        sys.exit()


    # calculate read identity
    df['identity'] = df.apply(
        lambda x: (x['bpTotal'] - x['snpSum'] - x['insertSum'] - x['deletionSum']) / (x['bpTotal']), axis = 1
    )

    # filter on read identity
    df = df.loc[df['identity'] >= min_identity, ]

    # rename columns
    df.columns = [c.replace('# ', '') for c in df.columns]
    
    # add sample name
    run_id = os.path.basename(f).split('.')[0]
    df['run_accession'] = run_id
    
    # get header
    p = subprocess.run(f"grep database {f}", shell=True, stdout=subprocess.PIPE)
    db = p.stdout.decode().split('\t')[-1].strip()
    header_data={}
    p = subprocess.run(f"grep '^## ' {f}", shell=True, stdout=subprocess.PIPE)
    for l in p.stdout.decode().split('\n'):
        ls = l.split('\t')
        if len(ls) == 2:
            k = ls[0].replace('## ', '')
            v = ls[1]
            header_data[k] = v.strip()


    db = header_data['database']
    if db == 'panres':
        db = 'PanRes'
        db_version = df['refSequence'].values[0].split('_')[-1]
        df['refSequence'] = df['refSequence'].str.replace('_' + db_version, '')
    elif db == 'db_mOTUs':
        db = 'db_mOTU'
        db_version = '20221205'
    df['db'] = db
    df['db_version'] = db_version

    kma_header = "CALL InsertKMAHeader('{}', '{}', '{}', '{}', '{}', '{}')".format(
        run_id, db, db_version, 
        header_data['version'],
        header_data['date'],
        header_data['fragmentCount']
    )

    return df, run_id, kma_header

def pad_taxonomy2(df, levels):
    levels_extended = [v + '_name' for v in levels]

    id_value = 'unannotated'


    to_drop_rows = []
    
    for rowid, row in df.iterrows():
        need_padding = row[levels_extended].isna()
        if need_padding.all():
            to_drop_rows.append(rowid)
        elif need_padding.sum() > 0:
            for i in np.where(need_padding)[0]:
                to_pad = levels_extended[i]
                # Find the last non-NaN value before the current NaN
                pad_name = None
                for k in range(i-1, -1,-1):
                    if pd.notna(row[levels_extended[k]]):
                        pad_name = row[levels_extended[k]]
                        break

                # if no non-NaN value is found, use first element
                if pad_name is None:
                    pad_name = row[levels_extended[0]]
                if pad_name.endswith(id_value):
                    row[to_pad] = pad_name
                else:
                    row[to_pad] = pad_name + ' ' + id_value
                
            df.loc[rowid] = row

    df = df.loc[~(df.index.isin(to_drop_rows))]
    return df 

def agg_taxa(mapstat, groupCols, aggDict={'fragmentCountAln_adj': 'sum', 'refSequence': 'nunique', 'fragmentCountAln': 'sum'}):
    
    aggData = mapstat.groupby(groupCols).agg(aggDict).reset_index()
    aggData['level'] = groupCols[-1].split('_')[0]
    aggData['path'] = ";".join([c.split('_')[0] for c in groupCols if 'name' in c])
    aggData['path_taxs'] = aggData[[c for c in groupCols if 'tax' in c]].astype(int).astype(str).apply(lambda x: ";".join(x), axis=1)
    aggData.rename(columns={c: c.split('_')[-1] for c in groupCols[-2:]}, inplace=True)
    aggData.drop(columns=groupCols[1:], errors='ignore', inplace=True)

    aggData['tax'] = aggData['tax'].astype(int)

    return aggData

def get_taxonomy(file='~/gs3/data/db_mOTU_taxonomy.tsv'):

    #query = "use taxonomy; select id, gene_len, kingdom_name, kingdom_tax, phylum_name, phylum_tax, class_name, class_tax, order_name, order_tax, family_name, family_tax, genus_name, genus_tax, mOTU_name, mOTU_tax from db_mOTU_20221205"
    #query = "use taxonomy; select * from db_mOTU_20221205"
    #cmd = f"mysql -e \"{query}\""
    #p = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
    #if p.returncode != 0:
    #    print(p)
    #o = p.stdout.decode()
    #taxa = pd.read_csv(StringIO(o), sep='\t')
    taxa = pd.read_csv(file, sep='\t')

    # spaces in some names?
    for strCol in taxa.select_dtypes(include='object'):
        taxa[strCol] = taxa[strCol].str.strip()
    return taxa


def to_sql(mapstat, is_motus, out, run=False):
    
    motus_cmd = lambda x: "CALL InsertmOTUsCounts('{}', '{}', '{}', '{}', '{}', '{}', '{}', '{}', '{}');".format(
        x['run_accession'], x['name'], x['tax'],
        x['level'], x['path'], x['path_taxs'],
        x['n_references'], x['fragmentCountAln'], x['fragmentCountAln_adj']
    )

    panres_cmd = lambda x: "CALL InsertKMAPanRes('{}', '{}','{}','{}','{}','{}','{}','{}','{}','{}','{}','{}','{}','{}','{}','{}','{}','{}');".format(
        x['run_accession'], x['db'], x['db_version'],
        x['refSequence'], x['readCount'], x['fragmentCount'],
        x['mapScoreSum'], x['refCoveredPositions'], x['refConsensusSum'],
        x['bpTotal'], x['depthVariance'], x['nucHighDepthVariance'],
        x['depthMax'], x['snpSum'], x['insertSum'],
        x['deletionSum'], x['readCountAln'], x['fragmentCountAln']
    )

    if is_motus:
        cmds = mapstat.apply(motus_cmd, axis=1).values.tolist()
    else:
        cmds = mapstat.apply(panres_cmd, axis=1).values.tolist()

    print("\n".join(cmds), file=out)
    

    

if __name__ == '__main__':
    args = parse_args()

    if not os.path.exists(args.file) or os.path.getsize(args.file) == 0:
        print(f"{args.file} is either emtpy or does not exist.")
        sys.exit(0)

    df, run_id, kma_header_sql = read_mapstat(args.file)


    db = 'mOTU' if args.motus else 'PanRes'
    outFilename = os.path.join(args.output, f"kma_{run_id}_{db}.sql")  
    out = open(outFilename, 'w')
    print(kma_header_sql, file=out)

    if args.motus:
        taxa = get_taxonomy()
    
        levels = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'mOTU']
        results = defaultdict(list)


        df2 = df.merge(taxa, left_on=['refSequence'], right_on=['refSequence'])
        df2['fragmentCountAln_adj'] = df2['fragmentCountAln'] / (df2['gene_length']/1e3)
        col = ['run_accession']

        if args.pad_taxa:
            df2 = pad_taxonomy2(df2, levels=levels)
        
        for taxlevel in levels:
            col += [taxlevel + '_name', taxlevel + '_tax']
            
            aggData = agg_taxa(mapstat=df2, groupCols=col)
            results[taxlevel].append(aggData)
    
        for k, v in results.items():
            r = pd.concat(v).rename(columns={'refSequence': 'n_references'})
            to_sql(r, is_motus=args.motus, run=args.run, out=out)
    else:
        to_sql(df, is_motus=args.motus, run=args.run, out=out)

    out.close()

    if args.run:
        cmd = f"mysql --database=AvA_2 < {outFilename}"
        p = subprocess.run(cmd, shell=True)
        print(p) 
