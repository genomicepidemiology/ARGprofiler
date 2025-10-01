#!/bin/python3

import os
import sys
import argparse
import pandas as pd
import subprocess
import re
from collections import defaultdict
import numpy as np
import glob

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
        '-n', '--name',
        type=str,
        required=True,
        dest='name'
    )

    parser.add_argument(
        '--run',
        action='store_true',
        help='Run the SQL command',
        dest='run'
    )

    return parser.parse_args()

def unpack(fname):
    p = subprocess.run(f"gzip -dc {fname}", shell=True, stdout=subprocess.PIPE)
    o = p.stdout.decode().split('\n')

    data = []
    for l in o:
        ll = l.split('\t')
        if len(ll) == 7:
            data.append(ll[1:])

    df = pd.DataFrame(data, columns=['equally_well_mapping_templates', 'alignment_score', 'start', 'end', 'template', 'query'])
    df = df.convert_dtypes()

    df['db'] = 'PanRes'
    run_id = os.path.basename(fname).split('.')[0]
    df['run_accession'] = run_id

    to_convert = ['equally_well_mapping_templates', 'alignment_score', 'start', 'end']
    if fname.endswith('frag.gz'):
        for c in to_convert:
            df[c] = df[c].astype(int)
        db_version = df['template'].values[0].split('_')[-1]
        df['db_version'] = db_version
        df['template'] = df['template'].str.replace('_' + db_version, '')
        df['contig_length'] = df['query'].str.extract(r"length_(\d+)_cov")
    elif fname.endswith('frag_raw.gz'):
        for c in to_convert + ['template']:
            df[c] = df[c].apply(lambda x: [int(v.strip()) for v in x.split(',')] if ',' in x else int(x))


    return df, run_id


def _convert(val):
    if isinstance(val, int):
        v = [str(val)]
    elif isinstance(val, list):
        v = [str(v) for v in val]
    elif isinstance(val, float) or val.isna().any(): #np.isnan(val):
        v =['']

    return ",".join(v)

def to_sql(row, out):


    start_multi = _convert(row['start_multi'])
    end_multi = _convert(row['end_multi'])
    multi_names = row['template_multi_name']



    # CALL InsertARGextender('run_accession','db','db_version','contig','contig_length','target','start_template','end_template','targets','start_templates','end_templates',)
    cmd="CALL InsertARGextender('{}', '{}', '{}', '{}', '{}', '{}', '{}', '{}', '{}', '{}', '{}');".format(
        row['run_accession'], row['db'], row['db_version'], 
        row['query'], row['contig_length'],
        row['template'], row['start'], row['end'],
        multi_names, start_multi, end_multi
    )

    print(cmd, file=out)

def get_panres_names(namefile):
    idx2name = {np.nan: ''}

    with open(namefile, 'r') as f:
        for i, l in enumerate(f.readlines()):
            idx2name[i+1] = "_".join(l.strip().split('_')[:-1])

    return idx2name

if __name__ == '__main__':
    args = parse_args()

    df, runid = unpack(args.file)
    fdf, _ = unpack(args.file.replace('frag', 'frag_raw'))
    
    df = df.merge(fdf, on=['query', 'equally_well_mapping_templates', 'db'], how='left', suffixes=('', '_multi'))

    idx2pan = get_panres_names(args.name)

    df['template_multi_name'] = df['template_multi'].apply(
        lambda x: ",".join([idx2pan[abs(i)] for i in x] if isinstance(x, list) else [idx2pan[x]])
    )

    outFile = os.path.join(args.output, 'argextender_' + runid + '.sql')
    print(outFile)
    out = open(outFile, 'w')
    df.apply(lambda x: to_sql(x, out), axis=1)

    out.close()

    if args.run:
        cmd = f"mysql --database=AvA_2 < {outFile}"
        p = subprocess.run(cmd, shell=True)
        print(p) 
