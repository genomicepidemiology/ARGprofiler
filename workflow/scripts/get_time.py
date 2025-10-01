import os
import argparse
import pandas as pd
import subprocess
import re

def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        '-f', '--file',
        type=str,
        help='Bench file',
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
        '--run',
        action='store_true',
        help='Run the SQL command',
        dest='run'
    )
    return parser.parse_args()

folder2rule = {
    'raw_reads': 'download_{}_reads',
    'trimmed_reads': 'trim_{}_reads',
    'kma_mOTUs': 'kma_{}_reads_mOTUs',
    'kma_panres': 'kma_{}_reads_panRes',
    'ARG_extender': 'ARG_extender_{}_reads',
    'mash_sketch': 'mash_sketch_{}_reads',
}

index2rename={
    'command_being_timed': 'command',
    'percent_of_cpu_this_job_got': 'percent_cpu_job',
    'system_time_seconds': 'system_time_sec',
    'elapsed_wall_clock_time_h:mm:ss_or_m:ss': 'elapsed_wallclock_time_sec',
    'major_requiring_i/o_page_faults': 'major_page_faults',
    'minor_reclaiming_a_frame_page_faults': 'minor_page_faults'
}

def read_bench(f):
    o = pd.read_csv(f, sep=': ', header=None, index_col=0)

    indices = [i.replace(' ', '_').replace('(', '').replace(')', '').lower() for i in o.index]
    o.index = indices
    o = o.rename(index=index2rename)

    run_accession = os.path.basename(f).split('.')[0]
    o.loc['run_accession', 1] = run_accession

    folders = f.replace(os.path.basename(f), '').split(os.sep)
    o.loc['rule', 1] = folder2rule[folders[1]].format(folders[2])
    o.loc['command', 1] = o.loc['command', 1].replace('"', '')
    
    return o

def convert_time(extracted_time):
    p = re.compile(r"(\d+:\d{2}:\d{2}|\d{1,2}:\d{2}\.\d{2})")

    if extracted_time.count(':') == 1:
        m, s = map(float, extracted_time.split(':'))
        h = 0
    elif extracted_time.count(':') == 2:
        h, m, s = map(int, extracted_time.split(':'))

    total_seconds = h * 3600 + m * 60 + s
    return total_seconds

def read_bench2(fname):
    indexes, values = [], []
    with open(fname, 'r') as f:
        #lines = [l.strip() for l in f.readlines()]
        for l in f.readlines():
            l = l.strip()
            ls = l.split(': ')
            indexes.append(ls[0].replace(' ', '_').replace('(','').replace(')', '').lower())
            values.append(ls[1].strip().replace("%", ""))

    o = pd.DataFrame(values, index=indexes)
    o = o.rename(index=index2rename)

    run_accession = os.path.basename(fname).split('.')[0]
    o.loc['run_accession', 0] = run_accession

    folders = fname.replace(os.path.basename(fname), '').split(os.sep)
    o.loc['rule', 0] = folder2rule[folders[1]].format(folders[2])
    o.loc['command', 0] = o.loc['command', 0].replace('"', '')


    o.loc['elapsed_wallclock_time_sec', 0] = convert_time(o.loc['elapsed_wallclock_time_sec', 0])

    return o


def parse_bench(df, out_dir):
    columns = ['run_accession','rule', 'exit_status','command','user_time_seconds','system_time_sec','percent_cpu_job','elapsed_wallclock_time_sec','average_shared_text_size_kbytes','average_unshared_data_size_kbytes','average_stack_size_kbytes','average_total_size_kbytes','maximum_resident_set_size_kbytes','average_resident_set_size_kbytes','major_page_faults','minor_page_faults','voluntary_context_switches','involuntary_context_switches','swaps','file_system_inputs','file_system_outputs','socket_messages_sent','socket_messages_received','signals_delivered','page_size_bytes']
    cmd = "CALL InsertRuleStatus({});".format(",".join(df.loc[columns, 0].apply(lambda x: f"\'{x}\'").values.tolist()))

    outName = "_".join(df.loc[columns[:2], 0].values.tolist()[::-1]) + '.sql'
    outFile = os.path.join(out_dir, outName)
    with open(outFile, 'w') as f:
        print(cmd, file=f)

    return outFile


if __name__ == "__main__":
    args = parse_args()
    os.makedirs(args.output, exist_ok=True)
    df = read_bench2(args.file)
    outFile = parse_bench(df, out_dir=args.output)

    if args.run:
        p = subprocess.run(f"mysql --database=AvA_2 < {outFile}", shell=True)



