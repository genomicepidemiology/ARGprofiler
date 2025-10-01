import argparse
import subprocess

def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        '--run_accession',
        type=str,
        help='Run accession',
        required=True,
        dest='run_accession'
    )

    parser.add_argument(
        '--rule',
        type=str,
        help='Rule name',
        required=True,
        dest='rule'
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

if __name__ == '__main__':

    args = parse_args()
    
    cmd = f"mysql -e \"use AvA_2; INSERT INTO pipeline_process VALUES('{args.run_accession}', '{args.rule}', 0) ON DUPLICATE KEY UPDATE status='0';\""
    subprocess.run(cmd, shell=True)
