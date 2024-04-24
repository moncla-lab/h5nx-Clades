import json
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--branch_lengths', type=str, help='path to branch lengths JSON produced by augur refine')
parser.add_argument('--output_branch_lengths', type=str, help='path to save branch lengths JSON without time information')
args = parser.parse_args()
file = args.branch_lengths
outfile = args.output_branch_lengths

with open(file,'r') as f:
    d = json.load(f)
    for k in d['nodes']:
        d['nodes'][k] = {'branch_length':d['nodes'][k]['mutation_length']}
        
with open(outfile,'w') as f:
    json.dump(d, f)