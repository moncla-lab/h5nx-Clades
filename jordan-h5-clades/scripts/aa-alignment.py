from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--nt_alignment', type=str, help='path to nt alignment fasta file')
parser.add_argument('--reference', type=str, help='path to reference genbank file')
parser.add_argument('--aa_alignment', type=str, help='path to output aa alignment fasta file')

args = parser.parse_args()

nt_alignment = args.nt_alignment
reference = args.reference
aa_alignment = args.aa_alignment


# determine CDS start and end from reference
ref_record = SeqIO.read(reference, 'gb')
cds_features = [feat for feat in ref_record.features if feat.type=='CDS']

if len(cds_features) == 1:
    cds_location = cds_features[0].location
    cds_start = cds_location.start
    cds_end = cds_location.end

elif len(cds_features) == 0:
    print('ERROR: no CDS feature found in reference')

elif len(cds_features) > 1:
    print('ERROR: multiple CDS features found in reference')


# translate CDS from aligned nt sequences
records = []

for record in SeqIO.parse(nt_alignment, 'fasta'):
    cds = record.seq.replace('-','n')[cds_start:cds_end]
    record.seq = cds.translate()
    records.append(record)

# write to output fasta
SeqIO.write(records, aa_alignment, 'fasta')