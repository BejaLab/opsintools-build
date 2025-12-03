import re, json
from Bio import SeqIO

pdb_file = snakemake.input['pdb']
opm_file = snakemake.input['opm']
fasta_file = snakemake.input['fasta']
output_file = str(snakemake.output)

seq = SeqIO.read(fasta_file, 'fasta')

def parse_opm(opm_file):
    with open(opm_file) as file:
        data = json.load(file)
    matches = re.findall('(\\d+)\\(\\s*(\\d+)-\\s*(\\d+)\\)', data['subunits'][0]['segment'])
    return { f"TM{tm}": [ int(start), int(end) ] for tm, start, end in matches }

data = {
    "id": snakemake.params['ref_pdb'],
    "lysine": snakemake.params['ref_lysine_pos'],
    "tms": parse_opm(opm_file),
    "seq": str(seq.seq)
}
with open(output_file, 'w') as file:
    json.dump(data, file)
