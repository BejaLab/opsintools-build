from Bio import SeqIO
from Bio.Seq import Seq
import re

input_a2m_file = str(snakemake.input)
output_a2m_file = snakemake.output['a2m']
output_a3m_file = snakemake.output['a3m']

records = SeqIO.parse(input_a2m_file, "fasta")
ss_pred = next(records)
ss_conf = next(records)
consensus = next(records)
assert ss_pred.id == "ss_pred", f"Expected ss_pred, got {ss_pred.id}"
assert ss_conf.id == "ss_conf", f"Expected ss_conf, got {ss_conf.id}"
assert consensus.id.endswith("_consensus"), f"Expected *_consensus, got {consensus.id}"
i, j = ss_pred.seq.find('H'), ss_pred.seq.rfind('H') + 1

def write_a2m(record, fh):
    SeqIO.write(record, fh, "fasta")

def write_a3m(record, fh):
    seq = str(record.seq).replace('.', '') 
    seq = re.sub(r'^[a-z]+', '', seq)
    seq = re.sub(r'[a-z]+$', '', seq)
    record.seq = Seq(seq)
    SeqIO.write(record, fh, "fasta")

with open(output_a2m_file, 'w') as a2m, open(output_a3m_file, 'w') as a3m:
    ss_pred.seq = ss_pred.seq[i:j]
    ss_conf.seq = ss_conf.seq[i:j]
    consensus.seq = consensus.seq[i:j]
    write_a3m(ss_pred, a3m)
    write_a3m(ss_conf, a3m)
    write_a3m(consensus, a3m)
    for record in records:
        record.seq = record.seq[i:j]
        write_a2m(record, a2m)
        write_a3m(record, a3m)
