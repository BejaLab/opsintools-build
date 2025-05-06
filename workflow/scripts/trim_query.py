from Bio import SeqIO, PDB

# Assign values from snakemake files
query = snakemake.wildcards.query
query_pdb = snakemake.input.query_p
query_fasta = snakemake.input.query_f
aln_file = snakemake.input.aln_file
trimmed_p = snakemake.output.trim_p
trimmed_f = snakemake.output.trim_f
max_query_len = snakemake.params.max_query_len
query_padding = snakemake.params.query_padding

# Define function for trimming the query
def trim_pdb(structure, output_pdb, start, stop):
    io = PDB.PDBIO()
    io.set_structure(structure)
    class TrimmedStructure(PDB.Select):
          def accept_residue(self, residue):
               res_index = residue.id[1] - 1
               return start <= res_index < stop
    io.save(output_pdb, TrimmedStructure())

# Parse query PDB, parse fasta and get length of sequence
parser = PDB.PDBParser(QUIET = True)
structure = parser.get_structure(query, query_pdb)
record = SeqIO.read(query_fasta, "fasta")
query_len = len(record.seq)

start = query_len
stop = 0
# Find earliest start and latest stop across the aligned files
with open(aln_file, 'r') as f:
    aln = list(SeqIO.parse(f, "fasta"))
    
    # Create an indexed dictionary for alinged sequences
    aln_dict = {i: (rep, query) for i, (rep, query) in enumerate(zip(aln[0].seq, aln[1].seq))}

    # Save first and last positions for query
    aln_start = next((pos for pos, (rep, _) in aln_dict.items() if rep != "-"), None)
    aln_stop = next((pos for pos, (rep, _) in reversed(aln_dict.items()) if rep != "-"), None)

    # Align start and stop index with the clean query sequence
    query_nogap = 0
    query_start = query_stop = None
    for pos, (rep, query) in aln_dict.items():
        if query != "-":
            query_nogap += 1
            if pos >= aln_start and query_start is None:
                query_start = query_nogap - 1
            if pos >= aln_stop:
                query_stop = query_nogap - 1
                break
    
    # Save the lowest start and highest stop
    start = min(start, query_start)
    stop = max(stop, query_stop)

# If there is aligned poistion add padding and trim prot
start = max(start - query_padding, 0)
stop = min(stop + query_padding, query_len)
    
# Trim pdb
trim_pdb(structure, trimmed_p, start, stop)

# Trim fasta
trimmed_seq = record.seq[start:stop]
trimmed_record = SeqIO.SeqRecord(trimmed_seq, id=record.id, description="")
with open(trimmed_f, 'w') as out_f:
    SeqIO.write(trimmed_record, out_f, "fasta")
