import os
from pathlib import Path
import math

num_records_file = snakemake.input['num_records']

chunk = int(snakemake.wildcards['chunk'])
chunk_of = int(snakemake.wildcards['of'])

input_prefix = snakemake.params['input_prefix']
output_prefix = snakemake.params['output_prefix']

a3m_ffindex = f"{input_prefix}_a3m.ffindex"
cs219_ffindex = f"{input_prefix}_cs219.ffindex"
hhm_ffindex = f"{input_prefix}_hhm.ffindex"

a3m_ffdata = f"{input_prefix}_a3m.ffdata"
cs219_ffdata = f"{input_prefix}_cs219.ffdata"
hhm_ffdata = f"{input_prefix}_hhm.ffdata"

a3m_chunk_ffindex = f"{output_prefix}_a3m.ffindex"
cs219_chunk_ffindex = f"{output_prefix}_cs219.ffindex"
hhm_chunk_ffindex = f"{output_prefix}_hhm.ffindex"

a3m_chunk_ffdata = f"{output_prefix}_a3m.ffdata"
cs219_chunk_ffdata = f"{output_prefix}_cs219.ffdata"
hhm_chunk_ffdata = f"{output_prefix}_hhm.ffdata"

Path(a3m_chunk_ffdata).parent.mkdir()

def rel_symlink(source, target):
    target_path = Path(target)
    target_path.symlink_to(os.path.relpath(source, target_path.parent))

rel_symlink(a3m_ffdata,   a3m_chunk_ffdata)
rel_symlink(cs219_ffdata, cs219_chunk_ffdata)
rel_symlink(hhm_ffdata,   hhm_chunk_ffdata)

with open(num_records_file) as f:
    num_records = int(f.read())

chunk_size = math.ceil(num_records / chunk_of)
chunk_start = chunk * chunk_size
chunk_stop = chunk_start + chunk_size

with open(a3m_ffindex) as a3m_in, open(cs219_ffindex) as cs219_in, open(hhm_ffindex) as hhm_in:
    with open(a3m_chunk_ffindex, 'w') as a3m_out, open(cs219_chunk_ffindex, 'w') as cs219_out, open(hhm_chunk_ffindex, 'w') as hhm_out:
        hhm_name = hhm_line = ''
        for line_num in range(num_records):
            if line_num == chunk_stop:
                break
            a3m_line = next(a3m_in)
            cs219_line = next(cs219_in)
            if line_num >= chunk_start:
                a3m_out.write(a3m_line)
                cs219_out.write(cs219_line)
                name = a3m_line.split('\t', maxsplit = 1)[0]
                if hhm_name < name:
                    for hhm_line in hhm_in:
                        hhm_name = hhm_line.split('\t', maxsplit = 1)[0]
                        if hhm_name >= name:
                            break
                if hhm_name == name:
                    hhm_out.write(hhm_line)
