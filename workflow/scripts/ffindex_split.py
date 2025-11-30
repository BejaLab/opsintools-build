import math

a3m_input = snakemake.input['a3m']
cs219_input = snakemake.input['cs219']
hhm_input = snakemake.input['hhm']

a3m_outputs = snakemake.output['a3m']
cs219_outputs = snakemake.output['cs219']
hhm_outputs = snakemake.output['hhm']

with open(a3m_input) as file:
    num_records = sum(1 for _ in file)

num_chunks = len(a3m_outputs)
chunk_size = math.ceil(num_records / num_chunks)
records_iter = iter(range(num_records))

with open(a3m_input) as a3m_in, open(cs219_input) as cs219_in, open(hhm_input) as hhm_in:
    for chunk, (a3m_output, cs219_output, hhm_output) in enumerate(zip(a3m_outputs, cs219_outputs, hhm_outputs)):
        stop = chunk_size * (chunk + 1) - 1
        with open(a3m_output, 'w') as a3m_out, open(cs219_output, 'w') as cs219_out, open(hhm_output, 'w') as hhm_out:
            hhm_number = hhm_line = ''
            offset = 0
            for record_num in records_iter:
                a3m_line = next(a3m_in)
                cs219_line = next(cs219_in)
                a3m_out.write(a3m_line)
                cs219_out.write(cs219_line)
                number, start, length = a3m_line.rstrip().split('\t')
                if hhm_number < number:
                    for hhm_line in hhm_in:
                        hhm_number = hhm_line.split('\t', maxsplit = 1)[0]
                        if hhm_number >= number:
                            break
                if hhm_number == number:
                    hhm_out.write(hhm_line)
                if record_num == stop:
                    break
