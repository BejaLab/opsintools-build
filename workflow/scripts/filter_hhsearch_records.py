from bioformats import parse_hhr

hhr_files = snakemake.input['hhr']

names_files = snakemake.input['names']
a3m_inputs = snakemake.input['a3m']
cs219_inputs = snakemake.input['cs219']
hhm_inputs = snakemake.input['hhm']

a3m_output = snakemake.output['a3m']
cs219_output = snakemake.output['cs219']
hhm_output = snakemake.output['hhm']

probab_threshold = snakemake.params['probab']

def write(input_file, output_fh, names):
    with open(input_file) as input_fh:
        for line in input_fh:
            name = line.split('\t', maxsplit = 1)[0]
            if name in names:
                output_fh.write(line)

def filter_hhsearch(hhr_file, names_file, probab_threshold):
    numbers = set()
    names = set()
    with open(hhr_file) as file:
        for record in parse_hhr.read_hhr(file):
            if record['Probab'] >= probab_threshold:
                names.add(record['template']['name'])
    with open(names_file) as file:
        for line in file:
            name, number = line.split()
            if name in names:
                numbers.add(number)
    return numbers

with open(a3m_output, 'w') as a3m_fh, open(cs219_output, 'w') as cs219_fh, open(hhm_output, 'w') as hhm_fh:
    for hhr_file, names_file, a3m_input, cs219_input, hhm_input in zip(hhr_files, names_files, a3m_inputs, cs219_inputs, hhm_inputs):
        numbers = filter_hhsearch(hhr_file, names_file, probab_threshold)
        write(a3m_input,   a3m_fh,   numbers)
        write(cs219_input, cs219_fh, numbers)
        write(hhm_input,   hhm_fh,   numbers)
