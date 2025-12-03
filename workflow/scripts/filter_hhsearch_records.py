from bioformats import parse_hhr

hhr_files = snakemake.input['hhr']
names_files = snakemake.input['names']
ffindex_files = snakemake.input['ffindex']

output_file = str(snakemake.output)

probab_threshold = snakemake.params['probab']
prefix = snakemake.params['prefix']

def write(input_fh, output_fh, names, prefix):
    for line in input_fh:
        name = line.split('\t', maxsplit = 1)[0]
        if name in names:
            output_fh.write(prefix + line)

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

with open(output_file, 'w') as out_fh:
    for hhr_file, names_file, ffindex_file in zip(hhr_files, names_files, ffindex_files):
        numbers = filter_hhsearch(hhr_file, names_file, probab_threshold)
        with open(ffindex_file) as ffindex_fh:
            write(ffindex_fh, out_fh, numbers, prefix)
