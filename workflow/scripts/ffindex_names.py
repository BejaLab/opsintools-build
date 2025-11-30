index_file = snakemake.input['index']
data_file = snakemake.input['data']
names_file = str(snakemake.output)

max_name_len = snakemake.params['max_name_len']

records = {}
with open(index_file) as index:
    for line in index:
        number, start, length = line.split()
        records[number] = int(start) + 1, min(max_name_len, int(length) - 1)

with open(data_file, 'rb') as data, open(names_file, 'w') as out:
    for number, (start, length) in sorted(records.items(), key = lambda x: x[1][0]):
        data.seek(start)
        header = data.read(length)
        name = header.split(maxsplit = 1)[0].decode("utf-8")
        out.write(f"{name}\t{number}\n")
