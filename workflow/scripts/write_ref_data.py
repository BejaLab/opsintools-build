import re, json

pdb_file = snakemake.input['pdb']
opm_file = snakemake.input['opm']
output_file = str(snakemake.output)

def parse_opm(opm_file):
    with open(opm_file) as file:
        data = json.load(file)
    matches = re.findall('(\\d+)\\(\\s*(\\d+)-\\s*(\\d+)\\)', data['subunits'][0]['segment'])
    return { f"TM{tm}": [ int(start), int(end) ] for tm, start, end in matches }

def get_res_pos(pdb_file):
    res_pos = {}
    with open(pdb_file) as file:
        for line in file:
            if line.startswith('ATOM'):
                res_pos[int(line[22:26])] = 1
    return list(res_pos.keys())

data = {
    "id": snakemake.params['ref_pdb'],
    "lysine": snakemake.params['ref_lysine_pos'],
    "res_pos": get_res_pos(pdb_file),
    "tms": parse_opm(opm_file)
}
with open(output_file, 'w') as file:
    json.dump(data, file)

