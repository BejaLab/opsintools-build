import requests, os, json

output_file = str(snakemake.output)
pdb_id = snakemake.wildcards['pdb']

host = "https://opm-back.cc.lehigh.edu/opm-backend/primary_structures"

# Define a function that gets the opm_id from the search request
def get_opm_id(pdb_id):
    response = requests.get(f"{host}?search={pdb_id}")
    data = response.json()
    matches = [ item for item in data["objects"] if item["pdbid"].lower() == pdb_id.lower() ]
    assert len(matches) > 0,  f"No OPM record found for {pdb_id}"
    assert len(matches) == 1, f"Found several OPM records for {pdb_id}"
    return matches[0]["id"]

def opm_data(pdb_id):
    # Get the opm_id
    opm_id = get_opm_id(pdb_id)
    # Perform the request with the found opm_id
    response = requests.get(f"{host}/{opm_id}")
    assert response.status_code == 200, f"Received code {response.status_code} for {host}/{opm_id}"
    return response.text

data = opm_data(pdb_id)

with open(output_file, 'w') as file:
    file.write(data)
