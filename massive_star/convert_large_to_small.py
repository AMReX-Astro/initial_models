
import numpy as np
import pynucastro
from pynucastro.networks.rate_collection import Composition

file_name = "20m_pre_cc_1335s_206_isos.data"

# I have a branch in mesa_reader that could be used to make this code shorter
# https://github.com/melilly/py_mesa_reader/tree/write-file

bulk_data = np.genfromtxt(
    file_name,
    skip_header=5,
    names=True,
    ndmin=1,
    dtype=None,
)
bulk_names = bulk_data.dtype.names
with open(file_name) as f:
    for i, line in enumerate(f):
        if i == 1:
            header_names = line.split()
        elif i == 2:
            header_data = [eval(datum) for datum in line.split()]
        elif i > 2:
            break
header_data = dict(zip(header_names, header_data))


# SHRINK WITH PYNUCASTRO
try:
    start, end = bulk_names.index('neut'), bulk_names.index('zn66')
    large_network = bulk_names[start:end+1]
except ValueError as e:
    print("The given start or end species is not listed in the MESA file.")
    raise e

aprox19_nuclei = ['h1', 'he3', 'he4', 'c12', 'n14', 'o16', 'ne20', 'mg24', 
                  'si28', 's32', 'ar36', 'ca40', 'ti44', 'cr48', 'fe52', 
                  'fe54', 'ni56', 'n']
# this excludes p_nse. We'll deal with that split later

small_network = aprox19_nuclei

network_col_names = []
network_col_bulk = []
c = Composition(large_network)
for row in bulk_data:
    c.set_array(list(row)[start:end+1])
    new_c = c.bin_as(aprox19_nuclei)
    dictionary = new_c.data

    if not network_col_names:
        network_col_names = dictionary.keys()

    network_col_bulk.append(list(dictionary.values()))
    
network_col_bulk = np.array(network_col_bulk)

# OUTPUT TO FILE
new_file_name = "re-"+file_name

with open(new_file_name, 'w') as f:
    def writeline(array):
        f.write(''.join([str(x).rjust(40) for x in array])+'\n')
    def replace_cols(data, replacement, start, end):
        return list(data)[:start]+list(replacement)+list(data)[end+1:]
    
    writeline(range(1,len(header_data)+1))
    writeline(header_data.keys())
    writeline(header_data.values())
    f.write('\n')

    col_names = replace_cols(bulk_names, network_col_names, start, end)
    writeline(range(1,len(col_names)+1))
    writeline(col_names)
    for i in range(bulk_data.size):
        writeline(replace_cols(bulk_data[i], network_col_bulk[i], start, end))
    f.write('\n')

