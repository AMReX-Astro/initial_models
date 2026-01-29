
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

c = Composition(large_network)
for row in bulk_data:
    c.set_array(list(row)[start:end+1])
    new_c = c.bin_as(aprox19_nuclei)
    print(new_c.get_nuclei())
    print(new_c.get_sum_X())
    print(new_c.get_molar())
    break

# OUTPUT TO FILE
new_file_name = "re-"+file_name

with open(new_file_name, 'w') as f:
    def writeline(array):
        f.write(''.join([str(x).rjust(40) for x in array])+'\n')
    
    writeline(range(1,len(header_data)+1))
    writeline(header_data.keys())
    writeline(header_data.values())
    f.write('\n')

    writeline(range(1,len(bulk_names)+1))
    writeline(bulk_names)
    for i in range(bulk_data.size):
        writeline(bulk_data[i])
    f.write('\n')

