
import numpy as np

file_name = "20m_pre_cc_1335s_206_isos.data"

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

