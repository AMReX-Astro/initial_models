import numpy as np
import argparse

'''
This script turns a file in MESA's log file format:
                           1                            2
                model_number                    num_zones
                       13700                         4055

                           1                            2
                          dm                            q
     1.6956238372097480E+022      1.0000000000000000E+000
     1.6956238372097480E+022      9.9999999999949996E-001

into Castro's new format:
# radius density temperature
        0.1           3         4.0
        0.2           4         3.9
        0.3           4         3.7
'''


def main(file_name, new_file_name):
    # read in MESA file
    bulk_data = np.genfromtxt(
        file_name,
        skip_header=5,
        names=True,
        ndmin=1,
        dtype=None,
    )
    bulk_names = bulk_data.dtype.names

    # identify needed columns and their aliases
    variable_names = {
        "radius": "radius_cm",
        "density": "density",
        "temperature": "temperature",
        "ye": "ye",
        "velx": "velocity",
        "pressure": "pressure",
        "entr": "entropy",
        "abar": "abar",
        "enuc": "eps_nuc"
    } # castro column name: mesa column name
    
    # add in the composition columns:
    # found in the range [neutron, opacity)
    neutron_name = {'n', 'neut', 'neutron'}.intersection(bulk_names)
    if len(neutron_name) < 1:
        raise ValueError("No 'n', 'neut', or 'neutron' in MESA file column names.")
    start = bulk_names.index(neutron_name.pop())
    end = bulk_names.index('opacity')

    # set up new columns
    col_names = list(variable_names.keys()) + list(bulk_names[start:end])
    associated_indices = [bulk_names.index(v) for v in variable_names.values()] + list(range(start, end))

    # write new columns to new file
    with open(new_file_name, 'w') as f:
        f.write(f"# {' '.join(col_names)}\n")
        for row in bulk_data:
            row = list(row)
            display_row = [str(row[i]).rjust(30) for i in associated_indices]
            f.write(f'{"".join(display_row)}\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reformat a 1D model file from MESA's output format to Castro's preferred input format.")
    parser.add_argument('file_name')
    parser.add_argument('new_file_name', nargs='?')
    args = parser.parse_args()
    if args.new_file_name is None:
        file_stem, file_ending = args.file_name.rsplit('.',1)
        new_file_name = f"{file_stem}.castro.{file_ending}"
    else:
        new_file_name = args.new_file_name
    main(args.file_name, new_file_name)

    print('Saved new Castro-readable file in', new_file_name)
