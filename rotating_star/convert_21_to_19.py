# take a raw inputs file for aprox21 and convert to the aprox19
# nuclei.  This means getting rid of Cr56 and Fe56 (by lumping them
# into Ni56)

import numpy as np

file = "15m_500_sec.txt"

SKIP_ROWS = 31

data = np.loadtxt(file, skiprows=SKIP_ROWS)

print(data.shape)

# now manually read to get the variable names

# the first column is position

names = ["r"]

with open(file) as f:
    for n, line in enumerate(f):
        if line.startswith("#"):
            continue

        if line.startswith("number"):
            continue

        if n >= SKIP_ROWS:
            break

        names.append(line.strip())


# create a new data

data_new = np.zeros((data.shape[0], data.shape[1]-2))

names_new = list(names)
names_new.remove("Cr56")
names_new.remove("Fe56")

extra = np.zeros((data.shape[0]))

# loop over the columns in the original file

for n, name in enumerate(names):

    # is this name present in the new file? if so, just copy,
    # otherwise we will add the data to the "extra"

    found = True
    try:
        nnew = names_new.index(name)
    except ValueError:
        print(f"not found {name}")
        found = False

    if found:
        data_new[:,nnew] = data[:,n]
    else:
        extra[:] += data[:,n]

# add the extra mass fractions to Ni56

data_new[:,names_new.index("Ni56")] += extra[:]

print(extra.min())
print(extra.max())


# output the new model

base = file[:file.find(".txt")]

with open(f"{base}.aprox19.dat", "w") as f:
    f.write(f"# conversion of {file} to aprox19 nuclei via convert_21_to_19.py\n")
    f.write(f"number of variables = {len(names_new)-1}\n")
    for name in names_new:
        if name == "r":
            continue
        f.write(f"{name}\n")
    for irow in range(data.shape[0]):
        l = [f"{q:30.20g}" for q in data[irow, :]]
        f.write(" ".join(l) + "\n")



