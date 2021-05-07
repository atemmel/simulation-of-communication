#!/usr/bin/python

import matplotlib.pyplot as plt

our_data = []
their_data = []

our_file = open("./our_results.txt", "r");
for line in our_file:
    vals = line.split(' ')
    for val in vals:
        if len(val) > 0:
            our_data.append(float(val))
our_file.close()

their_file = open("./their_results.txt", "r");
for line in their_file:
    vals = line.split(' ')
    for val in vals:
        if len(val) > 0:
            their_data.append(float(val))

n_bins = 20

# Generate a normal distribution, center at x=0 and y=5
fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=True)

# We can set the number of bins with the `bins` kwarg
axs[0].hist(our_data, bins=n_bins)
axs[1].hist(their_data, bins=n_bins)

plt.show()
