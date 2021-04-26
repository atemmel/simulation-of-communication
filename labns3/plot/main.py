#!/usr/bin/python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

df = pd.read_csv("./sixth.csv")

x = df["Error Rate"]
ys = df.transpose().iloc[1:]
titles = ["Throughput", "Delay", "Packet Loss"]

print(x)
print(ys)
assert(len(titles) == len(ys))

for i in range(len(titles)):
    plt.plot(x, ys.loc[titles[i]])
    plt.xlabel("Error rate")
    plt.ylabel("Average " + titles[i])
    plt.axis
    #plt.legend(loc=4)
    plt.show()

