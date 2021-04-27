#!/usr/bin/python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

df = pd.read_csv("./queue-1.tr", delim_whitespace=True)

warmup = 35
df = df.iloc[35:]

print(df["0"].mean())
