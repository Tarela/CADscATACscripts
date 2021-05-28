#!/bin/bash

import sys
import pandas as pd

peak_file = sys.argv[1]
outfile = sys.argv[2]

df = pd.read_csv(peak_file, header=None, sep='\t')
df = df[[0, 1, 2]]

df.to_csv(outfile, header=False, index=False, sep='\t')

