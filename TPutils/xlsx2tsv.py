#!/usr/bin/env python

import pandas as pd
import sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i", '--input_xls', type=str, help="input xls file")
parser.add_argument("-s", '--sheet', type=str, help="sheet name")


args = parser.parse_args()

data_xlsx = pd.read_excel(args.input_xls, args.sheet ,index_col=None)

data_xlsx.to_csv(sys.stdout, sep='\t', encoding='utf-8',  index=False)
