import json
import os, sys
import pandas as pd
import numpy as np

anatomy_data = json.loads(open("pathways_anatomy_factsheets_simplified.json").read())

columnNames = ["Connection", "From Cell", "From Layer", "From Type", "To Cell", "To Layer", "To Type"] + list(anatomy_data[list(anatomy_data.keys() )[0]].keys())
df = pd.DataFrame(columns=columnNames, data = [[k,  k.split(":")[0],  k.split(":")[0].split("_")[0], k.split(":")[0].split("_")[1], k.split(":")[1], k.split(":")[1].split("_")[0], k.split(":")[1].split("_")[1]] + list(v.values()) for k, v in anatomy_data.items()    ])
df = df.sort_values(by=['Connection'])

physiology_data = json.loads(open("pathways_physiology_factsheets_simplified.json").read())
physColumnNames = ["Connection", "From Cell", "From Layer", "From Type", "To Cell", "To Layer", "To Type"] + list(physiology_data[list(physiology_data.keys() )[0]].keys())
df2 = pd.DataFrame(columns=physColumnNames, data = [[k,  k.split(":")[0],  k.split(":")[0].split("_")[0], k.split(":")[0].split("_")[1], k.split(":")[1], k.split(":")[1].split("_")[0], k.split(":")[1].split("_")[1]] + list(v.values()) for k, v in physiology_data.items()    ])
df2 = df2.sort_values(by=['Connection'])

matrixinfo = df.merge(df2, how="outer", on = ["Connection", "From Cell", "From Layer", "From Type", "To Cell", "To Layer", "To Type"])

matrixoptions = []
for name in matrixinfo.keys():
    matrixoptions.append(str(name))

print(matrixoptions)
print(df.shape,df2.shape,matrixinfo.shape)
print(matrixinfo['Connection'][0],matrixinfo['connection_probability'][0])
print(matrixinfo[matrixoptions[0]][0],matrixinfo[matrixoptions[7]][0])
