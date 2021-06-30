import random
import csv
import datetime
import json
import os, sys, traceback
import pandas as pd
import numpy as np

def loadData():
    try:
        anatomy_data = json.loads(open("pathways_anatomy_s1.json").read())
        print(anatomy_data.keys() )

        columnNames = ["Connection", "From Cell", "From Layer", "From Type", "To Cell", "To Layer", "To Type"] + list(anatomy_data[list(anatomy_data.keys() )[0]].keys())
        df = pd.DataFrame(columns=columnNames, data = [[k,  k.split(":")[0],  k.split(":")[0].split("_")[0], k.split(":")[0].split("_")[1], k.split(":")[1], k.split(":")[1].split("_")[0], k.split(":")[1].split("_")[1]] + list(v.values()) for k, v in anatomy_data.items()    ])
        df = df.sort_values(by=['Connection'])
        df.to_csv("anatomy_out.csv", index=None)
        # print (df.head())
        physiology_data = json.loads(open("pathways_physiology_s1.json").read())
        physColumnNames = ["Connection", "From Cell", "From Layer", "From Type", "To Cell", "To Layer", "To Type"] + list(physiology_data[list(physiology_data.keys() )[0]].keys())
        df2 = pd.DataFrame(columns=physColumnNames, data = [[k,  k.split(":")[0],  k.split(":")[0].split("_")[0], k.split(":")[0].split("_")[1], k.split(":")[1], k.split(":")[1].split("_")[0], k.split(":")[1].split("_")[1]] + list(v.values()) for k, v in physiology_data.items()    ])
        df2 = df2.sort_values(by=['Connection'])
        df2.to_csv("physiology_out.csv", index=None)
        # print (df2.head())
        df_merged = df.merge(df2, how="outer", on = ["Connection", "From Cell", "From Layer", "From Type", "To Cell", "To Layer", "To Type"])
        print(df_merged.head())
        df_merged.to_csv("anatomy_physiology_s1.csv", index=None)
        df_merged.to_hdf('anatomy_physiology_s1.h5', key='df_merged', mode='w')
    except:
        traceback.print_exc(file=sys.stdout)
    return

loadData()
