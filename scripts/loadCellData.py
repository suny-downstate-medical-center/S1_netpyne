import random
import csv
import datetime
import json
import os, sys, traceback
import pandas as pd
import numpy as np
import zipfile
from shutil import copyfile, copytree

# the file you need to import from netpyne is template.hoc
# and just need these 3 other files (+ the subfolders):
# load_file(“morphology.hoc”)
# load_file(“biophysics.hoc”)
# load_file(“synapses/synapses.hoc”)

copyFileList = ['biophysics.hoc', 'template.hoc', 'morphology.hoc']
dirList = ['morphology', 'synapses']

def loadData():
    try:
        # fileList = os.listdir("/Users/mitras02/self/phd/thesis/S1/hoc_combos_syn.1_0_10.allzips")
        # for index, fileName in enumerate(fileList):
        #     print("/Users/mitras02/self/phd/thesis/S1/hoc_combos_syn.1_0_10.allzips/" + fileName)
        #     with zipfile.ZipFile("/Users/mitras02/self/phd/thesis/S1/hoc_combos_syn.1_0_10.allzips/" + fileName,"r") as zip_ref:
        #         zip_ref.extractall("/Users/mitras02/self/phd/thesis/S1/tempunzip")
            # if index > 2:
            #     break
        fileList = os.listdir("/Users/mitras02/self/phd/thesis/S1/tempunzip")
        for index, fileName in enumerate(fileList):
            print("/Users/mitras02/self/phd/thesis/S1/tempunzip/" + fileName)
            if os.path.isdir("/Users/mitras02/self/phd/thesis/S1/tempunzip/" + fileName):
                if not os.path.isdir("/Users/mitras02/self/phd/thesis/S1/cell_data/" + fileName):
                    os.mkdir("/Users/mitras02/self/phd/thesis/S1/cell_data/" + fileName)
                [ copyfile("/Users/mitras02/self/phd/thesis/S1/tempunzip/" + fileName + "/" + x, "/Users/mitras02/self/phd/thesis/S1/cell_data/" + fileName + "/" + x) for x in copyFileList]
                [ copytree("/Users/mitras02/self/phd/thesis/S1/tempunzip/" + fileName + "/" + x, "/Users/mitras02/self/phd/thesis/S1/cell_data/" + fileName + "/" + x) for x in dirList]
                # if index > 2:
                #     break
    except:
        traceback.print_exc(file=sys.stdout)
    return

loadData()
