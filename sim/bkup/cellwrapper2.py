# ~ from neuron import h, gui
# ~ import neuron
from neuron import h
import os

def loadCell(cellName, cellTemplateName):
    os.chdir('cell_data/'+cellName+'/')
    h.load_file("stdrun.hoc")
    h.load_file("import3d.hoc")
    print('Loading constants')
    h.load_file('constants.hoc')
    h.load_file("morphology.hoc")
    h.load_file("biophysics.hoc")
    h.load_file("template.hoc")
    # Instantiate the cell from the template
    add_synapses=False
    print ("Loading cell",cellTemplateName)
    cell = getattr(h, cellTemplateName)(1 if add_synapses else 0)    
    print (cell)
    return cell
