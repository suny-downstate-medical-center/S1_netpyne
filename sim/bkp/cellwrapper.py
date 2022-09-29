from neuron import h
import os
import sys

def loadCell(cellName, cellTemplateName):
    origDir = os.getcwd()
    os.chdir('cell_data/'+cellName+'/')
    h.load_file("stdrun.hoc")
    h.load_file('import3d.hoc')
    try:
        h.xopen("morphology.hoc")
    except:
        pass
    try:
        h.xopen("biophysics.hoc")
    except:
        pass
    try:
        h.xopen("synapses/synapses.hoc")
    except:
        pass
    h.xopen('template.hoc')
    cell = getattr(h, cellTemplateName)(0)
    i=0
    for secs in cell.somatic:
        sec = cell.soma[i]
        listmech = list(cell.soma[i](0.5))      
        for mech in listmech:
            if str(mech) == 'StochKv':
                sec.insert('StochKv_deterministic')
                cell.soma[i].gmax_StochKv_deterministic = 1e-4 * cell.soma[i].gkbar_StochKv
                # ~ print (sec, mech, i, cell.soma[i].gmax_StochKv_deterministic)
                sec.uninsert('StochKv')
        i=i+1

    i=0
    for secs in cell.basal:
        sec = cell.dend[i]
        listmech = list(cell.dend[i](0.5))      
        for mech in listmech:
            if str(mech) == 'StochKv':
                sec.insert('StochKv_deterministic')
                cell.dend[i].gmax_StochKv_deterministic = 1e-4 * cell.dend[i].gkbar_StochKv
                # ~ print (sec, mech, i, cell.dend[i].gmax_StochKv_deterministic)
                sec.uninsert('StochKv')
        i=i+1

    i=0
    for secs in cell.apical:
        sec = cell.apic[i]
        listmech = list(cell.apic[i](0.5))      
        for mech in listmech:
            if str(mech) == 'StochKv':
                sec.insert('StochKv_deterministic')
                cell.apic[i].gmax_StochKv_deterministic = 1e-4 * cell.apic[i].gkbar_StochKv
                # ~ print (sec, mech, i, cell.apic[i].gmax_StochKv_deterministic)
                sec.uninsert('StochKv')
        i=i+1

    i=0
    for secs in cell.axonal:
        sec = cell.axon[i]
        listmech = list(cell.axon[i](0.5))      
        for mech in listmech:
            if str(mech) == 'StochKv':
                sec.insert('StochKv_deterministic')
                cell.axon[i].gmax_StochKv_deterministic = 1e-4 * cell.axon[i].gkbar_StochKv
                # ~ print (sec, mech, i, cell.axon[i].gmax_StochKv_deterministic)
                sec.uninsert('StochKv')
        i=i+1     
    
    print (cell)
    os.chdir(origDir)
    return cell
