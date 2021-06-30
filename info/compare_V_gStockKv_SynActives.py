import numpy as np
import os
if 'DISPLAY' in os.environ:
    del os.environ['DISPLAY']
import sys
from matplotlib import pyplot as plt

rootFolder = '/home/fernando/S1detailed/'
os.chdir(rootFolder)

folder = os.listdir('cell_data/')


def compareTraces(cellnumber):

    os.chdir(outFolder)
    import neuron
    neuron.h.load_file("./init.hoc");
    neuron.h.create_cell(1); # argument 1 stands for 'load synapses'

    cell = neuron.h.cell
    soma = cell.soma[0]

    stimulus = neuron.h.IClamp(0.5, sec=soma)
    stimulus2 = neuron.h.IClamp(0.5, sec=soma)

    stimulus.dur = durationstim # ms
    stimulus.delay = delaystim  # ms     
    stimulus2.dur = timesimulation # ms
    stimulus2.delay = 0  # ms    

    step_number=1
    stimulus.amp = step1_current
    stimulus2.amp = holding_current

    recordings = {}

    recordings['time'] = neuron.h.Vector()
    recordings['soma(0.5)'] = neuron.h.Vector()

    recordings['time'].record(neuron.h._ref_t, 0.1)
    recordings['soma(0.5)'].record(cell.soma[0](0.5)._ref_v, 0.1)

    neuron.h.dt = 0.05
    neuron.h.cvode_active(0)
    neuron.h.tstop = timesimulation # ms

    with open('synapses/mtype_map.tsv') as mtype_map_file:
        mtype_map_content = mtype_map_file.read()

    mtype_map = {}
    for line in mtype_map_content.split('\n')[:-1]:
        n, mtype = line.split()
        mtype_map[mtype] = int(n)

    print (mtype_map)

    def init_synapses(enabled_mtypes=[]):
        """Enable all the synapses that are projected onto this cell from mtype listed in enabled_mtypes."""
        enabled_mtype_ints = [mtype_map[mtype] for mtype in enabled_mtypes]

        for i in range(0, int(cell.synapses.n_of_mtypes)): # Loop over all the m-type
            if i in enabled_mtype_ints: # Enable synapses
                #  The [were_]active_pre_mtypes is a NEURON vector 
                # (it uses the .x syntax to access the elements)
                # When the value in the vector is 1 all the presynaptic neurons
                # of a particular m-types are active (and inactive when it is 0)
                cell.synapses.were_active_pre_mtypes.x[i]= 0
                cell.synapses.active_pre_mtypes.x[i] = 1        
            else: # Disable synapses
                cell.synapses.were_active_pre_mtypes.x[i]= 1
                cell.synapses.active_pre_mtypes.x[i] = 0

        cell.synapses.update_synapses(neuron.h.synapse_plot); # Update the synapses


    # Enable incoming synapses from all mtypes
    init_synapses(enabled_mtypes=mtype_map.keys())

    exc_cells = ['L23_PC', 'L4_PC', 'L4_SS', 'L4_SP', 
                 'L5_TTPC1', 'L5_TTPC2', 'L5_STPC', 'L5_UTPC',
                 'L6_TPC_L1', 'L6_TPC_L4', 'L6_BPC', 'L6_IPC', 'L6_UTPC']
    for mtype in mtype_map:
        if mtype in exc_cells:
            freq = 10.0 # [Hz]
        else:
            freq = 10.0 # [Hz]
        cell.synapses.pre_mtype_freqs.x[mtype_map[mtype]]=freq

    cell.synapses.update_synapses(neuron.h.synapse_plot);

    recordings = {}

    recordings['time'] = neuron.h.Vector()
    recordings['soma(0.5)'] = neuron.h.Vector()

    recordings['time'].record(neuron.h._ref_t, 0.1)
    recordings['soma(0.5)'].record(cell.soma[0](0.5)._ref_v, 0.1)

    time = neuron.h.Vector()
    voltage = neuron.h.Vector()
    StochKvgk = neuron.h.Vector()
    StochKvN1 = neuron.h.Vector()
    ik = neuron.h.Vector()
    mtim = neuron.h.Vector()

    time.record(neuron.h._ref_t)
    voltage.record(soma(.5)._ref_v);
    StochKvgk.record(soma(.5)._ref_gk_StochKv);
    StochKvN1.record(soma(.5)._ref_N1_StochKv);
    ik.record(soma(.5)._ref_ik);
    mtim.record(soma(.5)._ref_m_Im);

    neuron.h.run()

    time1 = np.array(recordings['time'])
    soma_voltage = np.array(recordings['soma(0.5)'])

    N = cell.soma[0](0.5).N_StochKv

    determnistic = 1

    if determnistic:
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
        for secs in cell.dend:
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
        for secs in cell.axon:
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

    voltage = neuron.h.Vector()
    StochKvn_q = neuron.h.Vector()
    StochKvgion = neuron.h.Vector()

    voltage.record(soma(.5)._ref_v);
    StochKvn_q.record(soma(.5)._ref_n_q_StochKv_deterministic);
    StochKvgion.record(soma(.5)._ref_gion_StochKv_deterministic);

    neuron.h.run()

    figSize = (15,10)
    fig = plt.figure(figsize=figSize)  # Open a new figure
    fontsiz=18

    plt.subplot(3, 1, 1)
    plt.ylabel('Voltage (mV)', fontsize=fontsiz)
    plt.plot(time1,soma_voltage, linewidth=2.0, color='blue', label='BBP')
    plt.plot(time,voltage, linewidth=3.0, color='red', label='BBPdet')
    plt.xlim(0, timesimulation)
    plt.grid(True)
    plt.legend(loc='upper right', bbox_to_anchor=(1.1, 1.0))

    plt.subplot(3, 1, 2)
    plt.ylabel('g_Stoch (mS)', fontsize=fontsiz)
    plt.plot(time,StochKvgk*1e2, linewidth=2.0, color='blue', label='BBP')
    plt.plot(time,StochKvgion*1e6, linewidth=3.0, color='red', label='BBPdet')
    plt.xlim(0, timesimulation)
    plt.grid(True)
    plt.legend(loc='upper right', bbox_to_anchor=(1.1, 1.0))

    plt.subplot(3, 1, 3)
    plt.ylabel('n_q', fontsize=fontsiz)
    plt.plot(time,StochKvN1/N, linewidth=2.0, color='blue', label='BBP')
    plt.plot(time,StochKvn_q, linewidth=3.0, color='red', label='BBPdet')
    plt.xlim(0, timesimulation)
    plt.grid(True)
    plt.legend(loc='upper right', bbox_to_anchor=(1.1, 1.0))
    plt.xlabel('Time (ms)', fontsize=fontsiz)

    plt.ion()
    plt.tight_layout()
    plt.savefig(rootFolder+'Figures-comparation/comparison_InputSynActives_%s.png' % folder[cellnumber])
    
if __name__ == '__main__':
    if len(sys.argv) == 1:
        print ("Script need a cell number between 0 and 1034")
    elif int(sys.argv[1])>=0 and int(sys.argv[1])<=1034:
        print ("Comparing BBP and Netpyne Traces of:")
        cellnumber = int(sys.argv[1])
        cellName = folder[cellnumber]
        outFolder = rootFolder+'cell_data/'+folder[cellnumber]
        
        print ("CellNumber = %d" % cellnumber)
        print ("CellName = %s" % cellName)

        with open(outFolder + '/current_amps.dat') as current_file:
            current_content = current_file.read()

        holding_current, step1_current, step2_current, step3_current = [float(x) for x in current_content.split()]
        print ('load step1_current from current_amps.dat = %s' % step1_current)
        
        # ~ holding_current = 0
        step1_current = step3_current
                
        durationstim = 400
        delaystim = 200
        timesimulation = 700
        
        compareTraces(cellnumber)             
        
    else:
        raise Exception('Script need a cell number between 0 and 1034')
