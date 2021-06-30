This file holds the connectivity of an instance of a modeled microcircuit as described in (Markram et al., 2015; Cell). It uses the hdf5 file format.

THIS MICROCONNECTOME IS A STOCHASTIC INSTANCE OF A MICROCIRCUIT BASED ON AVERAGED BIOLOGICAL MEASUREMENTS.

Inside the file, two main data sources can be found, 'populations' and 'connectivity'. 'populations' holds summary information for each individual neuron in the microcircuit, like its position in 3d space, in- and out-degree of connectivity, etc. The information is split into 55 groups, one for each morphological type used in the model (for example 'L4_PC' = layer 4 pyramidal cell). For more information about the morphological types refer to (Markram et al., 2015; Cell). 

These individual groups are named for the morphological type they concern; they are listed in the overview below as $population_name1 to $population_name55.

NOTE: The in- and out-degree as well as the number of afferent and efferent synapses listed in the 'populations' data includes some connections from- or to neurons outside the modeled volume. In order to reduce the edge effect, i.e. the effect that neurons at the boundary of the modeled volume have fewer connections, a modeled microcircuit is surrounded byt 7 further microcircuits and connections from or to the surrounding microcircuits are also taken into account in this count.


The 'connectivity' group contains the actual binary connection matrices. They only include connections within the modeled volume. To make them easier to handle they are also split into submatrices labeled by the name of first the pre-, then the post-synaptic morphological type. 

For example, /connectivity/L4_MC/L4_PC/cMat contains the connections FROM L4_MCs TO L4_PCs:

In [0]: h5['connectivity/L4_MC/L4_PC/cMat'].shape
Out[0]: (118, 2674)

The shape is 118 by 2674 because there are 118 L4_MCs in this instance of the model and 2674 L4_PCs, i.e. the presynaptic neurons are along the first axis of the matrix, the postsynaptic neurons along the second axis.

In [1]: numpy.sum(h5['connectivity/L4_MC/L4_PC/cMat'])
Out[1]: 4381

There are 4381 connections from L4_MC to L4_PC neurons


Overview

cons_locs_pathways.h5
    /populations
        /$population_name1
        /...
        /$population_name55
            /locations : Nx3 dataset of neuron locations in 3d space (in um; one per neuron in the indicated population)
            /nCellAff : Nx1 dataset of number of afferent connections (one per neuron in the indicated population; NOTE: INCLUDES SOME CONNECTIONS FROM NEURONS OUTSIDE THE VOLUME, SEE ABOVE!)
            /nCellEff : Nx1 dataset of number of efferent connections (one per neuron in the indicated population; NOTE: INCLUDES SOME CONNECTIONS TO NEURONS OUTSIDE THE VOLUME, SEE ABOVE!)
            /nSynAff : Nx1 dataset of number of afferent synapses (multiple synapses per connection possible)
            /nSynEff : Nx1 dataset of number of efferent synapses
    /connectivity
        /$population_name1
        /...
        /$population_name55
            /$population_nameX
            /...
            /$population_nameY
                /cMat : NxM dataset holding the connection matrix from neurons of one population to neurons of another population. The population name appearing hierarchically higher in the file structure is the presynaptic type, the other the postsynaptic type.
