# GrowClust3D
Julia implementation of the GrowClust program for relative relocation of earthquake hypocenters based on waveform cross-correlation data. This is a testing version before later public release as a registered package. Apologies for the sparse readme; more to come. The test version of the package can be installed using the Julia Pkg manager:

` pkg> add https://github.com/dttrugman/GrowClust3D`

The examples/ directory has three different julia scripts, one to run a serial version of the program without uncertainty quantification, and two parallel options with 100 bootstrap resamples: multithreaded and multiprocessing. To run these codes, navigate to this directory and run:

`julia run_growclust-serial.jl example.serial.inp`

or 

`julia -t8 run_growclust-multithread.jl example.multithread.inp`

or

`julia -p8 run_growclust-multiprocess.jl example.multiprocess.inp`

or


In the second example, I specified the usage of 8 threads, while the third example uses 8 additional processors. Generally, if sufficient resources are available, multiprocessing will be faster than multithreading (especially for large datasets) due to memory considerations.

[Note, to download a local copy of this repository, try `git clone https://github.com/dttrugman/GrowClust3D`.]
