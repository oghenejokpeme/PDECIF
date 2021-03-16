Code for the paper:
A Simple Spatial Extension to the Extended Connectivity Interaction Features for Binding Affinity Prediction

There are two folders (a) `data_generation` and (b) `benchmarking`.
With (a), one can generate the datasets used in our experiments from the base mol and mol2 files.
To generate the ECIF and PDECIF datasets used in our experiments, one should download the zipped CASF benchmark dataset from `https://drive.google.com/file/d/1B7bUGcVTU-wBCSMs6TdPfF9o1WC_nhia/view?usp=sharing` and extract to the `data_generation\input\` folder. There should be three folders for each year discussed in the paper.
Running the files in `data_generation\src\` should generate the datasets and place them in the `data_generation\output\` folder.
Note: We downloaded the `PDB_Atom_Keys.csv` file used in the original ECIF paper from:
`https://github.com/DIFACQUIM/ECIF`.  

With (b), one can perform the benchmarking we did.
To run the exact benchmarks which generated the values in the paper, download the input folder from `https://drive.google.com/file/d/1Mwwq78UdvclHAQ4mBdjqSRtgAiO8JhvG/view?usp=sharing`, then extract and place the folder in `\benchmarking\` as `input`.
Executing the file in `benchmarking\code\` will run all the benchmarks and output the results to `benchmarking\output`. 