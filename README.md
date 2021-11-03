Code for the paper:
A Simple Spatial Extension to the Extended Connectivity Interaction Features for Binding Affinity Prediction

There are two folders (a) `data_generation` and (b) `benchmarking`.
With (a), one can generate the datasets used in our experiments from the base mol and mol2 files.
To generate the ligand, ECIF and PDECIF datasets used in our experiments, one should download the `data_generation_input.zip` file from `https://drive.google.com/drive/folders/16Pkvg9gx8BeFkCHnbzfuCRCJvmEU1QKH?usp=sharing` and place the extracted content (proteins and ligands folders) in the `data_generation\input\` folder.
Note: We downloaded the `PDB_Atom_Keys.csv` file used in the original ECIF paper from:
`https://github.com/DIFACQUIM/ECIF`.  

With (b), one can perform the benchmarking we did.
To run the exact benchmarks which generated the values in the paper, download the `benchmarking_input.zip` folder from `https://drive.google.com/drive/folders/16Pkvg9gx8BeFkCHnbzfuCRCJvmEU1QKH?usp=sharing`, then extract and place the folders (datasets and responses) in `\benchmarking\input\`.
Executing the file in `benchmarking\code\` will run all the benchmarks and output the results to `benchmarking\output`. 