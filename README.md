# EPIC450k_repro
### Reproducible analysis of the manuscript _Preprocessing, normalization and integration of the Illumina HumanMethylationEPIC array_
--------

**Authors**: Jean-Philippe Fortin, Tim Triche and Kasper Daniel Hansen


- Step 1: download raw data using the script `scripts/download.sh`. The raw IDAT files will be saved into the `data_raw` folder.
- Step 2: process the data using the script `processing/preprocess.R`. The processed data will be saved into the `data_processed` folder.
- Step 3: normalize the data with the 3 normalization scripts contained in the `processing` folder.
- Step 4: perform the analyses using the scripts contained in the `analysis` folder. The figures will be saved in the folder `figures`.

