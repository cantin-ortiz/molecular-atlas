# A Molecular Atlas of the Adult Mouse Brain
## Disclaimer
This code is shared for reproducibility purposes. Despite our efforts, there might still be some hardcoded paths and poorly commented parts left. We cannot guarantee this will run fine in your system. The sessionInfo() output is available in the bin directory and will provide guidance about the conditions in which scripts were originally run.

## 01-Registration
Data is available at [...]. Simply cd to the data directory and run the script. It will generate a meta table and a raw count expression table (this will take several hours).

## Figures
Figures were generated using many intermediary data files. These files are available for download at [...]. In order to render the figures, make sure to source the script 'includes.R' from the bin directory. It will load many functions required for plotting. You must update the two variables defined in the two first lines of includes.R : path.bin should be a string with the absolute path to the bin directory, path.matrices a string with the absolute path to the previously downloaded folder containing all the figures data.
