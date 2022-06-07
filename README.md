# Fractal Analysis on Time-Series Data
This project was developed as a part of a master thesis at the Norwegian University of Science and Technology

## Installation

1. Download and install Anaconda (https://www.anaconda.com/).
2. Create a new conda environment: `conda env create -f environment.yml`.
This will create a new environment called master-thesis  with the packages listed in `environment.yml`.
3. Activate the new environment: `conda activate fractal-analysis`.
4. Open the project with preferred IDE. Jupyter lab was is recommended:
```console
foo@bar:~$ jupyter lab .
```

## Examples

The examples folder contains three notebooks:
* divider_method_example: Shows a illustration of the underlying divider method.
* fractal_dimension_example: Shows how the analysis cleans the output and computes the fractal dimensions.
* sliding_window_example: Takes an input signal and runs a full fractal analysis on it using the sliding window method, and plots the results.

**NOTE: Not permitted to publish or share the dataset used, thus, the data folder is empty.**
