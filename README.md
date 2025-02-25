[![DOI](https://zenodo.org/badge/902773932.svg)](https://doi.org/10.5281/zenodo.14926086)
# Quantification of LAMP1 nearest distance to TDP43 clusters

This repository contains an ImageJ macro script for analyzing cellular microscopy data, specifically designed to process and analyze images containing nuclei, spots, and cluster markers across multiple cells.

## Requirements

- ImageJ version 1.54f
- Java 1.8.0_322 (64-bit)
- MorphoLibJ plugin version 1.6.3
- Clij2 plugin plugin version 2.5.3.5
- Bio-Formats plugin version 7.3.1

## Features

The script performs the following analyses:
- Cell segmentation from user defined ROIs
- Nuclei segmentation
- LAMP1 segmentation
- TDP43 segmentation
- Distance measurements between LAMP1 and TDP43 clusters
- Generation of histograms for spatial distributions

## Input Requirements

The script expects:
1. A directory containing multichannel stacks (supported formats: .nd2 or .tif)
2. A directory containing ROI files (.zip) corresponding to the images
3. Multi-channel images with at least:
   - Nuclei channel
   - TDP43 channel
   - LAMP1 channel

## Parameters

### User-Defined Parameters (GUI)
- File directory: Location of microscopy images
- ROI directory: Location of ROI files
- Batch process: Yes/No option
- Test file number: For testing (0-100)
- Test ROI number: For testing (0-10)
- File suffix: .nd2 or .tif

### Script Parameters
- Pixel size: 0.0693948 μm
- Z-step: 0.5 μm
- Nuclei filter size: 1000 μm³
- Cluster threshold: 7000
- Bin sizes for various measurements:
  - Separation: 0.5
  - Distance: 100
  - Scaled distance: 0.1

## Output

The script creates two main output directories:
1. `/Analysis`: Contains CSV files with measurements for each ROI:
   - Local separation measurements
   - Distance to cluster measurements
   - Scaled distance measurements
2. `/Binaries`: Contains processed binary images

## Usage

1. Install all required plugins in ImageJ
2. Open the script in ImageJ
3. Run the script and input the required parameters in the GUI
4. The script will process all images in the specified directory and generate corresponding analysis files

## Notes

- GPU acceleration is used for spot detection (requires NVIDIA GPU)
- All measurements are calibrated using the specified pixel size and Z-step
- The script includes extensive error handling and intermediate cleanup to manage memory usage

## Output File Format

The CSV output files contain the following columns:
- Values[um]: Measurement values in micrometers
- Counts[px]: Raw pixel counts
- Counts_fraction: Normalized counts
- Cell_id: Identifier for each analyzed cell
- Aggregate_ratio: Ratio of cluster to cell volume
- Filename: Source image filename
