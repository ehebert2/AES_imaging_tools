# AES_imaging_tools
Tools for processing imaging data produced using adaptive excitation. This includes three installable MATLAB apps written in app designer. A description of the binary file data output by the motion correction app and a matlab class for reading files are included in "aes_file_format.txt" and "AESFile.m" respectively alongside the matlab app packages in the /installation folder.

## Motion Correction App
This is the primary tool present in this repository. When given masks depicting the regions of the image which are illuminated, the program can calculate rigid x-y displacement for each frame, apply the inverse operation to the video, and write the displacement information and corrected video to file. Additionally, given masks segmenting regions of interested (i.e. neuron locations), it can save the pixel values after motion correction to file alongside average per frame background values (electrical background, not flourescence background) and determine when the segmented ROI drifts outside of the exposure region. A more complete description of the tool, necessary inputs, and potential outputs is included in the "Motion Correction App Manual". The code relies on the built in LibTiff library for fast reading and writing. Due to older LibTiff realeses in some Matlab versions, installations of the app on Matlab 2023a and earlier are unable to read Tiff stacks with more than 65000 frames.

## Mask App
A simple tool for generating binary masks over an image and saving them as uncompressed '.tif' files. It is helpful for producing the necessary inputs for the motion correction tool, but the masks can be generated in other ways.

## Trace Viewer App
A tool for getting a quick look at the output from the motion correction app before deciding whether a more thorough analysis is worthwhile. It is not intended for final processing of the data. A more complete description of this tool is contained in the "Trace Viewer App Manual" word document.

# Installation
The apps can either be run from the source code in the app designer environment or installed as MATLAB apps. Additionally, the Trace Viewer App and Motion Correction App can be installed as standalone applications using the MATLAB 2023b runtime, allowing them to run without opening MATLAB. The MATLAB app install files are located in installation\matlab_app while the standalone installation excecutables are located in installation\standalone_app. \for_redistribution subfolders contain the excecutable for installing the full application with the runtime while the \for_redistribution_files_only subfolders contain an application file that requires the MATLAB runtime to be already installed.

The code was tested on Windows 10 and Windows 11 machines running MATLAB 2022a and 2023b. Additionally, the code requires the Signal Processing and Image Processing Toolboxes. Due to the use of the LibTiff library integrated into MATLAB for accessing Tiff stacks, versions of MATLAB prior to 9.14 (2023a) are unable to read Tiff stacks with more than 2^16 frames. The code will truncate all data to 64000 frames if using prior releases.
