# AES_imaging_tools
Tools for processing imaging data produced using adaptive excitation. This includes three installable matlab apps. The ui was written in app designer and tested on Matlab release 2022a. A description of the binary file data output by the motion correction app and a matlab class for reading files are included in "aes_file_format.txt" and "AESFile.m" respectively alongside the matlab app packages in the /installation folder.
* __Motion Correction App:__ The primary application use for motion correcting videos with AES patterns and extracting trace data.
* __Mask App:__ Used for producing binary masks when generating AES patterns. This is useful for producing the masks used for segmentation in the Motion Correction App, but not necessary.
* __Trace Viewer App:__ Used for getting a quick overview of output from the Motion Correction App. This is not typically the final step of data analysis; however, it can be useful for understanding the data put out by the Motion Correction App before deciding whether a more thorough analysis will be useful.

## Motion Correction Tool
This is the primary tool present in this repository. When given masks depicting the regions of the image which are illuminated, the program can calculate rigid x-y displacement for each frame, apply the inverse operation to the video, and write the displacement information and corrected video to file. Additionally, given masks segmenting regions of interested (i.e. neuron locations), it can save the pixel values after motion correction to file alongside a average per frame background values (electrical background, not flourescence background) and tell when the segmented ROI drifts outside of the exposure region. A more complete description of the tool, necessary inputs, and required outputs is included in the "Motion Correction App Manual" word document. The code relies on the built in LibTiff library for fast reading and writing. Due to older LibTiff realeses in some Matlab versions, installations of the app on Matlab 2023a and earlier are unable to read Tiff stacks with more than 65000 frames.

## Mask App
A simple tool for generating binary masks over an image and saving them as uncompressed '.tif' files. It is helpful for producing the necessary inputs for the motion correction tool, but the masks can be generated in other ways.

## Trace App
A tool for getting a quick look at the output from the motion correction app. It is not intended for full processing of the data.
