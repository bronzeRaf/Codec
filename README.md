# Codec
This project contains a WAV to AAC encoder and an AAC to WAV decoder.
It is implemented in MATLAB and the project is divided into many components of the encoding and the decoding process.

## Run Codec
- Fire up the coding process with

	``` AACSeq3 = AACoder3('Path\name.wav', 'Path2\workspace.mat'); ```

	where **Path** is the Path to the wav file, **name.wav** is the name of the wav file, **Path2** is the output path for saving the workspace and **workspace.mat** the name of the output workspace. **AACSeq3** is the produced AAC sequence in the memory.
- Fire up the decoding process with

	``` iAACSeq3 = iAACoder3(AACSeq3, 'Path\out.wav'); ```

	where **AACSeq3** is the AAC sequence in the memory, **Path** the output path to save the wav file and **out.wav** the name of the output wav file. **iAACSeq3** is the wav file in the memory.

- Fire up the analysis
	``` [SNR, bitrate, compression] = demoAAC2(' Path1/name.wav', 'Path2/out.wav', Path3/workspace.mat' ) ```

	where **Path1** is the Path to the original wav file, **name.wav** is the name of the original wav file, **Path2** is the path of the decoded wav file, **out.wav** is the name of the decoded wav file, **Path3** is the output path for saving the workspace and **workspace.mat** the name of the output workspace. **SNR**, **bitrate** and **compression** are the stats of the encoding/decoding process.



## Implementation details

### decodeHuff.m
This file contains a Huffman Decoder function.

---
### encodeHuff.m
This file contains a Huffman Encoder set of functions.

---
### loadLUT.m
This file creates the Huffman look-up tables.

---
### huffCodebookSF.mat
This file is a MATLAB workspace, used to load the Huffman codebook.

---
### huffCodebook.mat
This file is a MATLAB workspace, used to load the Huffman codebook.

---
### SSC.m
This file contains functions that decide the type of a sound frame.

---
### filterbank.m
This file contains a set of functions that transform a sound frame into the MDCT after windowing. The windows could be either Kaizer windows either Sinusoid windows.

---
### ifilterbank.m
This file contains a set of functions that transform a window into a sound frame of. It contains the inverse proceses of the filterbank.m.

---
### TNS.m
This file contains a set of functions that filter a frame to apply temporal noise shaping. It transforms an MDCT to a new group of coeffs without periodicity.

---
### iTNS.m
This file contains a set of functions that transform a temporal noise shaping set of coeffs to an MDCT frame. It contains the inverse proceses of the TNS.m.

---
### AACquantizer.m
This file contains a function that quantizes an MDCT frame.

---
### iAACquantizer.m
This file contains a function that dequantizes a frame to restore the normal MDCT frame. It contains the inverse proceses of the AACquantizer.m.

---
### AACoder3.m
This file contains a set of functions that encode a WAV file to AAC.

---
### iAACoder3.m
This file contains a set of functions that decodes an AAC file to a WAV file. It contains the inverse proceses of the AACoder3.m.

---
### demoAAC3.m
This file contains a function that gives some stats about the coding/decoding process. It is able to provide the SNR, the bitrate and the compression rate of the process.
