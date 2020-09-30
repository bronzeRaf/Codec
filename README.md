# Codec
This project contains an WAV to AAC encoder and an AAC to WAV decoder.
It is implemented in MATLAB and the project is divided into many components of the encoding an the decoding process.

## decodeHuff.m
This file contains a Huffman Decoder function.

## encodeHuff.m
This file contains a Huffman Encoder set of functions.

## loadLUT.m
This file contains creates the Huffman look-up tables.

## huffCodebookSF.mat
This file is MATLAB workspace, used to load the Huffman codebook.

## huffCodebook.mat
This file is MATLAB workspace, used to load the Huffman codebook.

## SSC.m
This file contains function that decide the type of a sound frame.

## filterbank.m
This file contains a set of functions that produce the MDCT after windowing of a sound frame. The windowing could be either Kaizer windows either Sinusoid windows.

## ifilterbank.m
This file contains a set of functions that produce the sound frame of a window. It contains the inverse proceses of the filterbank.m.

## TNS.m
This file contains a set of functions that filters a frame to apply temporal noise shaping. It transforms a MDCT to a new group of coeffs without periodicity.

## iTNS.m
This file contains a set of functions that transforms a temporal noise shaping set of coeffs to a MDCT frame. It contains the inverse proceses of the TNS.m.

## AACquantizer.m
This file contains a function that quantizes a MDCT frame.

## iAACquantizer.m
This file contains a function that dequantizes a frame to restore the normal MDCT frame. It contains the inverse proceses of the AACquantizer.m.

## AACoder3.m
This file contains a set of functions that encodes a WAV file to AAC.

## iAACoder3.m
This file contains a set of functions that decodes an AAC file to a WAV file. It contains the inverse proceses of the AACoder3.m.
