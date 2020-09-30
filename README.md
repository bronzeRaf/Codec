# Codec
This project contains an WAV to AAC encoder and an AAC to WAV decoder.
It is implemented in MATLAB and the project is divided into many components of the encoding an the decoding process.

## decodeHuff.m
This file contains a Huffman Decoder function.

    ```decCoeffs = decodeHuff(huffSec, huffCodebook, huffLUT)```
#### Input
- huffSec: a string of '1' and '0' corresponding to the Huffman encoded stream
- huffCodebook: the index (0 to 12) of the codebook used, as outputted by encodeHuff
- huffLUT: the Huffman look-up tables to be loaded using loadLUT.m 

#### Output
- decCoeffs:the decoded quantised (integer) values


## encodeHuff.m
    [huffSec, huffCodebook] = encodeHuff(coeffSec, huffLUT, forcedCodebook)
bbl;ibn

    [huffSec,  huffCodebook]=encodeHuff1(coeffSec,huffLUT)
vuyvuykl

    [huffSec]=huffLUTCode1(huffLUT, coeffSec)
ytfcy

    [huffSec]=huffLUTCode0()
bhuilbi

    [huffSec]=huffLUTCodeESC(huffLUT, coeffSec)
ikjno
