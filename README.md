# Codec
This project contains an WAV to AAC encoder and an AAC to WAV decoder.
It is implemented in MATLAB and the project is divided into many components of the encoding an the decoding process.

## decodeHuff.m
This file contains a Huffman Decoder function.

```decCoeffs = decodeHuff(huffSec, huffCodebook, huffLUT)```
Performs huffman decoding.

#### Input
- huffSec: a string of '1' and '0' corresponding to the Huffman encoded stream
- huffCodebook: the index (0 to 12) of the codebook used, as outputted by encodeHuff
- huffLUT: the Huffman look-up tables to be loaded using loadLUT.m 

#### Output
- decCoeffs: the decoded quantised (integer) values


## encodeHuff.m
This file contains a Huffman Encoder set of functions.
```[huffSec, huffCodebook] = encodeHuff(coeffSec, huffLUT, forcedCodebook)```
Performs huffman coding, the codebook forcedCodebook to be used.
#### Input
- coeffSec: quantised (integer) values of a section 
- huffLUT: the Huffman look-up tables to be loaded using loadLUT.m 
- forcedCodebook: the codebook to be used

#### Output
- huffSec: string of '1' and '0' corresponding to the Huffman encoded stream
- huffCodebook: the number of the Huffman codebook used
---

```[huffSec,  huffCodebook]=encodeHuff1(coeffSec,huffLUT)```
Performs huffman coding.
#### Input
- coeffSec: quantised (integer) values of a section 
- huffLUT: the Huffman look-up tables to be loaded using loadLUT.m 

#### Output
- huffSec: string of '1' and '0' corresponding to the Huffman encoded stream
- huffCodebook: the number of the Huffman codebook used
---

```[huffSec]=huffLUTCode1(huffLUT, coeffSec)```
Performs the actual Huffman coding. It is not called by the user.

```[huffSec]=huffLUTCode0()```
Returns a blank sequence. It is not called by the user.

```[huffSec]=huffLUTCodeESC(huffLUT, coeffSec)```
Performs the actual Huffman coding with the escape sequence. It is not called by the user.


