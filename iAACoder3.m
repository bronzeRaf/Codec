%AA Decoder level 3

%Decodes the given aac frame sequence to create the wav sample vector (one
%for each channel)
function x = iAACoder3(AACSeq3, fNameOut)
    M = length(AACSeq3);                    %number of frames we have
    tempx = zeros(1024*(M+1), 2);            %preallocate and init tempx, x with zero
    x = zeros(1024*M, 2);

    for i = 1:M     %frame loop
        j = (i - 1) * 1024 + 1; %start of frame in the sample vector
        %decode Huffman
        S1 = decodeHuff(AACSeq3(i).chl.stream, AACSeq3(i).chl.codebook, huffLUT);
        S2 = decodeHuff(AACSeq3(i).chr.stream, AACSeq3(i).chr.codebook, huffLUT);
        %dequantize
        frameF1 = iAACquantizer(S1, AACSeq3(i).chl.sfc, AACSeq3(i).chl.G, AACSeq3(i).frameType);
        frameF2 = iAACquantizer(S2, AACSeq3(i).chr.sfc, AACSeq3(i).chr.G, AACSeq3(i).frameType);
        
        if strcmp(AACSeq3(i).frameType, 'ESH')  %EIGHT SHORT frame
            TNScoeffs(1:4,1:8,1) =  AACSeq3(i).chl.TNScoeffs(1:4,1:8);
            TNScoeffs(1:4,1:8,2) =  AACSeq3(i).chr.TNScoeffs(1:4,1:8);
            %iTNS
            frameFout = iTNS([frameF1, frameF2], AACSeq3(i).frameType, TNScoeffs);
            %get back time samples
            frameT8 = iFilterbank(frameFout, AACSeq3(i).frameType, AACSeq3(i).winType);
            %make frameT 2048x2 again (now is 256x8x2)
            frameT = isubframeMaker(frameT8);
        else                                    %other frame types
            TNScoeffs(1:4,1) = AACSeq3(i).chl.TNScoeffs(1:4,1);
            TNScoeffs(1:4,2) = AACSeq3(i).chr.TNScoeffs(1:4,1);
            %iTNS
            frameFout = iTNS([frameF1, frameF2], AACSeq3(i).frameType, TNScoeffs);
            %get back time samples
            frameT = iFilterbank(frameFout, AACSeq3(i).frameType, AACSeq3(i).winType);
        end
        %add half overlapped frames
        tempx(j:(j + 2047),1) = tempx(j:(j + 2047),1) + frameT(1:2048,1);
        tempx(j:(j + 2047),2) = tempx(j:(j + 2047),2) + frameT(1:2048,2);
    end
    
    x(1:end,1) = tempx(1025:end,1);     %ignore first 1024 frames without info
    x(1:end,2) = tempx(1025:end,2);
    
    fs = 48e3;     %fs = 48KHz
    wavwrite(x,fs,fNameOut);

end

%Gets back the frame matrix from the 8 subframes as we need had before 
%subframeMaker.
function frame = isubframeMaker(frameT)
    frame = zeros(2048,2);
    for i = 1:8
        j = (i - 1) * 128 + 448 + 1;    %start of subframes
        frame(j:(j + 255),1) = frame(j:(j + 255),1) + frameT(1:256,i,1);
        frame(j:(j + 255),2) = frame(j:(j + 255),2) + frameT(1:256,i,2);
    end
    frame(1:448,1)=0;
    frame(1:448,2)=0;
    frame(1601:2048,1)=0;
    frame(1601:2048,2)=0;


end