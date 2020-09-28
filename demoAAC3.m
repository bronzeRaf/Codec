%SNR

%Calculates the SNR, the compression and the bitrate of the signal after coding and decoding
function [SNR, bitrate, compression] = demoAAC3(fNameIn, fNameOut, frameAACoded)
    [x,fs1] = wavread(fNameIn);
    [tempy,fs2] = wavread(fNameOut);
    
    
    y(:,:)=tempy(1:length(x),:); %work off unwanted sample to equalize size of x,y 
    

    
    vare = mean((x-y).^2);      %variance of error
    varx = mean(x.^2);          %variance of the original signal
    SNR = 10 * log10(varx/vare);
    
    input=whos('y');        %get the list of variables for y (we need its bytes)
    load (fNameAACoded, 'AACSeq3'); %Load .mat workspace with the LEVEL 3 struct sequence
    output=whos('AACSeq3'); %get the list of variables for AACSeq3 (we need its bytes)
    
    %compression
    compression=input.bytes/output.bytes*100;
    
    time=input.size(1)/48000;
    bitrate(1)=input.bytes*8/time;  %Bitrate wav
    bitrate(2)=output.bytes*8/time; %Bitrate aac
end