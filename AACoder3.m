%AA Coder

%Makes a struct array with all the information for each object-frame.
function AACSeq3 = AACoder3(fNameIn, fnameAACoded)
    [tempy,fs] = wavread(fNameIn);      

    %start with 1024 zeros for better coding
    y(1:1024, 1:2) = 0;
    y(1025:1024 + length(tempy), 1:2) = tempy(:, 1:2);
        
%start timer
genStart = tic;
    
    N = length(y);                  %number of samples we have
   
    M = floor( (N) / 1024 ) - 1;    %number of frames we have
    ups =  N - M * 1024;            %numbers of unframed samples
    if ups ~=0
        M = M + 1;                  %frame with unframed and zero frame in start and stop
    end
    
    %initializing the frame1,2 with zeros for faster memory allocation
    frame1 = zeros(M,2048);
    frame2 = zeros(M,2048);
    
    % Divide into frames
    for i = 1:(M - 1)
        j = (i - 1) * 1024 + 1;
        frame1(i,:) = y( j : ( j + 2047 ) , 1 );    %split channels
        frame2(i,:) = y( j : ( j + 2047 ) , 2 );
    end
    
    if ups ~=0
        %add last frame with unframed samples
        frame1(M,1:ups) = y((N - ups + 1):N,1);
        frame2(M,1:ups) = y((N - ups + 1):N,2);
        frame1(M,(ups + 1):2048) = 0;
        frame2(M,(ups + 1):2048) = 0;
    end
    
    if (rem(M,2)~=0)
        M = M + 1;                      %force M to be even for MDCT requirement
        frame1(M,1:2048) = 0;
        frame2(M,1:2048) = 0;
    end
    AACSeq3(M).winType = 'UNK';
    for i = 1:M
        AACSeq3(i).winType = 'SIN';     %prealocate struct array and choose window type
    end
    
    eight1 = next8(frame1(2,:));
    eight2 = next8(frame2(2,:));
    
    if eight1 == 1 || eight2 ==1        %we dont expect big energy in the first frame 
        AACSeq3(1).frameType = 'LSS';   %so initial frame is long start
    else
        AACSeq3(1).frameType = 'OLS';   %or only long
    end        
    frameT = [frame1(1,:) ; frame2(1,:)]';
    %transform frame
    frameF = filterbank(frameT, AACSeq3(1).frameType, AACSeq3(1).winType);
    [frameF, TNScoeffs] = TNS(frameF, AACSeq3(1).frameType);
    %first frame doesn't have previous frames so SMR stays zero
    SMR1 = 0;
    SMR2 = 0;
    %quantize
    [S1, sfc1, G1] = AACquantizer(frameF(:,1) , AACSeq3(1).frameType , SMR1);
    [S2, sfc2, G2] = AACquantizer(frameF(:,2) , AACSeq3(1).frameType , SMR2);
    %encode huffman
    [huffSec1, huffcodebook1] = encodeHuff(S1, huffLUT, forceCodebook);
    [huffSec2, huffcodebook2] = encodeHuff(S2, huffLUT, forceCodebook);
    %store values
    AACSeq3(1).chl.frameF(1:1024) = frameF(1:1024,1);
    AACSeq3(1).chr.frameF(1:1024) = frameF(1:1024,2);
    AACSeq3(1).chl.TNScoeffs(1:4,1) = TNScoeffs(1:4,1);
    AACSeq3(1).chr.TNScoeffs(1:4,1) = TNScoeffs(1:4,2);
    AACSeq3(1).chl.T = SMR1;
    AACSeq3(1).chr.T = SMR2;
    AACSeq3(1).chl.G = G1;
    AACSeq3(1).chr.G = G2;
    AACSeq3(1).chl.sfc = sfc1;
    AACSeq3(1).chr.sfc = sfc2;
    AACSeq3(1).chl.stream = huffSec1;
    AACSeq3(1).chr.stream = huffSec2;
    AACSeq3(1).chl.codebook = huffcodebook1;
    AACSeq3(1).chr.codebook = huffcodebook2;
    
    prevFrameType = AACSeq3(1).frameType;
    
    %start coding the frames
    for i = 2:M-1
        frameT = [frame1(i,:) ; frame2(i,:)]';
        %get frameType
        AACSeq3(i).frameType = SSC(frameT,[frame1(i+1,:) ; frame2(i+1,:)]',prevFrameType);
        %transform frame
        frameF = filterbank(frameT, AACSeq3(i).frameType, AACSeq3(i).winType);
        if i>2
            %normaly take SMR
            SMR1 = psycho(frame1(i,:), AACSeq3(i).frameType, frame1(i-1,:),frame1(i-2,:));
            SMR2 = psycho(frame2(i,:), AACSeq3(i).frameType, frame2(i-1,:),frame2(i-2,:));
        else
            %second frame doesn't have previous 2 frames so SMR stays zero
            SMR1 = 0;
            SMR2 = 0;
        end
        %TNS
        [frameF, TNScoeffs] = TNS(frameF, AACSeq3(i).frameType);
        
        if strcmp(AACSeq3(i).frameType, 'ESH')      %EIGHT SHORT frame
            %quantize
            [S1, sfc1, G1] = AACquantizer(frameF(:,:,1) , AACSeq3(i).frameType , SMR1);
            [S2, sfc2, G2] = AACquantizer(frameF(:,:,2) , AACSeq3(i).frameType , SMR2);
        else    %NOT ESH frames
            %quantize
            [S1, sfc1, G1] = AACquantizer(frameF(:,1) , AACSeq3(i).frameType , SMR1);
            [S2, sfc2, G2] = AACquantizer(frameF(:,2) , AACSeq3(i).frameType , SMR2);
        end
        %encode huffman
        [huffSec1, huffcodebook1] = encodeHuff(S1, huffLUT, forceCodebook);
        [huffSec2, huffcodebook2] = encodeHuff(S2, huffLUT, forceCodebook);
        
        %store values (esh has different sizes of matrices)
        if strcmp(AACSeq3(i).frameType, 'ESH')      %EIGHT SHORT frame
            AACSeq3(i).chl.frameF(1:128,1:8) = frameF(1:128,1:8,1);
            AACSeq3(i).chr.frameF(1:128,1:8) = frameF(1:128,1:8,2);
            AACSeq3(i).chl.TNScoeffs(1:4,1:8) = TNScoeffs(1:4,1:8,1);
            AACSeq3(i).chr.TNScoeffs(1:4,1:8) = TNScoeffs(1:4,1:8,2);
        else                                        %other frame types
            AACSeq3(i).chl.frameF(1:1024) = frameF(1:1024,1);
            AACSeq3(i).chr.frameF(1:1024) = frameF(1:1024,2);
            AACSeq3(i).chl.TNScoeffs(1:4,1) = TNScoeffs(1:4,1);
            AACSeq3(i).chr.TNScoeffs(1:4,1) = TNScoeffs(1:4,2);
        end
        AACSeq3(i).chl.T = SMR1;
        AACSeq3(i).chr.T = SMR2;
        AACSeq3(i).chl.G = G1;
        AACSeq3(i).chr.G = G2;
        AACSeq3(i).chl.sfc = sfc1;
        AACSeq3(i).chr.sfc = sfc2;
        AACSeq3(i).chl.stream = huffSec1;
        AACSeq3(i).chr.stream = huffSec2;
        AACSeq3(i).chl.codebook = huffcodebook1;
        AACSeq3(i).chr.codebook = huffcodebook2;
        
        prevFrameType = AACSeq3(i).frameType;
    end
    %for the last frame we dont expect high energy if is not nessecary so
    %it not eight short probably. And sure it is not long start cause we
    %dont have energy in the next frame
    if strcmp(AACSeq3(M-1).frameType, 'LPS') || strcmp(AACSeq3(M-1).frameType, 'OLS')
        AACSeq3(M).frameType = 'OLS';               %LONG STOP / ONLY   --> ONLY LONG
    elseif strcmp(AACSeq3(M-1).frameType, 'LSS')
        AACSeq3(M).frameType = 'ESH';               %LONG START         --> EIGHT SHORT
    elseif strcmp(AACSeq3(M-1).frameType, 'ESH')
        AACSeq3(M).frameType = 'LPS';               %EIGHT SHORT        --> LONG STOP
    end
    frameT = [frame1(M,:) ; frame2(M,:)]';
    %transform frame
    frameF = filterbank(frameT, AACSeq3(M).frameType, AACSeq3(M).winType);
    [frameF, TNScoeffs] = TNS(frameF, AACSeq3(M).frameType);
    %get SMR normaly
    SMR1 = psycho(frame1(i,:), AACSeq3(M).frameType, frame1(M-1,:),frame1(M-2,:));
    SMR2 = psycho(frame2(i,:), AACSeq3(M).frameType, frame2(M-1,:),frame2(M-2,:));
    
    if strcmp(AACSeq3(M).frameType, 'ESH')      %EIGHT SHORT frame
        %quantize
    	[S1, sfc1, G1] = AACquantizer(frameF(:,:,1) , AACSeq3(M).frameType , SMR1);
        [S2, sfc2, G2] = AACquantizer(frameF(:,:,2) , AACSeq3(M).frameType , SMR2);
    else    %not ESH frame
        %quantize
        [S1, sfc1, G1] = AACquantizer(frameF(:,1) , AACSeq3(M).frameType , SMR1);
        [S2, sfc2, G2] = AACquantizer(frameF(:,2) , AACSeq3(M).frameType , SMR2);
    end
    %encode huffman
    [huffSec1, huffcodebook1] = encodeHuff(S1, huffLUT, forceCodebook);
    [huffSec2, huffcodebook2] = encodeHuff(S2, huffLUT, forceCodebook);
        
    %store values
    if strcmp(AACSeq3(M).frameType, 'ESH')      %EIGHT SHORT frame
            AACSeq3(M).chl.frameF(1:128,1:8) = frameF(1:128,1:8,1);
            AACSeq3(M).chr.frameF(1:128,1:8) = frameF(1:128,1:8,2);
            AACSeq3(M).chl.TNScoeffs(1:4,1:8) = TNScoeffs(1:4,1:8,1);
            AACSeq3(M).chr.TNScoeffs(1:4,1:8) = TNScoeffs(1:4,1:8,2);
    else                                        %other frame types
            AACSeq3(M).chl.frameF(1:1024) = frameF(1:1024,1);
            AACSeq3(M).chr.frameF(1:1024) = frameF(1:1024,2);
            AACSeq3(M).chl.TNScoeffs(1:4,1) = TNScoeffs(1:4,1);
            AACSeq3(M).chr.TNScoeffs(1:4,1) = TNScoeffs(1:4,2);
    end
    
    AACSeq3(M).chl.T = SMR1;
 	AACSeq3(M).chr.T = SMR2;
 	AACSeq3(M).chl.G = G1;
	AACSeq3(M).chr.G = G2;
   	AACSeq3(M).chl.sfc = sfc1;
  	AACSeq3(M).chr.sfc = sfc2;
   	AACSeq3(M).chl.stream = huffSec1;
  	AACSeq3(M).chr.stream = huffSec2;
  	AACSeq3(M).chl.codebook = huffcodebook1;
 	AACSeq3(M).chr.codebook = huffcodebook2;
    
%stop timer print
fullTime = toc(genStart)
%save workspace
save ('fnameAACoded.mat', 'AACSeq3');
end


%checks if the frame of matrix argument will be eight short or not
%used for predict the i+1 frame type (helps SSC)
function eight = next8(nextFrameT)
    eight = 1;                      %initialing in 1 to easy return
    %making filter parameters
    a = [1 -0.5095];
    b = [0.7548 -0.7548];
    sm = 0;                        %initializing the energy sum for the average
    %apply the filter to the samples of current frame for each channel (stereo)
    filtered1(1:2048) = filter(b,a,nextFrameT(1:2048));
    
    for i = 1:8
        count = 448 + i * 128 + 1;  %go to the start of the ith subframe
        %sum of the square of samples of current subframe
        s = sum(filtered1( count : (count + 127)).^2);
        
        %at last find the attack values of current subframe
        if i ~= 1
            ds = s /(sm /(i - 1));
            if (ds > 10)&& (s > 0.001) %if one of the channels is eight short the frame is eght short too
                return              %return cause current frame fixed as eight short

%if one channel is long start and the other is long stop the frame is eight short
%but we couldn't get for the same frame long start and long stop to its
%channels cause previous frame type is the same for both
            end
        end
        %current subframe will take part to the next average
        sm = sm + s;
    end
    
    eight=0;
end