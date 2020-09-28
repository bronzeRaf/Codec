%Inverter of filterbank part for the decoder

%Inverts the given MDCT coefficients to the normal signal samples
function frameT = iFilterbank(frameF, frameType, winType)
    clear frameT;               %clear to get sure the correct size matrix
    clear w;
    %get the wondow we need
    w = windowMaker(frameType,winType);
    if strcmp(frameType, 'ESH')
        frameF1(1:128,1:8) = frameF(1:128,1:8,1);       %split channel 1
        frameF2(1:128,1:8) = frameF(1:128,1:8,2);       %split channel 2
        frameT(1:256,1:8,1) = imdct4(frameF1);
        frameT(1:256,1:8,2) = imdct4(frameF2);
        frameT(1:256,1:8,1) = w(1:256,1:8).*frameT(1:256,1:8,1);   %multiply samples with the window (2 channels)
        frameT(1:256,1:8,2) = w(1:256,1:8).*frameT(1:256,1:8,2);            
    else
        frameF1(1:1024) = frameF(1:1024,1);             %split channel 1
        frameF2(1:1024) = frameF(1:1024,2);             %split channel 2
        frameT(1:2048,1) = imdct4(frameF1);
        frameT(1:2048,2) = imdct4(frameF2);
        frameT(1:2048,1) = w(1:2048).*frameT(1:2048,1);         %multiply samples with the window (2 channels)
        frameT(1:2048,2) = w(1:2048).*frameT(1:2048,2);
    end
end

%Prepares the window (or the matrix of windows in case of eight short
%frame) that we need to multiply our frame. Works for every frame type and
%every window type
function w = windowMaker(frameType,winType)
    clear w;        %clear w to have correct size matrix
    if strcmp(winType, 'KBD')
        big = kaiser(2048,6);
        small = kaiser(256,4);
    elseif strcmp(winType, 'SIN')
        big = sineWindow(2048);
        small = sineWindow(256);
    else
        error('unkown window type, unable to set window')
    end
    %check frame type to set window geometry 
    if strcmp(frameType, 'LSS')         %LONG START FRAME WINDOW
        w = big;
        w(1025:1472) = 1;
        w(1473:1600) = small(129:256);
        w(1601:2048) = 0;
    
    elseif strcmp(frameType, 'LPS')     %LONG STOP FRAME WINDOW
        w = big;
        w(1:448) = 0;
        w(449:576) = small(1:128);
        w(577:1024) = 1;
        
    elseif strcmp(frameType, 'OLS')     %ONLY LONG FRAME WINDOW
        w = big;%nothing to change
        
    elseif strcmp(frameType, 'ESH')     %EIGHT SHORT FRAME WINDOW
        w = zeros(256, 8);
        for i = 1:8
            w(1:256,i) = small(1:256);
        end
    else
        error('error!! unkown frame type, unable to set window')
    end
    
end

%Makes a N sized sinusoid window
%Taken from http://www.ee.columbia.edu/~marios/mdct/sinewin.m 
%Credits to Marios Athineos
function y = sineWindow(N)
    x = (0:(N-1)).';
    y = sin(pi*(x+0.5)/N);
end

%Converts the MDCT Coefficients of x input back to its time frame 
%Taken from http://www.ee.columbia.edu/~marios/mdct/imdct4.m
%Credits to Marios Athineos
function y = imdct4(x)

    [flen,fnum] = size(x);
    % Make column if it's a single row
    if (flen==1)
        x = x(:);
        flen = fnum;
        fnum = 1;
    end

    % We need these for furmulas below
    N     = flen;
    M     = N/2;
    twoN  = 2*N;
    sqrtN = sqrt(twoN);

    % We need this twice so keep it around
    t = (0:(M-1)).';
    w = diag(sparse(exp(-j*2*pi*(t+1/8)/twoN)));

    % Pre-twiddle
    t = (0:(M-1)).';
    c = x(2*t+1,:) + j*x(N-1-2*t+1,:);
    c = (0.5*w)*c;

    % FFT for N/2 points only !!!
    c = fft(c,M);

    % Post-twiddle
    c = ((8/sqrtN)*w)*c;

    % Preallocate rotation matrix
    rot = zeros(twoN,fnum);

    % Sort
    t = (0:(M-1)).';
    rot(2*t+1,:)   = real(c(t+1,:));
    rot(N+2*t+1,:) = imag(c(t+1,:)); 
    t = (1:2:(twoN-1)).';
    rot(t+1,:) = -rot(twoN-1-t+1,:);

    % Shift
    t = (0:(3*M-1)).';
    y(t+1,:) =  rot(t+M+1,:);
    t = (3*M:(twoN-1)).';
    y(t+1,:) = -rot(t-3*M+1,:);
    
end
