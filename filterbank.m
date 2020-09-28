%Filterbank

%Gives the trasposed (with MDCT after windowing) samples of the input frame
function frameF = filterbank(frameT , frameType , winType)
    clear frameF;               %clear to get sure the correct size matrix
    clear w;
    %get the wondow we need
    w = windowMaker(frameType,winType);
    if strcmp(frameType, 'ESH')
        frameT8 = subframeMaker(frameT);
        frameT8(1:256, 1:8, 1) =  w(1:256, 1:8).*frameT8(1:256, 1:8, 1);    %multiply window to our frame channel 1
        frameT8(1:256, 1:8, 2) =  w(1:256, 1:8).*frameT8(1:256, 1:8, 2);    %multiply window to our frame channel 2
        frameF(1:128,1:8,1) = mdct4(frameT8(:, :, 1));
        frameF(1:128,1:8,2) = mdct4(frameT8(:, :, 2));
    else
        frameT(1:2048, 1) = w(1:2048).*frameT(1:2048,1);                    %multiply window to our frame channel 1
        frameT(1:2048, 2) = w(1:2048).*frameT(1:2048,2);                    %multiply window to our frame channel 2
        frameF(1:1024,1) = mdct4(frameT(:, 1));
        frameF(1:1024,2) = mdct4(frameT(:, 2));
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
        disp('error!! unkown window type, unable to set window')
        return
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
        w = big;
        return %nothing to change
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

%Cuts the frame matrix into a matrix with the 8 subframes as we need for
%eight short frame type.
function frame = subframeMaker(frameT)
    count = 449;
    frame = zeros(256, 8, 2);
    
    for i = 1:8
        frame(1:256,i,1) = frameT(count:count + 255,1);
        frame(1:256,i,2) = frameT(count:count + 255,2);
        count = count + 128;
    end
end

%Transforms the given frame to its MDCT.
%Taken from http://www.ee.columbia.edu/~marios/mdct/mdct4.m
%Credits to Marios Athineos
function y = mdct4(x)

    [flen,fnum] = size(x);
    % Make column if it's a single row
    if (flen==1)
        x = x(:);
        flen = fnum;
        fnum = 1;
    end
    % Make sure length is multiple of 4
    if (rem(flen,4)~=0)
        error('MDCT4 defined for lengths multiple of four.');
    end

    % We need these for furmulas below
    N     = flen; % Length of window
    M     = N/2;  % Number of coefficients
    N4    = N/4;  % Simplify the way eqs look
    sqrtN = sqrt(N);

    % Preallocate rotation matrix
    % It would be nice to be able to do it in-place but we cannot
    % cause of the prerotation.
    rot = zeros(flen,fnum);

    % Shift
    t = (0:(N4-1)).';
    rot(t+1,:) = -x(t+3*N4+1,:);
    t = (N4:(N-1)).';
    rot(t+1,:) =  x(t-N4+1,:);
    clear x;

    % We need this twice so keep it around
    t = (0:(N4-1)).';
    w = diag(sparse(exp(-j*2*pi*(t+1/8)/N)));

    % Pre-twiddle
    t = (0:(N4-1)).';
    c =   (rot(2*t+1,:)-rot(N-1-2*t+1,:))...
    -j*(rot(M+2*t+1,:)-rot(M-1-2*t+1,:));
    % This is a really cool Matlab trick ;)
    c = 0.5*w*c;
    clear rot;

    % FFT for N/4 points only !!!
    c = fft(c,N4);

    % Post-twiddle
    c = (2/sqrtN)*w*c;

    % Sort
    t = (0:(N4-1)).';
    y(2*t+1,:)     =  real(c(t+1,:));
    y(M-1-2*t+1,:) = -imag(c(t+1,:));
    
end