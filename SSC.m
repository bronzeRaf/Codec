%Sequence Segmentation Control

%Calculates the frame type
function frameType = SSC(frameT, nextFrameT, prevFrameType)
    if strcmp(prevFrameType, 'LSS')     %LONG START --> EIGHT SHORT
        frameType = 'ESH';
    elseif strcmp(prevFrameType, 'LPS') %LONG STOP --> ONLY LONG
        frameType = 'OLS';
    elseif strcmp(prevFrameType, 'OLS') %ONLY LONG --> DECIDE
        eight1 = next8(nextFrameT(:,1));
        eight2 = next8(nextFrameT(:,2));
        if eight1 == 1                   %check next frame if is eight short
            frameType1 = 'LSS';
        else                            %if not
            frameType1 = 'OLS';
        end
        if eight2 == 1                  %check next frame if is eight short
            frameType2 = 'LSS';
        else                            %if not
            frameType2 = 'OLS';
        end
        if strcmp(frameType1, 'LSS') || strcmp(frameType2, 'LSS')
            frameType = 'LSS';
        else
            frameType = 'OLS';
        end
    elseif strcmp(prevFrameType, 'ESH') %EIGHT SHORT --> DECIDE
        eight1 = next8(nextFrameT(:, 1));
        eight2 = next8(nextFrameT(:, 2));
        if eight1 == 1                  %check next frame if is eight short
            frameType1 = 'ESH';
        else                            %if not
            frameType1 = 'LPS';
        end
        if eight2 == 1                  %check next frame if is eight short
            frameType2 = 'ESH';
        else                            %if not
            frameType2 = 'LPS';
        end
        if strcmp(frameType1, 'ESH') || strcmp(frameType2, 'ESH')
            frameType = 'ESH';
        else
            frameType = 'LPS';
        end
    else                                %unknown type of previous frame... ERROR occured
        frameType='OLS';
        disp('UNKNOWN PREVIOUS FRAME TYPE, Returning only long sequence...')
    end
    
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
        
        %at last find the attack values of the current subframe
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