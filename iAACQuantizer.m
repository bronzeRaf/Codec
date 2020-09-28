%Inverter of the quantizer

%Inverts the quntized given frame to its normal MDCT form
function frameF = iAACQuantizer( S , sfc , G ,frameType )

%bands in Barks for not ESH frames
%length(bl)=70 because 1024 was added at the end
bl=[0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 41 44 47 50 53 56 59 62 66 70 74 78 82 87 92 97 103 109 116 123 131 139 148 158 168 179 191 204 218 233 249 266 284 304 325 348 372 398 426 457 491 528 568 613 663 719 782 854 938 1024];

%bands in barks for ESH frames
%length(bl)=43 because 128 was added at the end
bs=[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 19 21 23 25 27 29 31 34 37 40 43 46 50 54 58 63 68 74 80 87 95 104 114 126 128];

    if(strcmp(frameType,'OLS') || strcmp(frameType,'LSS') || strcmp(frameType,'LPS'))   %NOT ESH frames
        %inverse DPCM
        for i = 2:69
            sfc(i) = sfc(i-1)+sfc(i);
        end
        %get a
        a(:) = G - sfc(:);
        for j = 1:69    %bands
            for i = (bl(j)+1):bl(j+1)   %samples
                frameF(i) = sign(S(i))*abs(S(i))^(4/3)*2^(a(j)/4);
            end
        end
    
    elseif(strcmp(frameType,'ESH'))
        for k = 0:7   %subframes
            %inverse DPCM
            for i = 2:42
                sfc(i,k+1) = sfc(i-1,k+1)+sfc(i,k+1);
            end
            %get a
            a((k*42+1):(k+1)*42) = G(k+1)-sfc(1:42,k+1);
            for j = 1:42  %bands
                for i = (bs(j)+1):bs(j+1) %samples
                    frameF(i,k+1) = sign(S(k*128+i))*abs(S(k*128+i))^(4/3)*2^(a(k*42+j)/4);
                end
            end
        end
    end


end