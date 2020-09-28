%Quantizer

%Quantizes the MDCT given frame and gives the parameters for the inversion
%of the quantization
function [ S , sfc , G ] = AACquantizer( frameF , frameType , SMR )

%bands in Barks for not ESH frames
%length(bl)=70 because 1024 (denotes the end of array) was added at the end
bl = [0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 41 44 47 50 53 56 59 62 66 70 74 78 82 87 92 97 103 109 116 123 131 139 148 158 168 179 191 204 218 233 249 266 284 304 325 348 372 398 426 457 491 528 568 613 663 719 782 854 938 1024];

%bands in barks for ESH frames
%length(bl)=43 because 128 was added at the end
bs = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 19 21 23 25 27 29 31 34 37 40 43 46 50 54 58 63 68 74 80 87 95 104 114 126 128];

    if(strcmp(frameType,'OLS') || strcmp(frameType,'LSS') || strcmp(frameType,'LPS'))   %NOT ESH frames
   
        %find f
        f = zeros(69);
        S = zeros(1024);
        Xhat = zeros(1024);
        for j = 1:69
            for i = (bl(j)+1):bl(j+1)
                f(j) = f(j)+sqrt(abs(frameF(i)));
            end
        end
        a = ones(69);
        %find a
        for j = 1:69
            a(j) = floor(8.8585*(log10(6.75*SMR(j))-log10(f(j))));
        end
        %find starting S and X(k)
        for j = 1:69
            for i = (bl(j)+1):bl(j+1)
                S(i) = sign(frameF(i))*fix(((abs(frameF(i))/(2^(a(j)/4)))^(3/4)+0.4054));
            
                Xhat(i) = sign(S(i))*abs(S(i))^(4/3)*2^(a(1)/4);
            end
        end
    
        for j = 1:69
            
            flag = 0;
            temp = zeros(2);
            while(1)
            
                if((max(a(:))-min(a(:))) > 60)
                    break
                end
                %noise for each band
                noise((bl(j)+1):bl(j+1)) = frameF((bl(j)+1):bl(j+1))-Xhat((bl(j)+1):bl(j+1));
                Pnoise = 0;
                for i = (bl(j)+1):bl(j+1)
                    Pnoise = Pnoise+noise(i)^2;     %Power of quantization error
                end
                PnoisedB = 10*log10(Pnoise);        %convert to dB
                temp(1) = sign(PnoisedB-SMR(j));    % temp(1) has the new value
                        
                if(temp(1) == -temp(2))
                    flag = flag+1;
                end
                if(flag == 2)
                    break
                end
            
                if(PnoisedB < SMR(j))
                    a(j) = a(j)+1;
                else
                    a(j) = a(j)-1;
                end
            
                temp(2) = temp(1);  %temp(2)has the old value
            
                %recalculate S and X(k)
                for i = (bl(j)+1):bl(j+1)
                    S(i) = sign(frameF(i))*fix((abs(frameF(i)/(2^(a(j)/4)))^(3/4)+0.4054));
                    Xhat(i) = sign(S(i))*abs(S(i))^(4/3)*2^(a(j)/4);
                end
            
            end
       
        end
    
    
        G = max(a(:));
    
        sfc(:) = G(1)-a(:);
        tmp(:) = sfc(:);
        %Code with DPCM
        for i = 2:69
            sfc(i) = tmp(i)-tmp((i-1));
        end
    
    elseif(strcmp(frameType,'ESH')) %ESH frame
        S = zeros(1024);
        Xhat = zeros(1024);
        a = ones(336);

        for k = 0:7
            %Calculate f
            f = zeros(42);
            %find f
            for j = 1:42
                for i = (bs(j)+1):bs(j+1)
                    f(j) = f(j)+sqrt(abs(frameF((k*128+i))));
                end
            end
            %find a
            a((k*42+1):((k+1)*42)) = floor(8.8585*(log10(6.75*SMR((k*42+1):((k+1)*42)))-log10(f(:))));
            %find starting S and X(k)
            for j = 1:42
                for i = (bs(j)+1):bs(j+1)
                    S((k*128+i)) = sign(frameF((k*128+i)))*fix(((abs(frameF((k*128+i)))/(2^(a((k*42+j))/4)))^(3/4)+0.4054));
            
                    Xhat((k*128+i)) = sign(S((k*128+i)))*abs(S((k*128+i)))^(4/3)*2^(a((k*42+j))/4);
                end
            end
        end
    
        for k = 0:7
            for j = 1:42
                
                flag = 0;
                temp = zeros(2);
                while(1)
                    if((max(a((k*42+1):((k+1)*42)))-min(a((k*42+1):((k+1)*42)))) > 60)
                        break
                    end
                    %noise for each band
                    noise((bs(j)+1):bs(j+1)) = frameF((k*128+bs(j)+1):(k*128+bs(j+1)))-Xhat((k*128+bs(j)+1):(k*128+bs(j+1)));
                    Pnoise = 0;
                    for i = (bs(j)+1):bs(j+1)
                        Pnoise = Pnoise+noise(i)^2;         %Power of quantization error
                    end
                    PnoisedB = 10*log10(Pnoise);            %convert to dB
                    temp(1) = sign(PnoisedB-SMR((k*42+j))); % temp(1) has the new value
            
                    if(temp(1) == -temp(2))
                        flag = flag+1;
                    end
                    if(flag == 2)
                        break
                    end
            
                    if(PnoisedB < SMR((k*42+j)))
                        a((k*42+j)) = a((k*42+j))+1;
                    else
                        a((k*42+j)) = a((k*42+j))-1;
                    end
            
                    temp(2) = temp(1);%temp(2)has the old value
            
                    %recalculate S and X(k)
                    for i = (bs(j)+1):bs(j+1)
                        S(k*128+i) = sign(frameF(k*128+i))*fix(((abs(frameF(k*128+i))/(2^(a((k*42+j))/4)))^(3/4)+0.4054));
                        Xhat(k*128+i) = sign(S(k*128+i))*abs(S(k*128+i))^(4/3)*2^(a((k*42+j))/4);
                    end
                end
            
            end
        end
    
        G = zeros(8);
        sfc = zeros(42,8);
    
        for k = 0:7
            G(k+1) = max(a((k*42+1):((k+1)*42),1));
        
            sfc(1:42,(k+1)) = G(k+1)-a((k*42+1):(k+1)*42);
        
            tmp(:) = sfc((k*42+1):(k+1)*42);
            %Code with DPCM
            for i = 2:42
                sfc((k*42+i)) = tmp(i)-tmp((i-1));
            end
        
        end
    
    end

end