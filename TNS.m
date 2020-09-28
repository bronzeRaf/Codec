%Temporal Noise Shaping

%Transforms MDCT to a new group of coeffs (with the same size) without
%periodicity and also returns the filter coeffs
function [ frameFout , TNScoeffs ] = TNS( frameFin , frameType)

    
%_____________________________________________________________________________________________
%_____________________________________________________________________________________________

    clear Sw;
    clear a;
    clear TNScoeffs;
    Sw = calcSw(frameType,frameFin);

    if(strcmp(frameType,'OLS') || strcmp(frameType,'LSS') || strcmp(frameType,'LPS'))   %NOT ESH
    
        TNScoeffs = zeros(8,2);
        %this is only for the final zero values in the whole frame that gives
        %NaN at Xw because of the 0/0 division
        Xw = frameFin./Sw;
    
        %Second step
        [a(1,:),g] = lpc(Xw(:,1),4);  %Calculates the coefficiency vector
        [a(2,:),g] = lpc(Xw(:,2),4);
        a=a';
        tmp=0;
        %if Xw(:,:)==0 or non finite a will have NaNs and we will receive error
        %this is gonna hapen in the last frame if its its covered only with 0
        for i=1:5
            if isnan(a(i,1))||isnan(a(i,2))
                a(2:5,:) = 0;
                break;
            end
        end
        
    
        while(tmp~=2)
            t=0;
        
            b(:,1) = roots([1 -a(2:end,1)']);
            b(:,2) = roots([1 -a(2:end,2)']);
        
            %checks if the filter will be stable
            for i = 1:4
                if(abs(b(i,1))>1)
                    b(i,1) = b(i,1)/(1.1*abs(b(i,1))); %moves the root inside the circle, radius=1 
                    t=t+1;
                end
                if(abs(b(i,2))>1)
                    b(i,2) = b(i,2)/(1.1*abs(b(i,2)));
                    t = t+1;
                end
            end
            %Get the new factors
            a(:,1) = poly(b(:,1));
            a(:,2) = poly(b(:,2));
            a(2:end,2) = -a(2:end,1);
            a(2:end,2) = -a(2:end,2);
    
            %convert to binary-return TNScoeffs
            for i = 2:5
                if(a(i,1) > 0.8) %moves the factors in the coding range
                    a(i,1) = 0.8;
                elseif(a(i,1) < -0.7)
                    a(i,1) = -0.7;
                end
                if(a(i,2) > 0.8) %moves the factors in the coding range
                a(i,2) = 0.8;
                elseif(a(i,2) < -0.7)
                a(i,2) = -0.7;
                end
            end
    
            if(t==0)
                tmp = tmp+1;
            else
                tmp = 0;
            end
    
        end
    
        %use step of 0.1
        a = a.*10;
        %convert to the int part
         a = fix(a);
    
    
        TNScoeffs(1:4,1) = a(2:5,1);
        TNScoeffs(1:4,2) = a(2:5,2);
    
        %devide by 10, get the original factors and filter the signal
        a = a./10;
        %filter the signal
        frameFout(:,1) = filter([1 -a(2,1) -a(3,1) -a(4,1) -a(5,1)],1,frameFin(:,1)); 
        frameFout(:,2) = filter([1 -a(2,2) -a(3,2) -a(4,2) -a(5,2)],1,frameFin(:,2)); 
    
    
    elseif(strcmp(frameType,'ESH'))     %ESH
    
        frameFout = zeros(128,8,2);
        TNScoeffs = zeros(4,8,2);
        Xw = frameFin./Sw;
        %Second step
        for k = 1:8
            [a(1,1:5),g] = lpc(Xw(1:128,k,1),4);
            [a(2,1:5),g] = lpc(Xw(1:128,k,2),4);
            a = a';
        
            tmp=0;
    
            while(tmp~=2)
                t=0;
                        
                b(:,1) = roots([1 -a(2:end,1)']);
                b(:,2) = roots([1 -a(2:end,2)']);
    
                %checks if the filter will be stable
                for i = 1:4
                    if(abs(b(i,1)) > 1)
                        b(i,1) = b(i,1)/(1.1*abs(b(i,1))); %moves the root inside the circle, radius=1 
                        t = t+1;
                    end
                    if(abs(b(i,2)) > 1)
                        b(i,2) = b(i,2)/(1.1*abs(b(i,2)));
                        t = t+1;
                    end
                end
                %Get the new factors
                a(:,1) = poly(b(:,1));
                a(:,2) = poly(b(:,2));
                a(2:end,1) = -a(2:end,1);
                a(2:end,2) = -a(2:end,2);
            
                %convert to binary-return TNScoeffs
                for i = 2:5
                    if(a(i,1) > 0.8) %moves the factors in the coding range
                        a(i,1) = 0.8;
                    elseif(a(i,1) < -0.7)
                        a(i,1) = -0.7;
                    end
                    if(a(i,2) > 0.8) %moves the factors in the coding range
                        a(i,2) = 0.8;
                    elseif(a(i,2) < -0.7)
                        a(i,2) = -0.7;
                    end
                end
    
                if(t==0)
                    tmp = tmp+1;
                else
                    tmp = 0;
                end
    
            end
        
            %use step of 0.1
            a = a.*10;
            %convert to the int part
            a = fix(a); 
        
            TNScoeffs(:,k,1) = a(2:5,1);
            TNScoeffs(:,k,2) = a(2:5,2);
        
            %devide by 10, get the original factors and filter the signal
            a = a./10;
            frameFout(:,k,1) = filter([1 -a(2,1) -a(3,1) -a(4,1) -a(5,1)],1,frameFin(:,k,1));
            frameFout(:,k,2) = filter([1 -a(2,2) -a(3,2) -a(4,2) -a(5,2)],1,frameFin(:,k,2));
        
        
        
        end
    
end

end

%Calculates and normalizes the Sw coeffs of the given frame
function Sww = calcSw(frameType, frameFin)
    %bands in Barks for long frames
    bl=[0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 41 44 47 50 53 56 59 62 66 70 74 78 82 87 92 97 103 109 116 123 131 139 148 158 168 179 191 204 218 233 249 266 284 304 325 348 372 398 426 457 491 528 568 613 663 719 782 854 938];
    %bands in barks for eight short frames
    bs=[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 19 21 23 25 27 29 31 34 37 40 43 46 50 54 58 63 68 74 80 87 95 104 114 126];


    clear Sww;
    clear P;
    P = calcP(frameType, frameFin);
    if(strcmp(frameType,'OLS') || strcmp(frameType,'LSS') || strcmp(frameType,'LPS'))   %NOT ESH
        Sww = zeros(1024,2);
        for j = 1:68  %number of bands
            for i = (bl(j)+1):bl(j+1)
                Sww(i,1) = sqrt(P(j,1));
                Sww(i,2) = sqrt(P(j,2));
            end
        end
        Sww(939:1024,1) = sqrt(P(69,1));
        Sww(939:1024,2) = sqrt(P(69,2));
    
        for i = 1023:-1:1   %Normalize data
            Sww(i,1) = (Sww(i,1)+Sww(i+1,1))/2;
            Sww(i,2) = (Sww(i,2)+Sww(i+1,2))/2;
        end
        for i = 2:1024      %Normalize data
            Sww(i,1) = (Sww(i,1) + Sww(i-1,1))/2;
            Sww(i,2) = (Sww(i,2) + Sww(i-1,2))/2;
        end
    else                %ESH case
        Sww = zeros(128,8,2);
        for k = 1:8
            for j = 1:41  %number of bands
                for i = (bs(j)+1):bs(j+1)
                    Sww(i,k,1) = sqrt(P(j,k,1));
                    Sww(i,k,2) = sqrt(P(j,k,2));
                end
            end
            Sww(128,k,1) = sqrt(P(42,k,1));
            Sww(128,k,2) = sqrt(P(42,k,2));
        end
    
        for k = 1:8 %Normalize data
            for i = 127:-1:1
                Sww(i,k,1) = (Sww(i,k,1)+Sww(i,k,1))/2;
                Sww(i,k,2) = (Sww(i,k,2)+Sww(i,k,2))/2;
            end
        end
    
        for k = 1:8 %Normalize data
            for i = 2:128
                Sww(i,k,1) = (Sww(i,k,1)+Sww(i-1,k,1))/2;
                Sww(i,k,2) = (Sww(i,k,2)+Sww(i-1,k,2))/2;
            end
        end
    
    end
end

%Calculating the energy spectrum Pi(j)
function P = calcP(frameType, frameF)
%bands in Barks for long frames
bl=[0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 41 44 47 50 53 56 59 62 66 70 74 78 82 87 92 97 103 109 116 123 131 139 148 158 168 179 191 204 218 233 249 266 284 304 325 348 372 398 426 457 491 528 568 613 663 719 782 854 938 1023];
%bands in barks for eight short frames
bs=[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 19 21 23 25 27 29 31 34 37 40 43 46 50 54 58 63 68 74 80 87 95 104 114 126 127];


    clear P;
    if strcmp(frameType,'ESH')  %eight short
        P = zeros(42,8,2);
        for i=1:8
            for j=1:42  %number of bands of each subframe
                for k=(bs(j)+1):1:bs(j+1)
                    P(j,i,1)=P(j,i,1)+ frameF(k,i,1)*frameF(k,i,1); %channel 1
                    P(j,i,2)=P(j,i,2)+ frameF(k,i,2)*frameF(k,i,2); %channel 2
                end
            end
        end
    else                        %other frameTypes
        P = zeros(69,2);
        for j=1:69      %number of bands
            for k=(bl(j)+1):1:bl(j+1)
                P(j,1)=P(j,1)+ frameF(k,1)*frameF(k,1);             %channel 1
                P(j,2)=P(j,2)+ frameF(k,2)*frameF(k,2);             %channel 2
            end
        end
    end
end