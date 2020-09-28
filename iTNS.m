%inverter of TNS part of the devoder

%Gives back the MCDT coeffs for each frame
function frameFout = iTNS(frameFin, frameType, TNScoeffs)

    if(strcmp(frameType,'OLS') || strcmp(frameType,'LSS') || strcmp(frameType,'LPS'))   %NOT ESH
        TNScoeffs = TNScoeffs./10; %TNScoeffs matrix contains the filter factors multiplied by 10
        %filter the frame
        frameFout(:,1) = filter(1 , [1 -TNScoeffs(1,1) -TNScoeffs(2,1) -TNScoeffs(3,1) -TNScoeffs(4,1)] , frameFin(:,1));
        frameFout(:,2) = filter(1 , [1 -TNScoeffs(1,2) -TNScoeffs(2,2) -TNScoeffs(3,2) -TNScoeffs(4,2)] , frameFin(:,2));
    
    elseif(strcmp(frameType,'ESH'))                     %ESH
        frameFout = zeros(128,8,2);
        for k = 1:8
            TNScoeffs = TNScoeffs./10;
            %filter the frame
            frameFout(1:128,k,1) = filter(1 , [1 -TNScoeffs(1,1) -TNScoeffs(2,1) -TNScoeffs(3,1) -TNScoeffs(4,1)] , frameFin(1:128,k,1));
            frameFout(1:128,k,2) = filter(1 , [1 -TNScoeffs(1,2) -TNScoeffs(2,2) -TNScoeffs(3,2) -TNScoeffs(4,2)] , frameFin(1:128,k,2));
        end
    end
end