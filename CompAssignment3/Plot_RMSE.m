function Plot_RMSE(cpRef, maxN, chord, vInf)
    % Size of of the reference vector
    refSize = length(cpRef);
    
    %% TODO:
    %       (1) Need to create vector of reference points that is the same
    %           size as current test. 
    %               (a) The selected control points must be close to the
    %                   points in the current vector.
    %       (2) Compute RMSE with the smaller version of the reference
    %           vector, and the current test vector.
    %%
    figure; hold on;
    panels = 10:10:1000;
    for N = panels
        
        % cpi is the test vector 
        [x0012i, y0012i] = NACA_Airfoils(0,0,0.12,chord,N+1);
        [~,cpi, ~] = Vortex_Panel(x0012i,y0012i,vInf,0,false);
        
        step = floor(refSize/N);
        cpRefControl = zeros(N,1);
        
        for index = 1:N
            cpRefControl(index) = cpRef(index*step);
        end
        
 
        
        RSME = sqrt(sum((cpRefControl-cpi).^2)/N);
        plot(N, RSME, 'ko', 'MarkerFaceColor', 'r');
    end
    





end