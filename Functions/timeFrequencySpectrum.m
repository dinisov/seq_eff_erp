function R = timeFrequencySpectrum(R, blocks, n_back)
%timeFrequencyAnalysis Summary of this function goes here
%   Detailed explanation goes here
    cwtERPs = [];
    
    %300 Hz seems like a sensible upper limit
    for s = 1:0.5*2^n_back
        try
            [wt,f] = cwt(R.meanERPs(:,s),blocks(1).resampleFreq,'FrequencyLimits',[0 300]);
            cwtERPs(:,s,:) = wt; %#ok<AGROW>
        catch
            ['fft calculation failure']
            crash = yes
        end
    end
    
    magnitudeSEs = abs(cwtERPs);
    phaseSEs = angle(cwtERPs);
     
    R.magnitudeSEs = magnitudeSEs;
    R.phaseSEs = phaseSEs;
    R.f = f;

    %Repeat for Isomers
    if isfield( R , 'ISOMER' )
        theseIsom = fieldnames(R.ISOMER);
        for iso = 1:size(theseIsom,1)
            thisIsom = theseIsom{iso};

            %Borrow above
            cwtERPs = [];
    
            %300 Hz seems like a sensible upper limit
            for s = 1:0.5*2^n_back
               %[wt,f] = cwt(R.meanERPs(:,s),blocks(1).resampleFreq,'FrequencyLimits',[0 300]);
               try
                    [wt,f] = cwt(R.ISOMER.(thisIsom).meanERPs(:,s),blocks(1).resampleFreq,'FrequencyLimits',[0 300]);
                    cwtERPs(:,s,:) = wt; %#ok<AGROW>
               catch
                    ['## second stage FFT failure ##']
                    crash = yes
               end
            end
            
            magnitudeSEs = abs(cwtERPs);
            phaseSEs = angle(cwtERPs);
             
            R.ISOMER.(thisIsom).magnitudeSEs = magnitudeSEs;
            R.ISOMER.(thisIsom).phaseSEs = phaseSEs;
            R.ISOMER.(thisIsom).f = f;
        end
    end


end