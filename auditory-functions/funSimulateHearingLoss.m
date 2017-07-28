
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%Simulating Hearing Loss functions %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ********** funSimulateHearingLoss **********
    function Threshold_dBSPL_FFT = funSimulateHearingLoss(frq)
    %
    %   usage: lcfSimulateHearingLoss(frq) 
    %      by: Ben Gurer
    %    date: 19/01/2017
    % purpose: Create hearing loss simulation masking noise
    % This masking noise aims to simualte hearing loss by raising the threshold of auditory perception.
    %   input: frequency in kHz
    
    Threshold_dBHL = createSteeplySlopingHearingLoss_dBHL(frq);
    
    RETSPL_int = getdBHLinSPL_inserts(frq*1000);
    
    Threshold_dBSPL_FFT = Threshold_dBHL + RETSPL_int;
    
    % set up variables from gramm
    x = [repmat(frq,1,3)];
    y = [Threshold_dBHL RETSPL_int Threshold_dBSPL_FFT];    
    names = [repmat({'tHL (dB HL)'},length(frq),1); repmat({'tH (dB SPL)'},length(frq),1); repmat({'tsHL (dB SPL)'},length(frq),1)];
    
%     figure
%     g = gramm('x',x,'y',y,'color',names);
%     g.geom_line()
%     % Set appropriate names for legends
%     g.set_names('x','Frequency (kHz)','y','Level (dB)','color','Threshold of:')
%     %Set figure title
%     g.set_title('Hearing loss simulation')
%     g.draw()
    end

    function RETSPL_int = getdBHLinSPL_inserts(f)
    %
    %   usage: getdBHLinSPL_inserts(f)
    %      by: Ben Gurer
    %    date: 19/01/2017
    % purpose: output the dB Hearling Level values in dB Sound Pressure Level for insert transducers for frequencies specified using the input arguement
    %
    % discription:
    % interpolates the dB HL values for input arguements of frequency from value specified by ISO 389.

    % BS EN ISO 389-2:1997
    % Acoustics - Reference zero for the calibration of audiometric equipment � Part 2: Reference equivalent threshold sound
    % pressure levels for pure tones and insert earphones.
    % Transducer: Etymotic Research ER-3A
    % Ear simulator: Occluded-ear simulator (IEC 711)
    
    RETSPL_int = zeros(1,length(f));

    RETSPL_125_8000Hz = [28.0, 24.5 21.5, 17.5, 15.5, 13.0, 9.5, 7.5, 6.0, 5.5, 5.5, 8.5, 9.5, 9.5, 11.5, 13.5, 13.0, 13.0, 15.0, 18.5, 16.0, 16.0, 15.5];

    f_hz_measured_125_8000Hz = [125, 160, 200, 250, 315, 400, 500, 630, 750, 800, 1000, 1250, 1500, 1600, 2000, 2500, 3000, 3150, 4000, 5000, 6000, 6300, 8000];

    % Acoustics - Reference zero for the calibration of audiometric equipment - Part 5: Reference equivalent threshold sound
    % pressure levels for pure tones in the frequency range 8 kHz to 16 kHz (ISO 389-5:2006)
    % Transducer: Etymotic Research ER-2b
    % Ear simulator: IEC 60711e
    % Adapter: ISO 389-2:1994

    f_hz_measured_8000_16000Hz =[8000, 9000, 10000, 11200, 12500, 14000, 16000];

    RETSPL_8000_1600Hz = [19, 16, 20, 30.5, 37, 43.5, 53];

%         figure
%         plot(f_hz_measured_125_8000Hz,RETSPL_125_8000Hz,'ko-')
%         hold on
%         plot(f_hz_measured_8000_16000Hz,RETSPL_8000_1600Hz,'ko-')

    RETSPL_int = interp1([f_hz_measured_125_8000Hz f_hz_measured_8000_16000Hz(2:end)],[RETSPL_125_8000Hz RETSPL_8000_1600Hz(2:end)],f,'spline');
%         plot(f,RETSPL_int,'r--')

%         legend('Threshold - 125 - 8000 Hz',...
%             'Threshold - 8000 - 16000Hz',...
%             'Threshold - Interpolated',...
%             'Location','best')
%         set(gca,'XLim',[min(f) max(f)])
%         title('Threshold of Hearing - Inserts')
%         xlabel('Frequency (kHz)')
%         ylabel('dB (SPL)')
    end

    function Threshold_dBHL = createSteeplySlopingHearingLoss_dBHL(f)
        
        Slope = 36.5; % dB/oct
        F_HearingLoss_lower = 3;
        F_HearingLoss_upper = 8;
        
        Threshold_HL = max(Slope*log2(f/F_HearingLoss_lower),0);
        Threshold_HL_upper = Threshold_HL(find(f>=F_HearingLoss_upper,1,'first'));
        
        Threshold_dBHL = min(Threshold_HL,Threshold_HL_upper);
        
%         figure
%         plot(f,Threshold_dBHL,'r--')
%         set(gca,'XLim',[min(f) max(f)])
%         title('Steeply sloping hearig loss according to screening crtiterion')
%         xlabel('Frequency (kHz)')
%         ylabel('dB (HL)')
    end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RP2-related functions
  