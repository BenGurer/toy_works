
% Filter
function EEG = kkFilt(EEG,SF,LPF,HPF)
Order = 8; % ==> produces order=16 (48 dB/oct filter roll-off at -3-dB down) for second pass of filter
Fc1 = HPF;
Fc2 = LPF;
% Calculate the zpk values using the BUTTER function.
if isempty(HPF)
    fprintf(1,'Filtering (Order = %d; SF = %g kHz; LPF = %g kHz; ...',Order,SF,LPF)
    [z,p,k] = butter(Order,Fc2/(SF/2),'low');
elseif isempty(LPF)
    fprintf(1,'Filtering (Order = %d; SF = %g kHz; HPF = %g kHz; ...',Order,SF,HPF)
    [z,p,k] = butter(Order,Fc1/(SF/2),'high');
else
    fprintf(1,'Filtering (Order = %d; SF = %g kHz; LPF = %g kHz; HPF = %g kHz; ...',Order,SF,LPF,HPF)
    [z,p,k] = butter(Order,[Fc1 Fc2]/(SF/2));
end  
% To avoid round-off errors, do not use the transfer function. Instead
% get the zpk representation and convert it to second-order sections.
[sos_var,g] = zp2sos(z,p,k);
Hd = dfilt.df2sos(sos_var,g);
NChans = size(EEG.data,1);
for I = 1:NChans
    EEG.data(I,:) = single(fliplr(filter(Hd,fliplr(filter(Hd,EEG.data(I,:))))));  
end
fprintf(1,' %d channels filtered!\n',I)
