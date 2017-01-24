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
% Acoustics - Reference zero for the calibration of audiometric equipment — Part 2: Reference equivalent threshold sound 
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

% f = linspace(20,20000,20);

figure
plot(f_hz_measured_125_8000Hz,RETSPL_125_8000Hz,'g')
hold on
plot(f_hz_measured_8000_16000Hz,RETSPL_8000_1600Hz,'b')

RETSPL_int = interp1([f_hz_measured_125_8000Hz f_hz_measured_8000_16000Hz(2:end)],[RETSPL_125_8000Hz RETSPL_8000_1600Hz(2:end)],f,'spline');
plot(f,RETSPL_int,'r--')
