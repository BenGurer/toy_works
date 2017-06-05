function plotStimuli
stimTR = 7500;
TR = 7500;
[params,stimulus] = corticalMagnification([],[],stimTR,TR);
scanner.frequency = 1;
scanner.bandwidth = 100;
scanner.duration = 1250;

% stimulus(1).frequency
% stimulus(1).duration
% stimulus(1).number

allFrequencies = funInvNErb(linspace(funNErb(params.lowFrequency),funNErb(params.highFrequency),params.nFrequencies));
lowCuttingFrequencies = funInvNErb(funNErb(allFrequencies)-params.bandwidthERB/2);
highCuttingFrequencies = funInvNErb(funNErb(allFrequencies)+params.bandwidthERB/2);
allFrequencies = (lowCuttingFrequencies+highCuttingFrequencies)/2;
allBandwidths = (highCuttingFrequencies-lowCuttingFrequencies);
lowCuttingFrequencies = round(lowCuttingFrequencies*1000);
highCuttingFrequencies = round(highCuttingFrequencies*1000);

% f = 50:10:10000;
f = logspace(log10(50),log10(10000),300);
silence = -10;
% stimuli = nan(length(f),length(allFrequencies));
stimuli = repmat(silence,length(f),length(allFrequencies));
for i = 1:length(allFrequencies)
    %     indexLower = f(f>=lowCuttingFrequencies(i) & f<=highCuttingFrequencies(i));
    index= f>=lowCuttingFrequencies(i) & f<=highCuttingFrequencies(i);
    stimuli(index,i) = i;
    
    stimuliFrequencies{i} = f(index);
end

c = 0;
for i = 1:6
    id = stimulus(i).number;
    %     c = c + scanner.duration;
    c_is = 0;
    
%     % scanner noise
%     scannerTimeOn{i} = repmat(c,1,length(f));
%     scannerIm{i} = ones(length(f),1);
%     scannerDur{i} = repmat(1.25,1,length(f));
%     scanNoise{i} = repmat(1,1,length(f));
%     scanCon{i} = repmat(0,1,length(f));
    
    % stimulus
    for ii = 1:length(stimulus(i).duration)
        if isnan(stimulus(i).frequency(ii))
            
%             if ii == 1
%                 stimIm{i}{ii} = ones(length(f),1);
%                 stimCon{i}{ii} = repmat(id,1,length(f));
%                 
%                 stimTimeOn{i}{ii}= repmat(c+c_is,1,length(f));
%                 stimTime{i}{ii} = repmat(stimulus(i).duration(ii)/1000+c+c_is,1,length(f));
%                 stimDur{i}{ii} = repmat(stimulus(i).duration(ii),1,length(f));
%                 stimDurStr{i}{ii} = repmat(mat2str(stimulus(i).duration(ii)),1,length(f));
%                 c_is = c_is + stimulus(i).duration(ii)/1000;
%                 stimNoise{i}{ii} = repmat(1,1,length(f));
%             else
                stimIm{i}{ii} = 0;
                stimNoise{i}{ii} = 1;
                stimDur{i}{ii} = 0;
                stimDurStr{i}{ii} = '0';
                stimCon{i}{ii} = id;
                stimTimeOn{i}{ii} = c +c_is;
                stimTime{i}{ii} = c +c_is + stimulus(i).duration(ii)/1000;
                c_is = c_is + stimulus(i).duration(ii)/1000;
%             end
            
        else
            stimIm{i}{ii} = stimuliFrequencies{id};
            stimCon{i}{ii} = repmat(id,1,length(stimuliFrequencies{id}));
            
            stimTimeOn{i}{ii}= repmat(c+c_is,1,length(stimuliFrequencies{id}));
            stimTime{i}{ii} = repmat(stimulus(i).duration(ii)/1000+c+c_is,1,length(stimuliFrequencies{id}));
            stimDur{i}{ii} = repmat(stimulus(i).duration(ii),1,length(stimuliFrequencies{id}));
            stimDurStr{i}{ii} = repmat(mat2str(stimulus(i).duration(ii)),1,length(stimuliFrequencies{id}));
            c_is = c_is + stimulus(i).duration(ii)/1000;
            stimNoise{i}{ii} = repmat(1,1,length(stimuliFrequencies{id}));
        end
    end
    c = c +(7500/1000);
    
end
unpack.stimTime = [];
unpack.stimTimeOn = [];
unpack.stimIm = [];
unpack.stimDur = [];
unpack.stimCon = [];
unpack.scanNoise = [];


for i = 1:6
    %     unpack.stimTime = [unpack.stimTime stimTime{i}];
    unpack.stimTimeOn = [unpack.stimTimeOn stimTimeOn{i}];
    unpack.stimIm = [unpack.stimIm stimIm{i}];
    unpack.stimDur = [unpack.stimDur stimDur{i}];
    unpack.stimCon = [unpack.stimCon stimCon{i}];
    unpack.scanNoise = [unpack.scanNoise stimNoise{i}];
end




figure
g = gramm('x',unpack.stimTimeOn,'y',unpack.stimIm,'color',unpack.stimCon);
% ,'color',unpack.stimCon,'lightness',unpack.stimDur);
g.geom_point()
% Set appropriate names for legends
g.set_names('x','Time (Seconds)','y','Frequency (Hertz)','color','Stimuli ID')
%Set figure title
g.set_title('Stimulus Presentation')
g.axe_property('YLim',[50 10000]);
g.axe_property('YScale','log');

g.axe_property('XLim',[0 40]);
g.axe_property('CLim',[1 32]);
g.draw()

figure
b = gramm('x',stimTimeOn{1},'y',stimIm{1},'color',stimDur{i});
b.geom_point()
b.axe_property('YLim',[3500 4000]);
b.axe_property('YScale','log');
b.axe_property('XLim',[0 7.5]);
b.draw()

