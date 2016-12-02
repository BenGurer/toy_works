function data = wrtIptNNMask_mv(partNam,opt)
% wrtIptNNMask_mv(partNam); partNam = e.g., 'OT'; opt = 'SIM' or 'FWD'

partNam = upper(partNam);
switch(upper(opt))
    case('SIM')
        expNam = 'SIMNNMask_sv';
    case('FWD')
        expNam = 'FWDNNMask_sv';
    otherwise
        fprintf(1,'==> ERROR: No valid option!')
        beep, return
end

iptFile = 'IptQuietThr_sv.mat';
ePars = lcfLoadIptFile(iptFile);
if isempty(ePars)
    fprintf(1,'==> ERROR: Condition parameter file missing or currupted!\n');
    beep, return
end
NConds = length(ePars);
history = zeros(NConds,1);

datFile = fullfile('DATA',expNam,sprintf('%s.mat',partNam));
if ~isempty(dir(datFile))
    load(datFile), pars = orderfields(pars); ePars = orderfields(ePars);
    allNams = fieldnames(pars);
    eNams = fieldnames(ePars);
    idx = true(size(allNams));
    for I = 1:length(eNams)
        idx = and(~strcmp(eNams{I},allNams),idx);
    end
    rPars = rmfield(pars,allNams(idx));
    
    data = cell(NConds+1,length(fieldnames(ePars))+2);
    data(1,:) = [eNams' {'me'} {'se'}];
    for I = 1:NConds
        for II = 1:length(eNams)
            data(I+1,II) = {eval(sprintf('ePars(I).%s',eNams{II}))};
        end
        
        idx = [];
        for II = 1:length(rPars)            
            if isequalwithequalnans(rPars(II),ePars(I))
                idx = [idx II];
            end
        end
        if ~isempty(idx)
            idx = idx(find(cellfun(@(x,y)length(x)==y+1,rvsls(idx),num2cell(arrayfun(@(x)x.NRvsls,pars(idx))))));
        end
        history(I) = length(idx);
        
        if ~isempty(idx)
            jwd = [];
            for II = 1:length(idx)
                jwd = [jwd mean(reshape(rvsls{idx(II)}(4:end),[2 length(rvsls{idx(II)}(4:end))/2]))];
            end

            data(I+1,length(fieldnames(ePars))+1) = {mean(jwd)};
            data(I+1,length(fieldnames(ePars))+2) = {std(jwd)/sqrt(length(jwd))};
        else
            data(I+1,length(fieldnames(ePars))+1) = nan;
            data(I+1,length(fieldnames(ePars))+2) = nan;
        end    
    end
    
    if ~any(isnan(cellfun(@(x)x,data(2:end,length(fieldnames(ePars))+1))))
        if any(cellfun(@(x)x>5,data(2:end,length(fieldnames(ePars))+2)))
            fprintf(1,'==> WARNING: High variance in condition %d! Consider remeasuring ...',find(data(:,length(fieldnames(ePars))+2)>5))
        end

        fieldNams = {'Ear','MaskG','SigF','SigL'};
        SNR1dB = 10*log10(10^(1/10)-1);
        maskG = [0 0.1 0.2 0.3]; N1 = length(maskG); N2 = size(data,1)-1;
        Idx = find(strcmp('Ear',data(1,:)));
        jwd = repmat(data(2:end,Idx)',1,N1); % Ear (1 = left, 2 = right);       
        jwd = [jwd;eval(['[' sprintf('repmat({maskG(%d)},1,N2) ',1:length(maskG)) ']'])]; % MaskG;
        Idx = find(strcmp('SigF',data(1,:)));
        jwd = [jwd;repmat(data(2:end,Idx)',1,N1)]; % SigF;
        jwd = [jwd;repmat(num2cell(cellfun(@(x)x+10,data(2:end,length(fieldnames(ePars))+1)))',1,N1)]; % SigL (10 dB SL);
        jwd           
        ePars = orderfields(cell2struct(jwd,fieldNams,1)); 

        iptNam = sprintf('Ipt%sNNMask_%s_mv',upper(opt),partNam);
        save([iptNam '.mat'],'ePars')

        FId = fopen([iptNam '.txt'],'wt');
        for I = 1:length(ePars)
            fprintf(FId,'ePars(%d) = struct(',I);
            for II = 1:length(fieldNams)-1
                if ischar(eval(sprintf('ePars(%d).%s',I,fieldNams{II})))
                    fprintf(FId,'''%s'',''%s'',',fieldNams{II},eval(sprintf('ePars(%d).%s',I,fieldNams{II})));
                else
                    if length(eval(sprintf('ePars(%d).%s',I,fieldNams{II})))==1
                        fprintf(FId,'''%s'',%g,',fieldNams{II},eval(sprintf('ePars(%d).%s',I,fieldNams{II})));
                    else
                        fprintf(FId,'''%s'',[',fieldNams{II}); fprintf(FId,' %g',eval(sprintf('ePars(%d).%s',I,fieldNams{II}))); fprintf(FId,'],');
                    end
                end
            end
            if ischar(eval(sprintf('ePars(%d).%s',I,fieldNams{end})))
                fprintf(FId,'''%s'',''%s'');\n',fieldNams{end},eval(sprintf('ePars(%d).%s',I,fieldNams{end})));
            else
                if length(eval(sprintf('ePars(%d).%s',I,fieldNams{end})))==1
                    fprintf(FId,'''%s'',%g);\n',fieldNams{end},eval(sprintf('ePars(%d).%s',I,fieldNams{end})));
                else
                    fprintf(FId,'''%s'',[',fieldNams{end}); fprintf(FId,' %g',eval(sprintf('ePars(%d).%s',I,fieldNams{end}))); fprintf(FId,']);\n\n');
                end
            end
        end
        fclose(FId);
    else
        fprintf('==> ERROR: Data set incomplete!');
        beep, return
    end
else
    fprintf('==> ERROR: No data found!');
    data = [];
end

% ***** lcfLoadIptFile *****
function ePars = lcfLoadIptFile(iptFile) 

parNams = sort({'Ear';'MaskG';'MaskL';'SigF'});
if ~isempty(dir(iptFile))
    load(iptFile)
    if length(fieldnames(ePars))~=length(parNams)||~all(strcmp(sort(fieldnames(ePars)),parNams))
        ePars = ([]);
    end
else
    ePars = ([]);
end





