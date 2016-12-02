function wrtIptQuietThr_sv

sigF = {1 4}; N = length(sigF);

jwd = repmat({1},1,N); % Ear (1 = left, 2 = right, 3 = binaural);
jwd = [jwd;repmat({nan},1,N)]; % MaskG;
jwd = [jwd;repmat({0},1,N)]; % MaskL;
jwd = [jwd;sigF]; % SigF;
jwd

fieldNams = {'Ear','MaskG','MaskL','SigF'};
ePars = orderfields(cell2struct(jwd,fieldNams,1)); 

iptNam = 'IptQuietThr_sv';
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

