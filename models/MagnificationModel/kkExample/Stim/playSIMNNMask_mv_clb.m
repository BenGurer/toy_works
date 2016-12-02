function playSIMNNMask_mv_clb(item)

global HMain H1 H2 H3 ...
    partNam iptFile NTracks Hdphs ...
    IPars ear ...
    ItrTrack ItrSess AbrtFlag ...   
    ScrPos dGrey lGrey ...
    ePars RP2 HRsBx

switch(item)
case('Create')
    rand('state',sum(100*clock));
    randn('state',sum(100*clock));
    IPars.Gate = 10; 
    IPars.StimDur = 180; 
    IPars.MaxL = 90;
    IPars.NRvsls = 8;
    IPars.NPractice = 0;
    IPars.ISI = 250;
    IPars.ITI = 500;  
    ear = struct('strg',{'left','right'},'chan',{1 2});
    lcfIniVars
case('PartNam')
    partNam = upper(get(gcbo,'String'));
    set(gcbo,'String',partNam)
    if ~isempty(partNam)
        if ~isempty(iptFile)
            datFile = fullfile('DATA',IPars.expNam,sprintf('%s.mat',partNam));
            history = lcfHstry(datFile,IPars,ePars);       
            set(H2.Hstry,'String',char(['history:' sprintf(' %d',history(1,:))],[blanks(11) sprintf(' %d',history(2,:))]))
        end
    else
        set(H2.Hstry,'String','history:')
    end   
case('IptFile')       
    iptFile = get(gcbo,'String');
    if ~isnumeric(iptFile)&&~isempty(dir(iptFile))
        ePars = lcfLoadIptFile(iptFile);
        if isempty(ePars)
            iptFile = ''; set(gcbo,'String',iptFile);
            set(H1.Mssg,'String','==> WARNING: Condition parameter file missing or currupted!'), beep
        else
            set(H2.Trgt,'String',['[' sprintf('%d ',ones(length(ePars)-1,1)) '1]'])
            if ~isempty(partNam)
                datFile = fullfile('DATA',IPars.expNam,sprintf('%s.mat',partNam));
                history = lcfHstry(datFile,IPars,ePars);
                set(H2.Hstry,'String',char(['history:' sprintf(' %d',history(1,:))],[blanks(11) sprintf(' %d',history(2,:))]))
            end
        end
    end
case('Browse')
    iptFile = uigetfile('IptSIM*_mv.mat','Condition parameter file');
    if ~isnumeric(iptFile)&&~isempty(dir(iptFile))
        ePars = lcfLoadIptFile(iptFile);
        if isempty(ePars)
            iptFile = ''; set(H2.IptFile,'String',iptFile);
            set(H1.Mssg,'String','==> WARNING: Condition parameter file missing or currupted!'), beep
        else
            set(H2.IptFile,'String',iptFile);
            set(H2.Trgt,'String',['[' sprintf('%d ',ones(length(ePars)-1,1)) '1]'])
            if ~isempty(partNam)
                datFile = fullfile('DATA',IPars.expNam,sprintf('%s.mat',partNam));
                history = lcfHstry(datFile,IPars,ePars);
                set(H2.Hstry,'String',char(['history:' sprintf(' %d',history(1,:))],[blanks(11) sprintf(' %d',history(2,:))]))
            end
        end
    end
case('Trgt')
    if length(eval(get(gcbo,'String')))~=length(ePars)
        set(H2.Trgt,'String',['[' sprintf('%d ',ones(length(ePars)-1,1)) '1]'])
        set(H1.Mssg,'String','==> WARNING: Target has to contain %d elements!',length(ePars)), beep
    end
case('MaxL')
    IPars.MaxL = min(eval(get(gcbo,'String')),100); 
    set(gcbo,'String',sprintf('%g',IPars.MaxL))
case('NRvsls')
    IPars.NRvsls = eval(get(gcbo,'String'));
    IPars.NRvsls = max(IPars.NRvsls+mod(IPars.NRvsls,2),4);
    set(gcbo,'String',sprintf('%d',IPars.NRvsls))
case('ItrTrack')
    ItrTrack = 1;
case('ItrSess')
    ItrSess = 1;
case('Quit')
    lcfQuit(HMain,HRsBx)
case('Start')
    set(H1.Mssg,'String','')
    if isempty(partNam)
        set(H1.Mssg,'String','==> ERROR: specify participant''s name!')
        beep, return
    end
    if isempty(ePars)
        set(H1.Mssg,'String','==> ERROR: specify valid parameter file!')
        beep, return
    end
    if isempty(NTracks)
        set(H1.Mssg,'String','==> ERROR: specify NTracks!')
        beep, return
    end
    if get(H2.Hdphs,'Value')==1
        set(H1.Mssg,'String','==> ERROR: Select headphones!')
        beep, return
    end
    
    target = eval(get(H2.Trgt,'String'));
    datFile = fullfile('DATA',IPars.expNam,sprintf('%s.mat',partNam));
    history = lcfHstry(datFile,IPars,ePars);   
    if any(isnan(history))
        set(H1.Mssg,'String','==> ERROR: Parameters don''t match existing data file; remove existing file to start new one!')            
        beep, return
    end            
    if sum(max(target-history(1,:),0))==0
        set(H1.Mssg,'String','==> ERROR: Nothing to measure!')            
        beep, return
    end        
    NTracks = min(NTracks,sum(max(target-history(1,:),0)));

    Hdphs = get(H2.Hdphs,'Value');
    switch(Hdphs)
        case(2)
            NomL = 102; % Sennheiser HD 600: 1 Veff = 102 dB SPL;
        case(3)
            NomL = 113; % Sennheiser HD 280 pro: 1 Veff = 113 dB SPL;
    end
    
    if ~isdir('DATA')
        mkdir('DATA')
        mkdir(fullfile('DATA',IPars.expNam))
    else
        if ~isdir(fullfile('DATA',IPars.expNam))
            mkdir(fullfile('DATA',IPars.expNam))
        end
    end
    
    seq = [];
    for I = max(target-history(1,:)):-1:1
        jwd = find(target-history(1,:)>=I);
        jwd = jwd(randperm(length(jwd)));
        seq = [seq jwd];
    end
    seq = seq(1:NTracks);

    % ~~~~~ Response box ~~~~~
    HRsBx = figure('Position',[(0.5-1/3)*ScrPos(3) (0.5-1/6)*ScrPos(4) 2/3*ScrPos(3) 1/3*ScrPos(4)],...
        'Color',lGrey,...
        'DeleteFcn','global RP2, delete(RP2)',...
        'Menubar','none',...
        'Name','Response box',...
        'NumberTitle','off',...
        'UserData',0);
    XGap = 0.025; YGap = 0.05;
    XLen = (1-5*XGap)/4; YLen = 1-2*YGap; 
    RsBx.StrtTrack = uicontrol('Parent',HRsBx,...
        'Units','normalized',...
        'Position',[XGap YGap XLen YLen],...
        'BackgroundColor',dGrey,...
        'Callback','global HRsBx, set(HRsBx,''UserData'',4)',...
        'Enable','off',...
        'FontName','Arial',...
        'FontSize',18,...
        'Style','pushbutton',...
        'String','Start');
    RsBx.But1 = uicontrol('Parent',HRsBx,...
        'Units','normalized',...
        'Position',[2*XGap+XLen YGap XLen YLen],...
        'BackgroundColor',0.5*[1 1 0]+0.5*dGrey,...
        'Callback','global HRsBx, set(HRsBx,''UserData'',1)',...
        'Enable','off',...
        'FontName','Arial',...
        'FontSize',18,...
        'Style','pushbutton',...
        'String','Button 1');  
    RsBx.But2 = uicontrol('Parent',HRsBx,...
        'Units','normalized',...
        'Position',[3*XGap+2*XLen YGap XLen YLen],...
        'BackgroundColor',0.5*[1 0 1]+0.5*dGrey,...
        'Callback','global HRsBx, set(HRsBx,''UserData'',2)',...
        'Enable','off',...
        'FontName','Arial',...
        'FontSize',18,...
        'Style','pushbutton',...
        'String','Button 2');
    RsBx.But3 = uicontrol('Parent',HRsBx,...
        'Units','normalized',...
        'Position',[4*XGap+3*XLen YGap XLen YLen],...
        'BackgroundColor',0.5*[0 1 1]+0.5*dGrey,...
        'Callback','global HRsBx, set(HRsBx,''UserData'',3)',...
        'Enable','off',...
        'FontName','Arial',...
        'FontSize',18,...
        'Style','pushbutton',...
        'String','Button 3');

    set(H1.Mssg,'String','Creating TDT ActiveX interface...')
    RP2 = actxcontrol('RPco.x',get(HRsBx,'Position'),HRsBx);
    invoke(RP2,'ConnectRP2','USB',1);

    chain = 'NNMask';
    invoke(RP2,'ClearCOF');
    invoke(RP2,'LoadCOF',chain);
    invoke(RP2,'Run');
    Status = uint8(invoke(RP2,'GetStatus'));
    if bitget(Status,1)==0;
        set(H1.Mssg,'String',char(get(H1.Mssg,'String'),'==> ERROR: Problem connecting to RP2!'))
        beep, delete(HRsBx), HRsBx = []; return
    elseif bitget(Status,2)==0;
        set(H1.Mssg,'String',char(get(H1.Mssg,'String'),'==> ERROR: Problem loading circuit!'))
        beep, delete(HRsBx), HRsBx = []; return
    elseif bitget(Status,3)==0;
        set(H1.Mssg,'String',char(get(H1.Mssg,'String'),'==> ERROR: Problem running circuit!'))
        beep, delete(HRsBx), HRsBx = []; return
    else
        set(H1.Mssg,'String',char(get(H1.Mssg,'String'),'RP2 circuit loaded and running!'))
    end
  
    set(H2.PartNam,'Enable','off'), set(H2.IptFile,'Enable','off'), set(H2.Trgt,'Enable','off'), set(H2.NTracks,'Enable','off'), set(H2.Hdphs,'Enable','off')
    set(H2.MaxL,'Enable','off'), set(H2.NRvsls,'Enable','off'), set(H2.NPractice,'Enable','off'), set(H2.ISI,'Enable','off'), set(H2.ITI,'Enable','off')
    set(H3.StrtSess,'Enable','off'), set(H3.Quit,'Enable','off'), set(H3.ItrTrack,'Enable','on'), set(H3.ItrSess,'Enable','on')

    SF = invoke(RP2,'GetSFreq')/1000;
    HB7Gain = -(max(NomL-IPars.MaxL,0)-mod(max(NomL-IPars.MaxL,0),3)); 
    set(H1.Mssg,'String',char(get(H1.Mssg,'String'),sprintf('Set HB7Gain to %d dB! ==> Press any key to continue ...',HB7Gain))), pause
    Nfft = 2^ceil(log2(SF/2*1000));

    set(H1.Mssg,'String',['seq: ' sprintf('%d ',seq)])
    set(H1.Mssg,'String',char(get(H1.Mssg,'String'),sprintf('Data will be saved to file %s ...',datFile)))
    
    CTrack = 1;
    while CTrack<=NTracks        
        jwd = clock; IPars.Date = sprintf('%s %.2d:%.2d',date,jwd(4),jwd(5));

        Pars = ePars(seq(CTrack));
        jwd = fieldnames(IPars);
        for I = 1:length(jwd)
            Pars = setfield(Pars,jwd{I},eval(sprintf('IPars.%s',jwd{I})));
        end
        Pars = orderfields(Pars); Pars.SigL = min(Pars.SigL,Pars.MaxL);

        strg = get(H2.Ear,'String'); set(H2.Ear,'String',[strg(1:strfind(strg,':')) sprintf(' %s',ear(Pars.Ear).strg)])
        strg = get(H2.MaskG,'String'); set(H2.MaskG,'String',[strg(1:strfind(strg,':')) sprintf(' %g',Pars.MaskG)]) ;      
        strg = get(H2.SigF,'String'); set(H2.SigF,'String',[strg(1:strfind(strg,':')) sprintf(' %g',Pars.SigF)])
        strg = get(H2.SigL,'String'); set(H2.SigL,'String',[strg(1:strfind(strg,':')) sprintf(' %g',Pars.SigL)]) ;      

        strg = get(get(H1.Fig1,'Title'),'String'); set(get(H1.Fig1,'Title'),'String',[strg(1:strfind(strg,':')) sprintf(' %d',CTrack)]);    
       
        lcfTrack(CTrack,Pars,datFile,SF,NomL,HB7Gain,Nfft,RsBx)
        figure(HRsBx)
        for I = 1:3
            invoke(RP2,'SetTagVal','Feedback',15); 
            set(RsBx.But1,'BackgroundColor',[1 1 0]), drawnow 
            set(RsBx.But2,'BackgroundColor',[1 0 1]), drawnow
            set(RsBx.But3,'BackgroundColor',[0 1 1]), drawnow
            pause(0.25) 
            invoke(RP2,'SetTagVal','Feedback',0); 
            set(RsBx.But1,'BackgroundColor',0.5*[1 1 0]+0.5*dGrey), drawnow 
            set(RsBx.But2,'BackgroundColor',0.5*[1 0 1]+0.5*dGrey), drawnow
            set(RsBx.But3,'BackgroundColor',0.5*[0 1 1]+0.5*dGrey), drawnow
            pause(0.25)      
        end
        invoke(RP2,'SetTagVal','Feedback',15); 
        set(RsBx.But1,'BackgroundColor',[1 1 0]), drawnow 
        set(RsBx.But2,'BackgroundColor',[1 0 1]), drawnow
        set(RsBx.But3,'BackgroundColor',[0 1 1]), drawnow
        pause(1)
        invoke(RP2,'SetTagVal','Feedback',0); 
        set(RsBx.But1,'BackgroundColor',0.5*[1 1 0]+0.5*dGrey), drawnow 
        set(RsBx.But2,'BackgroundColor',0.5*[1 0 1]+0.5*dGrey), drawnow
        set(RsBx.But3,'BackgroundColor',0.5*[0 1 1]+0.5*dGrey), drawnow
        CTrack = CTrack+1;

        if AbrtFlag
            AbrtFlag = false;
            strg1 = get(H1.Mssg,'String'); strg2 = [deblank(strg1(end,:)) ' aborted']; set(H1.Mssg,'String',char(strg1(1:end-1,:),strg2))
            set(RsBx.But1,'Enable','off'), set(RsBx.But1,'Enable','off'), set(RsBx.But3,'Enable','off')
        elseif ItrTrack
            ItrTrack = false;
            strg1 = get(H1.Mssg,'String'); strg2 = [deblank(strg1(end,:)) ' interrupted']; set(H1.Mssg,'String',char(strg1(1:end-1,:),strg2))
            set(RsBx.But1,'Enable','off'), set(RsBx.But1,'Enable','off'), set(RsBx.But3,'Enable','off')
        elseif ItrSess
            ItrSess = false;
            strg1 = get(H1.Mssg,'String'); strg2 = [deblank(strg1(end,:)) ' interrupted']; set(H1.Mssg,'String',char(strg1(1:end-1,:),strg2,'==> Session was interrupted!'))
            set(H3.ItrTrack,'Enable','off'), set(H3.ItrSess,'Enable','off'), set(H3.StrtSess,'Enable','on'), set(H3.Quit,'Enable','on')        
            lcfIniVars, lcfIniGUI
            set(H2.PartNam,'Enable','on'), set(H2.IptFile,'Enable','on'), set(H2.Trgt,'Enable','on'), set(H2.NTracks,'Enable','on'), set(H2.Hdphs,'Enable','on')
            set(H2.MaxL,'Enable','on'), set(H2.NRvsls,'Enable','on'), set(H2.NPractice,'Enable','on'), set(H2.ISI,'Enable','on'), set(H2.ITI,'Enable','on')
            invoke(RP2,'Halt'); delete(HRsBx)
            beep, return
        else            
            strg1 = get(H1.Mssg,'String'); strg2 = [deblank(strg1(end,:)) ' completed']; set(H1.Mssg,'String',char(strg1(1:end-1,:),strg2))
        end
    end
    set(H3.ItrTrack,'Enable','off'), set(H3.ItrSess,'Enable','off'), set(H3.StrtSess,'Enable','on'), set(H3.Quit,'Enable','on')        
    lcfIniVars, lcfIniGUI
    set(H2.PartNam,'Enable','on'), set(H2.IptFile,'Enable','on'), set(H2.Trgt,'Enable','on'), set(H2.NTracks,'Enable','on'), set(H2.Hdphs,'Enable','on')
    set(H2.MaxL,'Enable','on'), set(H2.NRvsls,'Enable','on'), set(H2.NPractice,'Enable','on'), set(H2.ISI,'Enable','on'), set(H2.ITI,'Enable','on')

    set(H1.Mssg,'String',char(get(H1.Mssg,'String'),'==> Session complete!'))
    invoke(RP2,'Halt'); delete(HRsBx)
    beep
end

% ********** lcfTrack **********
function lcfTrack(CTrack,Pars,datFile,SF,NomL,HB7Gain,Nfft,RsBx)

global HMain H1 dGrey ...
    ItrTrack ItrSess AbrtFlag ...  
    RP2 ear HRsBx

if ~isempty(dir(datFile))
    load(datFile)
else
    pars = ([]);
    rvsls = {};
    rslts = struct('maskL',{},'idxC',{},'idxF',{},'idxR',{});
end
parNams = fieldnames(Pars);
for I = 1:length(parNams)
    eval(sprintf('%s = Pars.%s;',parNams{I},parNams{I}))
end

CRvsl = 0;
CTrial = 0;
Lee = lcfLee(1,SF);
SNR1dB = 10*log10(10^(1/10)-1); P = 4*SigF/lcfErb(SigF);  
MaskL = min(SigL-SNR1dB-5*log10((1+P*MaskG)*exp(-P*MaskG))-25,MaxL);
dMaskL = zeros(1,NRvsls+1); dMaskL(1) = 10; dMaskL(2) = 5; dMaskL(3:end) = 2.5;

invoke(RP2,'SetTagVal','Amp',1);
invoke(RP2,'SetTagVal','Ear',ear(Ear).chan);    
invoke(RP2,'SetTagVal','StimSize',round(StimDur*SF)+2*round(Gate*SF));    

flag = zeros(1,2);
resp = zeros(1,2); % 3A3IFC, 2-down 1-up;
Rvsls = MaskL;
Rslts.maskL = [];
Rslts.idxC = [];
Rslts.idxF = [];
Rslts.idxR = [];

Resp = 0; lcfClrResp(HRsBx,RP2);
set(H1.Mssg,'String',char(get(H1.Mssg,'String'),sprintf('Track %d: ready to start ...',CTrack))) 
set(RsBx.StrtTrack,'Enable','on'), invoke(RP2,'SetTagVal','Feedback',1);
figure(HRsBx), while Resp<4&&~or(ItrSess,ItrTrack) 
    pause(0.05), Resp = lcfChckResp(HRsBx,RP2);
end
set(RsBx.StrtTrack,'Enable','off'), invoke(RP2,'SetTagVal','Feedback',0); 
if or(ItrSess,ItrTrack)
    return
end

figure(HMain)
cla(H1.Fig1), hold(H1.Fig1,'on') 
cla(H1.Fig2), hold(H1.Fig2,'on')
set(H1.Fig1,'XLim',[0 SF/2])

I0 = 1e-12; f = SF/2*linspace(0,1,Nfft); DF = diff(f(1:2));
ee = sqrt(1./lcfErb(f));
F1(1) = max(round(SigF*(1-MaskG-0.4)/DF),1); F1(2) = round(SigF*(1-MaskG)/DF); 
F2(1) = round(SigF*(1+MaskG)/DF); F2(2) = min(round(SigF*(1+MaskG+0.4)/DF),Nfft); 
if MaskG>0
    bstp = ones(1,Nfft); bstp(1:F1(1)) = 0; bstp(F1(2):F2(1)) = 0; bstp(F2(2):end) = 0; 
else
    bstp = ones(1,Nfft); bstp(1:F1(1)) = 0; bstp(F2(2):end) = 0;
end
ntch = bstp.*ee;
mask = randn(1,2*Nfft); mask = real(ifft([ee fliplr(ee)].*fft(mask)));
mask = mask/sqrt(mean(mask.^2)); mask = real(ifft([bstp fliplr(bstp)].*fft(mask)));
sig = cos(2*pi*SigF*(0:2*Nfft-1)/SF)*sqrt(2);

F1 = round(lcfInvNErb(lcfNErb(SigF)-0.5)/DF); F2 = round(lcfInvNErb(lcfNErb(SigF)+0.5)/DF);
plot(H1.Fig1,f,10*log10(max(ntch.^2*10^(MaskL/10)*I0/sum(ee(F1:F2).^2*DF),I0)/I0),'k-')
plot(H1.Fig1,repmat(SigF,1,2),[0 SigL],'r-')

set(RsBx.But1,'Enable','on'), set(RsBx.But2,'Enable','on'), set(RsBx.But3,'Enable','on')

NC = 0; NF = 0;
while NC+NF<NPractice 
    CTrial = CTrial+1;
    invoke(RP2,'WriteTagV','Mask',0,lcfGate(mask,StimDur,Gate,SF)*10^((MaskL-Lee-NomL-HB7Gain)/20));
    invoke(RP2,'WriteTagV','SigMask',0,lcfGate(sig,StimDur,Gate,SF)*10^((SigL-NomL-HB7Gain)/20)+lcfGate(mask,StimDur,Gate,SF)*10^((MaskL-Lee-NomL-HB7Gain)/20));
    
    figure(HRsBx), Int = ceil(3*rand); % 3A3IFC;
    if Int==1          
        set(RsBx.But1,'BackgroundColor',[1 1 0]), drawnow, invoke(RP2,'SetTagVal','Feedback',2);
        T0 = clock; invoke(RP2,'SoftTrg',3); % play mask+sig;
        while etime(clock,T0)<(StimDur+2*Gate)/1000, end   
        set(RsBx.But1,'BackgroundColor',0.5*[1 1 0]+0.5*dGrey), drawnow, invoke(RP2,'SetTagVal','Feedback',0);
        T0 = clock; while etime(clock,T0)<ISI/1000, end   

        set(RsBx.But2,'BackgroundColor',[1 0 1]), drawnow, invoke(RP2,'SetTagVal','Feedback',4);
        T0 = clock; invoke(RP2,'SoftTrg',2); % play mask alone;
        while etime(clock,T0)<(StimDur+2*Gate)/1000, end   
        set(RsBx.But2,'BackgroundColor',0.5*[1 0 1]+0.5*dGrey), drawnow, invoke(RP2,'SetTagVal','Feedback',0);
        T0 = clock; while etime(clock,T0)<ISI/1000, end   

        set(RsBx.But3,'BackgroundColor',[0 1 1]), drawnow, invoke(RP2,'SetTagVal','Feedback',8);
        T0 = clock; invoke(RP2,'SoftTrg',2); % play mask alone;
        while etime(clock,T0)<(StimDur+2*Gate)/1000, end   
        set(RsBx.But3,'BackgroundColor',0.5*[0 1 1]+0.5*dGrey), drawnow, invoke(RP2,'SetTagVal','Feedback',0);
        T0 = clock; while etime(clock,T0)<ISI/1000, end   
    elseif Int==2        
        set(RsBx.But1,'BackgroundColor',[1 1 0]), drawnow, invoke(RP2,'SetTagVal','Feedback',2);
        T0 = clock; invoke(RP2,'SoftTrg',2); % play mask alone;
        while etime(clock,T0)<(StimDur+2*Gate)/1000, end   
        set(RsBx.But1,'BackgroundColor',0.5*[1 1 0]+0.5*dGrey), drawnow, invoke(RP2,'SetTagVal','Feedback',0);
        T0 = clock; while etime(clock,T0)<ISI/1000, end   

        set(RsBx.But2,'BackgroundColor',[1 0 1]), drawnow, invoke(RP2,'SetTagVal','Feedback',4);
        T0 = clock; invoke(RP2,'SoftTrg',3); % play mask+sig;
        while etime(clock,T0)<(StimDur+2*Gate)/1000, end   
        set(RsBx.But2,'BackgroundColor',0.5*[1 0 1]+0.5*dGrey), drawnow, invoke(RP2,'SetTagVal','Feedback',0);
        T0 = clock; while etime(clock,T0)<ISI/1000, end   

        set(RsBx.But3,'BackgroundColor',[0 1 1]), drawnow, invoke(RP2,'SetTagVal','Feedback',8);
        T0 = clock; invoke(RP2,'SoftTrg',2); % play mask alone;
        while etime(clock,T0)<(StimDur+2*Gate)/1000, end   
        set(RsBx.But3,'BackgroundColor',0.5*[0 1 1]+0.5*dGrey), drawnow, invoke(RP2,'SetTagVal','Feedback',0);
        T0 = clock; while etime(clock,T0)<ISI/1000, end  
    else
        set(RsBx.But1,'BackgroundColor',[1 1 0]), drawnow, invoke(RP2,'SetTagVal','Feedback',2);
        T0 = clock; invoke(RP2,'SoftTrg',2); % play mask alone;
        while etime(clock,T0)<(StimDur+2*Gate)/1000, end   
        set(RsBx.But1,'BackgroundColor',0.5*[1 1 0]+0.5*dGrey), drawnow, invoke(RP2,'SetTagVal','Feedback',0);
        T0 = clock; while etime(clock,T0)<ISI/1000, end   

        set(RsBx.But2,'BackgroundColor',[1 0 1]), drawnow, invoke(RP2,'SetTagVal','Feedback',4);
        T0 = clock; invoke(RP2,'SoftTrg',2); % play mask alone;
        while etime(clock,T0)<(StimDur+2*Gate)/1000, end   
        set(RsBx.But2,'BackgroundColor',0.5*[1 0 1]+0.5*dGrey), drawnow, invoke(RP2,'SetTagVal','Feedback',0);
        T0 = clock; while etime(clock,T0)<ISI/1000, end   

        set(RsBx.But3,'BackgroundColor',[0 1 1]), drawnow, invoke(RP2,'SetTagVal','Feedback',8);
        T0 = clock; invoke(RP2,'SoftTrg',3); % play mask+sig;
        while etime(clock,T0)<(StimDur+2*Gate)/1000, end   
        set(RsBx.But3,'BackgroundColor',0.5*[0 1 1]+0.5*dGrey), drawnow, invoke(RP2,'SetTagVal','Feedback',0);
        T0 = clock; while etime(clock,T0)<ISI/1000, end   
    end
    
    Resp = 0; lcfClrResp(HRsBx,RP2)
    invoke(RP2,'SetTagVal','Feedback',14);
    while ~and(Resp>=1,Resp<=3)&&~or(ItrSess,ItrTrack) 
        pause(0.05), Resp = lcfChckResp(HRsBx,RP2);
    end
    invoke(RP2,'SetTagVal','Feedback',0);
    if or(ItrSess,ItrTrack) 
        return
    end       
    
    T0 = clock; Rslts.maskL = [Rslts.maskL MaskL];
    if Resp==Int
        Rslts.idxC = [Rslts.idxC CTrial]; NC = NC+1;
        set(RsBx.StrtTrack,'BackgroundColor','g'), drawnow, invoke(RP2,'SetTagVal','Feedback',1);
    else
        Rslts.idxF = [Rslts.idxF CTrial]; NF = NF+1;
        set(RsBx.StrtTrack,'BackgroundColor','r'), drawnow, invoke(RP2,'SetTagVal','Feedback',8);            
    end   
    while etime(clock,T0)<ITI/1000, end   
    set(RsBx.StrtTrack,'BackgroundColor',dGrey), drawnow, invoke(RP2,'SetTagVal','Feedback',0);
end

while CRvsl<NRvsls
    CTrial = CTrial+1;
    invoke(RP2,'WriteTagV','Mask',0,lcfGate(mask,StimDur,Gate,SF)*10^((MaskL-Lee-NomL-HB7Gain)/20));
    invoke(RP2,'WriteTagV','SigMask',0,lcfGate(sig,StimDur,Gate,SF)*10^((SigL-NomL-HB7Gain)/20)+lcfGate(mask,StimDur,Gate,SF)*10^((MaskL-Lee-NomL-HB7Gain)/20));
   
    figure(HRsBx), Int = ceil(3*rand); % 3A3IFC;
    if Int==1          
        set(RsBx.But1,'BackgroundColor',[1 1 0]), drawnow, invoke(RP2,'SetTagVal','Feedback',2);
        T0 = clock; invoke(RP2,'SoftTrg',3); % play mask+sig;
        while etime(clock,T0)<(StimDur+2*Gate)/1000, end   
        set(RsBx.But1,'BackgroundColor',0.5*[1 1 0]+0.5*dGrey), drawnow, invoke(RP2,'SetTagVal','Feedback',0);
        T0 = clock; while etime(clock,T0)<ISI/1000, end   

        set(RsBx.But2,'BackgroundColor',[1 0 1]), drawnow, invoke(RP2,'SetTagVal','Feedback',4);
        T0 = clock; invoke(RP2,'SoftTrg',2); % play mask alone;
        while etime(clock,T0)<(StimDur+2*Gate)/1000, end   
        set(RsBx.But2,'BackgroundColor',0.5*[1 0 1]+0.5*dGrey), drawnow, invoke(RP2,'SetTagVal','Feedback',0);
        T0 = clock; while etime(clock,T0)<ISI/1000, end   

        set(RsBx.But3,'BackgroundColor',[0 1 1]), drawnow, invoke(RP2,'SetTagVal','Feedback',8);
        T0 = clock; invoke(RP2,'SoftTrg',2); % play mask alone;
        while etime(clock,T0)<(StimDur+2*Gate)/1000, end   
        set(RsBx.But3,'BackgroundColor',0.5*[0 1 1]+0.5*dGrey), drawnow, invoke(RP2,'SetTagVal','Feedback',0);
        T0 = clock; while etime(clock,T0)<ISI/1000, end   
    elseif Int==2        
        set(RsBx.But1,'BackgroundColor',[1 1 0]), drawnow, invoke(RP2,'SetTagVal','Feedback',2);
        T0 = clock; invoke(RP2,'SoftTrg',2); % play mask alone;
        while etime(clock,T0)<(StimDur+2*Gate)/1000, end   
        set(RsBx.But1,'BackgroundColor',0.5*[1 1 0]+0.5*dGrey), drawnow, invoke(RP2,'SetTagVal','Feedback',0);
        T0 = clock; while etime(clock,T0)<ISI/1000, end   

        set(RsBx.But2,'BackgroundColor',[1 0 1]), drawnow, invoke(RP2,'SetTagVal','Feedback',4);
        T0 = clock; invoke(RP2,'SoftTrg',3); % play mask+sig;
        while etime(clock,T0)<(StimDur+2*Gate)/1000, end   
        set(RsBx.But2,'BackgroundColor',0.5*[1 0 1]+0.5*dGrey), drawnow, invoke(RP2,'SetTagVal','Feedback',0);
        T0 = clock; while etime(clock,T0)<ISI/1000, end   

        set(RsBx.But3,'BackgroundColor',[0 1 1]), drawnow, invoke(RP2,'SetTagVal','Feedback',8);
        T0 = clock; invoke(RP2,'SoftTrg',2); % play mask alone;
        while etime(clock,T0)<(StimDur+2*Gate)/1000, end   
        set(RsBx.But3,'BackgroundColor',0.5*[0 1 1]+0.5*dGrey), drawnow, invoke(RP2,'SetTagVal','Feedback',0);
        T0 = clock; while etime(clock,T0)<ISI/1000, end  
    else
        set(RsBx.But1,'BackgroundColor',[1 1 0]), drawnow, invoke(RP2,'SetTagVal','Feedback',2);
        T0 = clock; invoke(RP2,'SoftTrg',2); % play mask alone;
        while etime(clock,T0)<(StimDur+2*Gate)/1000, end   
        set(RsBx.But1,'BackgroundColor',0.5*[1 1 0]+0.5*dGrey), drawnow, invoke(RP2,'SetTagVal','Feedback',0);
        T0 = clock; while etime(clock,T0)<ISI/1000, end   

        set(RsBx.But2,'BackgroundColor',[1 0 1]), drawnow, invoke(RP2,'SetTagVal','Feedback',4);
        T0 = clock; invoke(RP2,'SoftTrg',2); % play mask alone;
        while etime(clock,T0)<(StimDur+2*Gate)/1000, end   
        set(RsBx.But2,'BackgroundColor',0.5*[1 0 1]+0.5*dGrey), drawnow, invoke(RP2,'SetTagVal','Feedback',0);
        T0 = clock; while etime(clock,T0)<ISI/1000, end   

        set(RsBx.But3,'BackgroundColor',[0 1 1]), drawnow, invoke(RP2,'SetTagVal','Feedback',8);
        T0 = clock; invoke(RP2,'SoftTrg',3); % play mask+sig;
        while etime(clock,T0)<(StimDur+2*Gate)/1000, end   
        set(RsBx.But3,'BackgroundColor',0.5*[0 1 1]+0.5*dGrey), drawnow, invoke(RP2,'SetTagVal','Feedback',0);
        T0 = clock; while etime(clock,T0)<ISI/1000, end   
    end     

    Resp = 0; lcfClrResp(HRsBx,RP2)
    invoke(RP2,'SetTagVal','Feedback',14);     
    while ~and(Resp>=1,Resp<=3)&&~or(ItrSess,ItrTrack) 
        pause(0.05), Resp = lcfChckResp(HRsBx,RP2);
    end
    invoke(RP2,'SetTagVal','Feedback',0);
    if ~or(ItrSess,ItrTrack)     
        T0 = clock; Rslts.maskL = [Rslts.maskL MaskL];
        resp(1:end-1) = resp(2:end); flag(1:end-1) = flag(2:end);
        if Resp==Int
            Rslts.idxC = [Rslts.idxC CTrial];
            set(RsBx.StrtTrack,'BackgroundColor','g'), drawnow, invoke(RP2,'SetTagVal','Feedback',1);            
            resp(end) = 1;  
            if all(resp)
                resp(:) = 0; flag(end) = 1; 
                if abs(diff(flag))==2
                    CRvsl = CRvsl+1; Rslts.idxR = [Rslts.idxR CTrial]; Rvsls(CRvsl+1) = MaskL;
                end
                MaskL = MaskL+dMaskL(CRvsl+1);
                if MaskL>MaxL
                    AbrtFlag = true; 
                end                    
            end
        else
            Rslts.idxF = [Rslts.idxF CTrial];
            set(RsBx.StrtTrack,'BackgroundColor','r'), drawnow, invoke(RP2,'SetTagVal','Feedback',8);            
            resp(end) = 0; flag(end) = -1;
            if abs(diff(flag))==2
                CRvsl = CRvsl+1; Rslts.idxR = [Rslts.idxR CTrial]; Rvsls(CRvsl+1) = MaskL;
            end
            MaskL = MaskL-dMaskL(CRvsl+1);
        end
        figure(HMain)
        plot(H1.Fig2,1:CTrial,Rslts.maskL,'k-','LineWidth',1.3)
        plot(H1.Fig2,Rslts.idxC,Rslts.maskL(Rslts.idxC),'kv','MarkerSize',6)
        plot(H1.Fig2,Rslts.idxF,Rslts.maskL(Rslts.idxF),'k^','MarkerSize',6)
        plot(H1.Fig2,Rslts.idxR,Rslts.maskL(Rslts.idxR),'o','MarkerSize',8,'MarkerFaceColor','none')

        while etime(clock,T0)<ITI/1000, end   
        set(RsBx.StrtTrack,'BackgroundColor',dGrey), drawnow, invoke(RP2,'SetTagVal','Feedback',0);            
        if AbrtFlag
            lcfSave(datFile,pars,Pars,rslts,Rslts,rvsls,Rvsls)
            return
        end
    else  
        return
    end
end
lcfSave(datFile,pars,Pars,rslts,Rslts,rvsls,Rvsls)                   
set(RsBx.But1,'Enable','off'), set(RsBx.But2,'Enable','off'), set(RsBx.But3,'Enable','off')

% ***** lcfIniVars *****
function lcfIniVars
global partNam iptFile ePars NTracks Hdphs ...
    ItrTrack ItrSess AbrtFlag

partNam = ''; iptFile = ''; ePars = ([]); NTracks = []; Hdphs = 1;
ItrTrack = false; ItrSess = false; AbrtFlag = false;

% ***** lcfIniGUI *****
function lcfIniGUI
global partNam iptFile NTracks Hdphs H2

set(H2.PartNam,'String',partNam), set(H2.IptFile,'String',iptFile), 
set(H2.NTracks,'String',sprintf('%d',NTracks)),set(H2.Hdphs,'Value',Hdphs)
set(H2.Hstry,'String','history:'), set(H2.Trgt,'String','')

% ***** lcfLoadIptFile *****
function ePars = lcfLoadIptFile(iptFile) 

load(iptFile), ePars = orderfields(ePars);
parNams = sort({'Ear';'MaskG';'SigF';'SigL'});
if length(fieldnames(ePars))~=length(parNams)||~all(strcmp(fieldnames(ePars),parNams))
    ePars = ([]);
end

% ***** lcfQuit *****    
function lcfQuit(HMain,HRsBx)

delete(HMain), if isobject(HRsBx), delete(HRsBx), end
fprintf(1,'\nGoodbye ...\n\n\n');

% ***** lcfHstry *****
function history = lcfHstry(datFile,IPars,ePars)

NConds = length(ePars); history = zeros(2,length(ePars));
if ~isempty(dir(datFile))
    load(datFile), pars = rmfield(pars,'Date'); ePars = orderfields(ePars);
    
    Pars = ePars(1);
    jwd = fieldnames(IPars);
    for I = 1:length(jwd)
        Pars = setfield(Pars,jwd{I},eval(sprintf('IPars.%s',jwd{I})));
    end
    Pars = orderfields(Pars); if isfield(Pars,'Date'), Pars = rmfield(Pars,'Date'); end
    if length(fieldnames(Pars))~=length(fieldnames(pars))||~all(strcmp(fieldnames(Pars),fieldnames(pars)))
        history = nan(2,length(ePars)); 
    else            
        history = zeros(2,length(ePars));
        allNams = fieldnames(pars);
        eNams = fieldnames(ePars);
        idx = true(size(allNams));
        for I = 1:length(eNams)
            idx = and(~strcmp(eNams{I},allNams),idx);
        end
        rPars = rmfield(pars,allNams(idx));
        for I = 1:NConds
            idx = [];
            for II = 1:length(rPars)            
                if isequalwithequalnans(rPars(II),ePars(I))
                    idx = [idx II];
                end
            end
            if ~isempty(idx)
                history(1,I) = length(find(cellfun(@(x,y)length(x)==y+1,rvsls(idx),num2cell(arrayfun(@(x)x.NRvsls,pars(idx))))));
                history(2,I) = length(find(cellfun(@(x,y)length(x)<y+1,rvsls(idx),num2cell(arrayfun(@(x)x.NRvsls,pars(idx))))));
            end
        end
    end             
end  

% ***** lcfClrResp ***** 
function lcfClrResp(HRsBx,RP2)

invoke(RP2,'SoftTrg',1);
set(HRsBx,'UserData',0)

% ***** lcfChckResp ***** 
function Resp = lcfChckResp(HRsBx,RP2)

But= double(invoke(RP2,'GetTagVal','Button'));
if But==1
    Resp = 4; invoke(RP2,'SoftTrg',1);
elseif But==2
    Resp = 1; invoke(RP2,'SoftTrg',1);
elseif But==4
    Resp = 2; invoke(RP2,'SoftTrg',1);
elseif But==8
    Resp = 3; invoke(RP2,'SoftTrg',1);
else
    Resp = get(HRsBx,'UserData');
end

% ***** lcfSave *****
function lcfSave(datFile,pars,Pars,rslts,Rslts,rvsls,Rvsls)

pars = [pars;Pars];
rslts = [rslts;Rslts];
rvsls = [rvsls;{Rvsls}];   
save(datFile,'pars','rslts','rvsls')

parNams = fieldnames(Pars);
[pth,nam] = fileparts(datFile); FId = fopen(fullfile(pth,[nam '.dat']),'at');
fprintf(FId,'Pars = struct(');
for I = 1:length(parNams)-1
    if ischar(eval(sprintf('Pars.%s',parNams{I})))
        fprintf(FId,'''%s'',%s,',parNams{I},eval(sprintf('Pars.%s',parNams{I})));
    else
        if numel(eval(sprintf('Pars.%s',parNams{I})))==1
            fprintf(FId,'''%s'',%g,',parNams{I},eval(sprintf('Pars.%s',parNams{I})));
        else
            fprintf(FId,'''%s'',[',parNams{I});fprintf(FId,' %g',eval(sprintf('Pars.%s',parNams{I}))); fprintf(FId,'],');
        end
    end
end
if ischar(eval(sprintf('Pars.%s',parNams{end})))
    fprintf(FId,'''%s'',%s);\n',parNams{end},eval(sprintf('Pars.%s',parNams{end})));
else
    if numel(eval(sprintf('Pars.%s',parNams{end})))==1
        fprintf(FId,'''%s'',%g);\n',parNams{end},eval(sprintf('Pars.%s',parNams{end})));
    else
        fprintf(FId,'''%s'',[',parNams{end});fprintf(FId,' %g',eval(sprintf('Pars.%s',parNams{end}))); fprintf(FId,']);\n');                                
    end
end

fprintf(FId,'Rslts = struct(''MaskL'',[');fprintf(FId,' %g',Rslts.maskL);fprintf(FId,'],...\n');
fprintf(FId,'\t''idxC'',[');fprintf(FId,' %g',Rslts.idxC);fprintf(FId,'],...\n');
fprintf(FId,'\t''idxF'',[');fprintf(FId,' %g',Rslts.idxF);fprintf(FId,'],...\n');
fprintf(FId,'\t''idxR'',[');fprintf(FId,' %g',Rslts.idxR);fprintf(FId,'];\n');

fprintf(FId,'Rvsls = ['); fprintf(FId,' %g',Rvsls); fprintf(FId,'];\n\n');
fclose(FId);

% ***** lcfLee ***** 
function Lee = lcfLee(F,SF)
% Lee = level of ee-noise within 1 ERB around F of wideband ee-noise relative to overall level;

F1 = lcfInvNErb(lcfNErb(F)-0.5); F2 = lcfInvNErb(lcfNErb(F)+0.5);
Lee = 10*log10(lcfIntErb(F1,F2)/lcfIntErb(0,SF/2));

% ***** lcfIntErb *****
function I = lcfIntErb(F1,F2)

A = 24.7/1000; B = 4.37;
I = (log(A*(B*F2+1))-log(A*(B*F1+1)))/(A*B);

% ***** lcfErb *****
function erb = lcfErb(f)

A = 24.7/1000; B = 4.37;
erb = A*(B*f+1);

% ***** lcfNErb *****
function nerb = lcfNErb(f)

A = 24.7/1000; B = 4.37;
nerb = 1/(A*B)*log(B*f+1);

% ***** lcfInvNErb *****
function f = lcfInvNErb(nerb)

A = 24.7/1000; B = 4.37;
f = 1/B*(exp(A*B*nerb)-1);

% ***** lcfGate *****
function stim = lcfGate(stim,StimDur,Gate,SF)

Gate = round(Gate*SF); StimSize = round(StimDur*SF); env = cos(pi/2*(0:Gate-1)/(Gate-1));

stim = circshift(stim,ceil(length(stim)*rand)); stim = stim(1:StimSize+2*Gate);
stim(1:Gate) = stim(1:Gate).*sqrt(1-env.^2); stim(end-Gate+1:end) = stim(end-Gate+1:end).*env;







