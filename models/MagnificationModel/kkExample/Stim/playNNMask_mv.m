function playNNMask_mv(opt)
% playNNMask_mv(opt); opt = 'SIM' or 'FWD'

global HMain H1 H2 H3 ...
    partNam iptFile NTracks Hdphs...
    IPars ...
    ScrPos dGrey lGrey

ScrPos = get(0,'ScreenSize');
lGrey = 0.75*ones(1,3); dGrey = 0.375*ones(1,3); white = 0.95*[1 1 1];
FigFnt = 10; ButFnt = 16; EdtFnt = 13;

switch(upper(opt))
    case('SIM')
        IPars = struct('expNam','SIMNNMask_mv');
        titleStrg = 'Simultaneous notched-noise masking - signal vary/Quiet threshold long';
    case('FWD')
        IPars = struct('expNam','FWDNNMask_mv');
        titleStrg = 'Forward notched-noise masking - signal vary/Quiet threshold short';
    otherwise
        fprintf(1,'==> ERROR: No valid option!')
        beep, return
end
clbStrg = sprintf('play%s_clb',IPars.expNam);

% ~~~~~ Main GUI ~~~~~
HMain = figure('Position',[1 1 ScrPos(3) ScrPos(4)],...
    'Color',lGrey,...
    'CreateFcn',sprintf('%s(''Create'')',clbStrg),...
    'Menubar','none',...
    'Name',titleStrg,...
    'NumberTitle','off');

XGap = 0.0125; YGap = 0.025;
XLen = 1-2*XGap; YLen = 1-3*YGap;
xPos = zeros(1,2); yPos = zeros(1,2); xLen = zeros(1,2); yLen = zeros(1,2);
xPos(2) = XGap; yPos(2) = 2*YGap;
xLen(2) = XLen*0.45; yLen(2) = YLen*0.35;
H1.Mssg = uicontrol('Parent',HMain,...
    'Units','normalized',...
    'Position',[xPos(2) yPos(2) xLen(2) yLen(2)],...
    'BackgroundColor',0.15*[1 0 0]+0.85*[1 1 1],...
    'FontName','Arial',...
    'FontSize',FigFnt,...
    'HorizontalAlignment','left',...
    'Style','text');
xPos = lcfShift(xPos); yPos = lcfShift(yPos); xLen = lcfShift(xLen); yLen = lcfShift(yLen);
xPos(2) = 2.5*XGap; yPos(2) = yPos(1)+yLen(1)+2*YGap;
xLen(2) = xLen(1)+xPos(1)-xPos(2); yLen(2) = (YLen-yPos(2)+yPos(1)-2*YGap)/2;
H1.Fig2 = axes('Parent',HMain,...
    'Units','normalized',...
    'Position',[xPos(2) yPos(2) xLen(2) yLen(2)],...
    'FontName','Arial',...
    'FontSize',FigFnt,...
    'Box','on');
xlabel('#Decision','FontName','Arial','FontSize',FigFnt)
ylabel('MaskL (dB SPL/ERB)','FontName','Arial','FontSize',FigFnt)
xPos = lcfShift(xPos); yPos = lcfShift(yPos); xLen = lcfShift(xLen); yLen = lcfShift(yLen);
xPos(2) = xPos(1); yPos(2) = yPos(1)+yLen(1)+2*YGap;
xLen(2) = xLen(1); yLen(2) = yLen(1);
H1.Fig1 = axes('Parent',HMain,...
    'Units','normalized',...
    'Position',[xPos(2) yPos(2) xLen(2) yLen(2)],...
    'FontName','Arial',...
    'FontSize',FigFnt,...
    'Box','on');
xlabel('Frequency in kHz','FontName','Arial','FontSize',FigFnt)
ylabel('Magnitude in dB','FontName','Arial','FontSize',FigFnt)
title('Current track:')

% ~~~~~ Control buttons ~~~~~
xPos = lcfShift(xPos); yPos = lcfShift(yPos); xLen = lcfShift(xLen); yLen = lcfShift(yLen);
xPos(2) = xPos(1)+xLen(1)+XGap; yPos(2) = 2*YGap; xLen(2) = XLen-xPos(1)-xLen(1); yLen(2) = YLen*0.2;
pos = [xPos(2) yPos(2) xLen(2) yLen(2)];
uicontrol('Parent',HMain,...
    'Units','normalized',...
    'Position',[xPos(2) yPos(2) xLen(2) yLen(2)],...
    'BackgroundColor',white,...
    'Style','frame');
xPos = lcfShift(xPos); yPos = lcfShift(yPos); xLen = lcfShift(xLen); yLen = lcfShift(yLen);
xPos(2) = xPos(1)+XGap; yPos(2) = yPos(1)+0.5*YGap; xLen(2) = (xLen(1)-3.5*XGap)/4; yLen(2) = yLen(1)-YGap;
H3.StrtSess = uicontrol('Parent',HMain,...
    'Units','normalized',...
    'Position',[xPos(2) yPos(2) xLen(2) yLen(2)],...
    'BackgroundColor',dGrey,...
	'Callback',sprintf('%s(''Start'')',clbStrg),...
    'FontName','Arial',...
    'FontSize',ButFnt,...
    'Interruptible','on',...
    'Style','pushbutton',...
	'String','Start session');
xPos = lcfShift(xPos); yPos = lcfShift(yPos); xLen = lcfShift(xLen); yLen = lcfShift(yLen);
xPos(2) = xPos(1)+xLen(1)+0.5*XGap; yPos(2) = yPos(1); xLen(2) = xLen(1); yLen(2) = yLen(1);
H3.Quit = uicontrol('Parent',HMain,...
    'Units','normalized',...
    'Position',[xPos(2) yPos(2) xLen(2) yLen(2)],...
    'BackgroundColor',dGrey,...
    'BusyAction','queue',...
	'Callback',sprintf('%s(''Quit'')',clbStrg),...
    'FontName','Arial',...
    'FontSize',ButFnt,...
    'Style','pushbutton',...
	'String','Quit');
xPos = lcfShift(xPos); yPos = lcfShift(yPos); xLen = lcfShift(xLen); yLen = lcfShift(yLen);
xPos(2) = xPos(1)+xLen(1)+0.5*XGap; yPos(2) = yPos(1); xLen(2) = xLen(1); yLen(2) = yLen(1);
H3.ItrTrack = uicontrol('Parent',HMain,...
    'Units','normalized',...
    'Position',[xPos(2) yPos(2) xLen(2) yLen(2)],...
    'BackgroundColor',dGrey,...
    'Callback',sprintf('%s(''ItrTrack'')',clbStrg),...
    'Enable','off',...
    'FontName','Arial',...
    'FontSize',ButFnt,...
    'Style','pushbutton',...
	'String','Interrupt track');
xPos = lcfShift(xPos); yPos = lcfShift(yPos); xLen = lcfShift(xLen); yLen = lcfShift(yLen);
xPos(2) = xPos(1)+xLen(1)+0.5*XGap; yPos(2) = yPos(1); xLen(2) = xLen(1); yLen(2) = yLen(1);
H3.ItrSess = uicontrol('Parent',HMain,...
    'Units','normalized',...
    'Position',[xPos(2) yPos(2) xLen(2) yLen(2)],...
    'BackgroundColor',dGrey,...
    'Callback',sprintf('%s(''ItrSess'')',clbStrg),...
    'Enable','off',...
    'FontName','Arial',...
    'FontSize',ButFnt,...
    'Style','pushbutton',...
	'String','Interrupt session');

% ~~~~~ Parameters ~~~~~
xPos(2) = pos(1); yPos(2) = pos(2)+pos(4)+2*YGap; xLen(2) = pos(3); yLen(2) = YLen-yPos(2)+pos(2); 
pos = [xPos(2) yPos(2) xLen(2) yLen(2)];
uicontrol('Parent',HMain,...
    'Units','normalized',...
    'Position',[xPos(2) yPos(2) xLen(2) yLen(2)],...
    'BackgroundColor',white,...
    'Style','frame');
xLen(1) = (pos(3)-2*XGap)*0.3; xLen(2) = xLen(1)*2/3; 
N = 18; yLen(2) = (pos(4)-(2+(N-1)*0.25)*YGap)/N;
xPos(1) = pos(1)+XGap; xPos(2) = xPos(1)+xLen(1)+0.5*XGap; 
yPos(2) = YLen+YGap-yLen(2); 
uicontrol('Parent',HMain,...
    'Units','normalized',...
    'Position',[xPos(1) yPos(2) xLen(1) yLen(2)],...
    'BackgroundColor',0.75*lGrey+0.25*[1 0 0],...
    'FontName','Arial',...
    'FontSize',EdtFnt,...
    'HorizontalAlignment','left',...
    'Style','text',...
    'String','Participant''s initials:');
H2.PartNam = uicontrol('Parent',HMain,...
	'Units','normalized', ...
    'Position',[xPos(2) yPos(2) xLen(2) yLen(2)],...
    'BackgroundColor',[1 1 1],...
	'Callback',sprintf('%s(''PartNam'')',clbStrg),...
    'Enable','on',...
    'FontName','Arial',...
    'FontSize',EdtFnt,...
    'HorizontalAlignment','left',...
    'String',partNam,...
    'Style','edit');
yPos = lcfShift(yPos); yPos(2) = yPos(1)-yLen(2)-0.25*YGap; 
uicontrol('Parent',HMain,...
    'Units','normalized',...
    'Position',[xPos(1) yPos(2) xLen(1) yLen(2)],...
    'BackgroundColor',0.75*lGrey+0.25*[1 0 0],...
    'FontName','Arial',...
    'FontSize',EdtFnt,...
    'HorizontalAlignment','left',...
    'Style','text',...
    'String','Condition parameter file:');
xPos(2) = xPos(1)+xLen(1)+0.5*XGap; xLen(2) = (pos(3)-3*XGap-xLen(1))*0.7;
H2.IptFile = uicontrol('Parent',HMain, ...
	'Units','normalized', ...
    'Position',[xPos(2) yPos(2) xLen(2) yLen(2)],...
    'BackgroundColor',[1 1 1],...
	'Callback',sprintf('%s(''IptFile'')',clbStrg),...
    'Enable','on',...
    'FontName','Arial',...
    'FontSize',EdtFnt,...
    'HorizontalAlignment','left',...
    'String',iptFile,...
    'Style','edit');
xPos(2) = xPos(1)+sum(xLen)+XGap; xLen(2) = (pos(3)-3*XGap-xLen(1))*0.3;
H2.Browse = uicontrol('Parent',HMain, ...
	'Units','normalized', ...
    'Position',[xPos(2) yPos(2) xLen(2) yLen(2)],...
    'BackgroundColor',lGrey,...
	'Callback',sprintf('%s(''Browse'')',clbStrg),...
    'Enable','on',...
    'FontName','Arial',...
    'FontSize',EdtFnt,...
    'HorizontalAlignment','left',...
    'String','Browse',...
    'Style','pushbutton');
xLen(2) = pos(3)-2*XGap;
yPos = lcfShift(yPos); yPos(2) = yPos(1)-yLen(2)*1.5-0.25*YGap; 
H2.Hstry = uicontrol('Parent',HMain,...
    'Units','normalized',...
    'Position',[xPos(1) yPos(2) xLen(2) yLen(2)*1.5],...
    'BackgroundColor',lGrey,...
    'FontName','Arial',...
    'FontSize',EdtFnt,...
    'HorizontalAlignment','left',...
    'Style','text',...
    'String','history:');
xPos(2) = xPos(1)+xLen(1)+0.5*XGap; xLen(2) = pos(3)-2.5*XGap-xLen(1);
yPos = lcfShift(yPos); yPos(2) = yPos(1)-yLen(2)-0.25*YGap; 
uicontrol('Parent',HMain,...
    'Units','normalized',...
    'Position',[xPos(1) yPos(2) xLen(1) yLen(2)],...
    'BackgroundColor',lGrey,...
    'FontName','Arial',...
    'FontSize',EdtFnt,...
    'HorizontalAlignment','left',...
    'Style','text',...
    'String','target:');
H2.Trgt = uicontrol('Parent',HMain, ...
	'Units','normalized', ...
    'Position',[xPos(2) yPos(2) xLen(2) yLen(2)],...
    'BackgroundColor',[1 1 1],...
	'Callback',sprintf('%s(''Trgt'')',clbStrg),...
    'Enable','on',...
    'FontName','Arial',...
    'FontSize',EdtFnt,...
    'HorizontalAlignment','left',...
    'Style','edit');
xLen(2) = xLen(1)*2/3; 
yPos = lcfShift(yPos); yPos(2) = yPos(1)-yLen(2)-YGap; 
uicontrol('Parent',HMain,...
    'Units','normalized',...
    'Position',[xPos(1) yPos(2) xLen(1)*1.5 yLen(2)],...
    'BackgroundColor',0.75*lGrey+0.25*[1 0 0],...
    'FontName','Arial',...
    'FontSize',EdtFnt,...
    'HorizontalAlignment','left',...
    'Style','text',...
    'String','Number of tracks (NTracks):');
xPos(2) = xPos(1)+xLen(1)*1.5+0.5*XGap; 
H2.NTracks = uicontrol('Parent',HMain, ...
	'Units','normalized', ...
    'Position',[xPos(2) yPos(2) xLen(2) yLen(2)],...
    'BackgroundColor',[1 1 1],...
	'Callback','global NTracks, NTracks = eval(get(gcbo,''String''));',...
    'Enable','on',...
    'FontName','Arial',...
    'FontSize',EdtFnt,...
    'HorizontalAlignment','left',...
    'String',sprintf('%d',NTracks),...
    'Style','edit');
yPos = lcfShift(yPos); yPos(2) = yPos(1)-yLen(2)-0.25*YGap; 
uicontrol('Parent',HMain,...
    'Units','normalized',...
    'Position',[xPos(1) yPos(2) xLen(1)*1.5 yLen(2)],...
    'BackgroundColor',0.75*lGrey+0.25*[1 0 0],...
    'FontName','Arial',...
    'FontSize',EdtFnt,...
    'HorizontalAlignment','left',...
    'Style','text',...
    'String','Headphones (for NomL):');
hdphs = {'please select' 'HD 600' 'HD 280 pro'};
H2.Hdphs = uicontrol('Parent',HMain, ...
	'Units','normalized', ...
    'Position',[xPos(2) yPos(2) xLen(2) yLen(2)],...
    'BackgroundColor',[1 1 1],...
    'Enable','on',...
    'FontName','Arial',...
    'FontSize',EdtFnt,...
    'HorizontalAlignment','left',...
    'Style','popupmenu',...
    'String',hdphs,...
    'Value',Hdphs);

xPos(1) = pos(1)+XGap; xLen(2) = xLen(1)*2/3;
yPos = lcfShift(yPos); yPos(2) = yPos(1)-yLen(2)-0.25*YGap;
uicontrol('Parent',HMain,...
    'Units','normalized',...
    'Position',[xPos(1) yPos(2) xLen(1)*1.5 yLen(2)],...
    'BackgroundColor',lGrey,...
    'FontName','Arial',...
    'FontSize',EdtFnt,...
    'HorizontalAlignment','left',...
    'Style','text',...
    'String','Maximum signal level (MaxL) in dB SPL:');
H2.MaxL = uicontrol('Parent',HMain, ...
	'Units','normalized', ...
    'Position',[xPos(2) yPos(2) xLen(2) yLen(2)],...
    'BackgroundColor',[1 1 1],...
	'Callback',sprintf('%s(''MaxL'')',clbStrg),...
    'Enable','on',...
    'FontName','Arial',...
    'FontSize',EdtFnt,...
    'HorizontalAlignment','left',...
    'Style','edit',...
    'String',sprintf('%g',IPars.MaxL));
yPos = lcfShift(yPos); yPos(2) = yPos(1)-yLen(2)-0.25*YGap;
uicontrol('Parent',HMain,...
    'Units','normalized',...
    'Position',[xPos(1) yPos(2) xLen(1)*1.5 yLen(2)],...
    'BackgroundColor',lGrey,...
    'FontName','Arial',...
    'FontSize',EdtFnt,...
    'HorizontalAlignment','left',...
    'Style','text',...
    'String','Number of reversals (NRvsls):');
H2.NRvsls = uicontrol('Parent',HMain, ...
	'Units','normalized', ...
    'Position',[xPos(2) yPos(2) xLen(2) yLen(2)],...
    'BackgroundColor',[1 1 1],...
	'Callback',sprintf('%s(''NRvsls'')',clbStrg),...
    'Enable','on',...
    'FontName','Arial',...
    'FontSize',EdtFnt,...
    'HorizontalAlignment','left',...
    'Style','edit',...
    'String',sprintf('%d',IPars.NRvsls));
yPos = lcfShift(yPos); yPos(2) = yPos(1)-yLen(2)-0.25*YGap; 
uicontrol('Parent',HMain,...
    'Units','normalized',...
    'Position',[xPos(1) yPos(2) xLen(1)*1.5 yLen(2)],...
    'BackgroundColor',lGrey,...
    'FontName','Arial',...
    'FontSize',EdtFnt,...
    'HorizontalAlignment','left',...
    'Style','text',...
    'String','Number of practice trials (NPractice):');
H2.NPractice = uicontrol('Parent',HMain, ...
	'Units','normalized', ...
    'Position',[xPos(2) yPos(2) xLen(2) yLen(2)],...
    'BackgroundColor',[1 1 1],...
	'Callback','global IPars, IPars.NPractice = eval(get(gcbo,''String''));',...
    'Enable','on',...
    'FontName','Arial',...
    'FontSize',EdtFnt,...
    'HorizontalAlignment','left',...
    'Style','edit',...
    'String',sprintf('%d',IPars.NPractice));
yPos = lcfShift(yPos); yPos(2) = yPos(1)-yLen(2)-0.25*YGap; 
uicontrol('Parent',HMain,...
    'Units','normalized',...
    'Position',[xPos(1) yPos(2) xLen(1)*1.5 yLen(2)],...
    'BackgroundColor',lGrey,...
    'FontName','Arial',...
    'FontSize',EdtFnt,...
    'HorizontalAlignment','left',...
    'Style','text',...
    'String','Inter-stimulus interval (ISI) in ms:');
H2.ISI = uicontrol('Parent',HMain, ...
	'Units','normalized', ...
    'Position',[xPos(2) yPos(2) xLen(2) yLen(2)],...
    'BackgroundColor',[1 1 1],...
	'Callback','global IPars, IPars.ISI = eval(get(gcbo,''String''));',...
    'Enable','on',...
    'FontName','Arial',...
    'FontSize',EdtFnt,...
    'HorizontalAlignment','left',...
    'Style','edit',...
    'String',sprintf('%d',IPars.ISI));
yPos = lcfShift(yPos); yPos(2) = yPos(1)-yLen(2)-0.25*YGap; 
uicontrol('Parent',HMain,...
    'Units','normalized',...
    'Position',[xPos(1) yPos(2) xLen(1)*1.5 yLen(2)],...
    'BackgroundColor',lGrey,...
    'FontName','Arial',...
    'FontSize',EdtFnt,...
    'HorizontalAlignment','left',...
    'Style','text',...
    'String','Inter-trial interval (ITI) in ms:');
H2.ITI = uicontrol('Parent',HMain, ...
	'Units','normalized', ...
    'Position',[xPos(2) yPos(2) xLen(2) yLen(2)],...
    'BackgroundColor',[1 1 1],...
	'Callback','global IPars, IPars.ITI = eval(get(gcbo,''String''));',...
    'Enable','on',...
    'FontName','Arial',...
    'FontSize',EdtFnt,...
    'HorizontalAlignment','left',...
    'Style','edit',...
    'String',sprintf('%d',IPars.ITI));

xPos = fliplr(xPos); xLen(2) = pos(3)-2*XGap;
yPos = lcfShift(yPos); yPos(2) = yPos(1)-yLen(2)-YGap; 
H2.Ear = uicontrol('Parent',HMain,...
    'Units','normalized',...
    'Position',[xPos(2) yPos(2) xLen(2) yLen(2)],...
    'BackgroundColor',lGrey,...
    'FontName','Arial',...
    'FontSize',EdtFnt,...
    'HorizontalAlignment','left',...
    'Style','text',...
    'String','Ear of presentation (Ear):');
yPos = lcfShift(yPos); yPos(2) = yPos(1)-yLen(2)-0.25*YGap; 
H2.MaskG = uicontrol('Parent',HMain,...
    'Units','normalized',...
    'Position',[xPos(2) yPos(2) xLen(2) yLen(2)],...
    'BackgroundColor',lGrey,...
    'FontName','Arial',...
    'FontSize',EdtFnt,...
    'HorizontalAlignment','left',...
    'Style','text',...
    'String','Masker notch (MaskG) rel. SigF:');
yPos = lcfShift(yPos); yPos(2) = yPos(1)-yLen(2)-0.25*YGap; 
H2.SigF = uicontrol('Parent',HMain,...
    'Units','normalized',...
    'Position',[xPos(2) yPos(2) xLen(2) yLen(2)],...
    'BackgroundColor',lGrey,...
    'FontName','Arial',...
    'FontSize',EdtFnt,...
    'HorizontalAlignment','left',...
    'Style','text',...
    'String','Signal frequency (SigF) in kHz:');
yPos = lcfShift(yPos); yPos(2) = yPos(1)-yLen(2)-0.25*YGap; 
H2.SigL = uicontrol('Parent',HMain,...
    'Units','normalized',...
    'Position',[xPos(2) yPos(2) xLen(2) yLen(2)],...
    'BackgroundColor',lGrey,...
    'FontName','Arial',...
    'FontSize',EdtFnt,...
    'HorizontalAlignment','left',...
    'Style','text',...
    'String','Signal level (SigL) in dB SPL:');

% ***** lcfShift(x) *****
function y = lcfShift(x)

y = [x(2:end) 0];























