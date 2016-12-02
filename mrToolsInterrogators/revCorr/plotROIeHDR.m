
function plotROIeHDR(filename)
% load n files, loop and plot all
load(filename)
f = figure;
numberEVs = length(e.contrastBetas);
% select the window to plot into
fignum = get(f,'Number');
initializeFigure(fignum,numberEVs)
set(fignum,'Name',['glmPlot: ' filename]);
ehdrAxes = axes('parent',fignum);
for iSte = 1:size(e.betaSte,3)
    [h,hEhdrSte]=plotEhdr(ehdrAxes,e.time,e.hdr,e.hdrSte(:,:,iSte),[],[],iSte~=1);
    % hold on
end
ylabel(ehdrAxes,{'Scaled HRF','% Signal change'});
lhandle = legend(hHdr,EVnames,'position',legendEhdrPosition);
set(lhandle,'Interpreter','none','box','off');
xlabel(ehdrAxes,'Time (sec)');
figure
pehdr=permute(e.hdr,[2,1]);
waterfall(pehdr);

function initializeFigure(fignum,numberColors)

lineWidth = 2;
fontSize = 15;

%set default values for figure aspect
set(fignum,'DefaultLineLineWidth',lineWidth);
set(fignum,'DefaultAxesFontSize',fontSize);
%set the colors
colors = color2RGB;
colors = colors([7 5 6 8 4 3 2 1]); %remove white and black and re-order
for i_color = 1:length(colors)
    colorOrder(i_color,:) = color2RGB(colors{i_color});
end
if numberColors>size(colorOrder,1)
    colorOrder(end+1:numberColors,:) = randomColors(numberColors-size(colorOrder,1));
end
colorOrder = colorOrder(1:numberColors,:);


set(fignum,'DefaultAxesColorOrder',colorOrder);
%for bars, need to set the colormap
set(fignum,'colormap',colorOrder);

% turn off menu/title etc.
set(fignum,'NumberTitle','off');
function [h,hSte] = plotEhdr(hAxes,time,ehdr,ehdrSte,lineSymbol,drawSymbols, steOnly)

colorOrder = get(hAxes,'colorOrder');
% whether to plot the line inbetween points or not
if ieNotDefined('lineSymbol'),lineSymbol = '-';end
if ieNotDefined('drawSymbols'),drawSymbols = 1;end
if ieNotDefined('steOnly'),steOnly = 0;end

% and display ehdr
if steOnly
    h=[];
else
    %      figure;
    %      pehdr=permute(ehdr,[2,1]);
    %      [nrow ncol] = size(pehdr);
    %      imagesc(pehdr);
    
    h=plot(hAxes,time,ehdr,lineSymbol);
    
    if drawSymbols
        for iEv = 1:size(ehdr,2)
            set(h(iEv),'Marker',getsymbol(iEv),'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',colorOrder(iEv,:));
        end
    end
end

if ~ieNotDefined('ehdrSte')
    hold on
    %if ~ishold(hAxes),hold(hAxes);end;
    hSte=errorbar(hAxes,repmat(time',1,size(ehdr,2)),ehdr,ehdrSte,ehdrSte,'lineStyle','none','LineWidth',1)';
end
