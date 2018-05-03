function scnsize = getScnSize
f = figure ('Visible','off',...
    'Menu', 'None', ...
    'Color', [0.5, 0.5, 0.5]);

set(0,'Units','pixels') ;
p = get(0, 'MonitorPositions');
H = size (p);

% switch H(1)
%     case 1
        data.scnsize = p(1,:);
%     case 2
%         m1 = p(1,:);
%         m2 = p(2,:);
%         x = m2(3) - m1(3);
%         data.scnsize = [m2(1), m2(2) - (m2(4).*0.165), x, m2(4)];
% end

position = get(f,'Position');
outerpos = get(f,'OuterPosition');
borders = outerpos - position;
data.edge = -borders(1)/2;

close (f)

 scnsize = expData.scnsize;
            edge = expData.edge;
            
            butSizeX = scnsize(3)./10;
            butSizey = scnsize(4)./20;
            txtQSizex = scnsize(3).*0.75;
            txtQSizey = scnsize(4).*0.4;
            
            
            pos.FIG = [scnsize(1),...
                scnsize(2),...
                scnsize(3) - edge,...
                scnsize(4)];
            
            pos.SLD = [scnsize(3)./2-((scnsize(3)./2)/2),...
                scnsize(4)./3,...
                scnsize(3)./2,...
                scnsize(4)./5];
            
            txtPosx = pos.SLD(1);
            txtPosy = pos.SLD(2) + (scnsize(4).*0.13);
            
            txtSizex = pos.SLD(3).*0.2;
            txtSizey = pos.SLD(4).*0.1;
            
            pos.SelectBTN = [(scnsize(3).*0.5) - (butSizeX./2),...
                (scnsize(4).*0.25) - (butSizey./2),...
                butSizeX,...
                butSizey];
            
            pos.AbortBTN = [(scnsize(3).*0.75) - (butSizeX./2),...
                (scnsize(4).*0.15) - (butSizey./2),...
                butSizeX,...
                butSizey];
            
            pos.QTXT = [(scnsize(3).*0.5) - (txtQSizex./2),...
                (scnsize(4).*0.65) - (txtQSizey./2),...
                txtQSizex,...
                txtQSizey];
            
            %     0, ‘extremely weak’; 30, ‘moderate’; 50, ‘strong’; 70, ‘very strong’; and 100, ‘extremely strong’ [15].
            pos.ExWeakTXT = [(txtPosx - (txtSizex./2))+(pos.SLD(3).*0.08),...
                (txtPosy - (txtSizey./2)),...
                txtSizex,...
                txtSizey];
            