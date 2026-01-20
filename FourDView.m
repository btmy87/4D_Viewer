function hf = FourDView(gridVarNames,data)
% set default inputs if required
if ~exist('gridVarNames','var')
    gridVarNames = {'X1','X2','X3','X4'};
end
if ~exist('data','var')
    x1 = linspace(0,1, 16);
    x2 = linspace(0,1.2,19);
    x3 = linspace(0,1.3,21);
    x4 = linspace(0,1.4,22);
    [X1,X2,X3,X4] = ndgrid(x1,x2,x3,x4);
    data.X1 = X1;
    data.X2 = X2;
    data.X3 = X3;
    data.X4 = X4;
    data.Y1 = X1 + X2 + X3 + X4;
    data.Y2 = X1 + sin(X2) + cos(X3) + X4.^2;
    clear x1 x2 x3 x4
end

% Constants
TYP_HEIGHT = 25; % Pixels
INPUT_WIDTH = 250; % Pixels
AX_PROPS.Units = 'Pixels';
AX_PROPS.FontSize = 11;
AX_PROPS.Layer = 'top';
AX_PROPS.Box = 'on';
AX_PROPS.NextPlot = 'add';
LABEL_PROPS.Style = 'text';
LABEL_PROPS.Units = 'Pixels';
LABEL_PROPS.FontSize = 11;
LABEL_PROPS.HorizontalAlignment = 'left';
LABEL_PROPS.Position = [10,0,INPUT_WIDTH-20,TYP_HEIGHT];
EDIT_PROPS = LABEL_PROPS;
EDIT_PROPS.Style = 'edit';
LIST_PROPS = LABEL_PROPS;
LIST_PROPS.Style = 'list';
LIST_PROPS.Position(4) = 3*TYP_HEIGHT;
LIST_PROPS.Callback = @update_plots;
POPUP_PROPS = LABEL_PROPS;
POPUP_PROPS.Style = 'popupmenu';

% defaults and initialization
DEFAULT_NUM_LEVELS = 21;
levels = [];
p = [1,2,3,4]; % permutation vector.
hg1 = gobjects(1,1); % contour object in axes 1
hg2 = gobjects(1,1); % contour object in axes 2
hg3 = gobjects(1,1); % contour object in axes 3
hg4 = gobjects(3,1); % vector of slice object handles
pTextOld = '[1,2,3,4]';
yVarName = '';

% Make figure and axes
hf = figure('Name','FourDView',...
            'Visible','off',...
            'Units','Normalized',...
            'SizeChangedFcn',@resize,...
            'Position',[0.1,0.1,0.8,0.8]);
hf.Units = 'Pixels'; % set units back to pixels after initial size

ha1 = axes(AX_PROPS);
xlabel('xVar-1')
ylabel('xVar-3')

ha2 = axes(AX_PROPS);
xlabel('xVar-2')
ylabel('xVar-3')

ha3 = axes(AX_PROPS);
xlabel('xVar-1')
ylabel('xVar-2')

ha4 = axes(AX_PROPS);
xlabel('xVar-1')
ylabel('xVar-2')
zlabel('xVar-3')
view(3)

hin = uipanel('Title','Inputs',...
              'Units','pixels',...
              'FontSize',AX_PROPS.FontSize);
hlOrder1 = uicontrol(hin,'String','Var. Order');
hlOrder2 = uicontrol(hin,'String',strjoin(gridVarNames,', '));
set([hlOrder1,hlOrder2],LABEL_PROPS);
ttString = sprintf(['Enter a new permutation\n',...
                    'for the variable order.\n',...
                    'ex: [3,2,1,4]']);
heOrder = uicontrol(hin,'TooltipString',ttString,...
                        'String','[1,2,3,4]',...
                        'Callback',@update_permute);
set(heOrder,EDIT_PROPS);

ht1 = uicontrol(hin,'String',['xVar-1: ',gridVarNames{1}]);
ht2 = uicontrol(hin,'String',['xVar-2: ',gridVarNames{2}]);
ht3 = uicontrol(hin,'String',['xVar-3: ',gridVarNames{3}]);
ht4 = uicontrol(hin,'String',['xVar-4: ',gridVarNames{4}]);
hty = uicontrol(hin,'String','yVar:');
set([ht1,ht2,ht3,ht4,hty],LABEL_PROPS)

x1 = strsplit(num2str(reshape(unique(data.(gridVarNames{1})),1,[])));
x2 = strsplit(num2str(reshape(unique(data.(gridVarNames{2})),1,[])));
x3 = strsplit(num2str(reshape(unique(data.(gridVarNames{3})),1,[])));
x4 = strsplit(num2str(reshape(unique(data.(gridVarNames{4})),1,[])));
hl1 = uicontrol(hin,'String',x1);
hl2 = uicontrol(hin,'String',x2);
hl3 = uicontrol(hin,'String',x3);
hl4 = uicontrol(hin,'String',x4);
set([hl1,hl2,hl3,hl4],LIST_PROPS);
hl4.Callback = @update_clim;

varNames = fieldnames(data);
hpy = uicontrol(hin,'String',varNames,...
                    'Callback',@update_yvar);
set(hpy,POPUP_PROPS)

% build menu
hmenu = uimenu(hf,'Label','4DView');
uimenu(hmenu,'Label','Load Workspace Data',...
                     'Callback',@load_data_workspace);
uimenu(hmenu,'Label','Load Data from File ...',...
                     'Callback',@load_data_file);

% setup stuff and show figure
update_permute([],[])
update_yvar([],[])
hf.Visible = 'on';

    function resize(~,~)
        % enforce min dimensions
        if hf.Position(3) < 4*INPUT_WIDTH
            hf.Position(3) = 4*INPUT_WIDTH;
        end
        if hf.Position(4) < 2*INPUT_WIDTH
            hf.Position(4) = 2*INPUT_WIDTH;
        end
        if hf.Position(4) < 20*TYP_HEIGHT
            hf.Position(4) = 20*TYP_HEIGHT;
        end
        
        % determine basic axes dims
        W = hf.Position(3);
        H = hf.Position(4);
        wAxes = (W - INPUT_WIDTH)/2;
        hAxes = H/2;
        
        % set axes positions
        scalar = 1.5;
        ha1.Position = [    0,hAxes,wAxes,hAxes] + ...
                       scalar*(ha1.TightInset.*[1,1,-1,-1] - ...
                       [0,0,ha1.TightInset(1:2)]);
        ha2.Position = [wAxes,hAxes,wAxes,hAxes] + ...
                       scalar*(ha2.TightInset.*[1,1,-1,-1] - ...
                       [0,0,ha2.TightInset(1:2)]);
        ha3.Position = [    0,    0,wAxes,hAxes] + ...
                       scalar*(ha3.TightInset.*[1,1,-1,-1] - ...
                       [0,0,ha3.TightInset(1:2)]);
        ha4.Position = [wAxes,    0,wAxes,hAxes] + ...
                       scalar*(ha4.TightInset.*[1,1,-1,-1] - ...
                       [0,0,ha4.TightInset(1:2)]);
                   
        % set input panel properties
        hin.Position = [2*wAxes,0,INPUT_WIDTH,H];
        top = hin.Position(2)+hin.Position(4)-1*TYP_HEIGHT;
        hlOrder1.Position(2) = top -  1*TYP_HEIGHT;
        hlOrder2.Position(2) = top -  2*TYP_HEIGHT;
        heOrder.Position(2)  = top -  3*TYP_HEIGHT;
        ht1.Position(2)      = top -  5*TYP_HEIGHT;
        ht2.Position(2)      = top - 10*TYP_HEIGHT;
        ht3.Position(2)      = top - 15*TYP_HEIGHT;
        ht4.Position(2)      = top - 20*TYP_HEIGHT;
        hty.Position(2)      = top - 25*TYP_HEIGHT;
        hl1.Position(2) = ht1.Position(2) - hl1.Position(4);
        hl2.Position(2) = ht2.Position(2) - hl2.Position(4);
        hl3.Position(2) = ht3.Position(2) - hl3.Position(4);
        hl4.Position(2) = ht4.Position(2) - hl4.Position(4);
        hpy.Position(2) = hty.Position(2) - hty.Position(4);
    end

    function update_permute(hObj,~)
        % read string
        pText = heOrder.String;
        
        % split into elements, verify 4 are found
        pCellString = regexp(pText,'\d','match');
        if length(pCellString) ~= 4
            errordlg('Variable Order must be a 4 element array, ex. [1,2,3,4]')
            uiwait % wait to allow user to inspect error
            heOrder.String = pTextOld;
            return;
        end
        
        % convert to numbers, and verify that it's a permutation of 1,2,3,4
        pTemp = str2double(pCellString);
        if ~all(sort(pTemp) == [1,2,3,4])
            errordlg('Variable Order must be a permutation of [1,2,3,4]')
            uiwait % wait to allow user to inspect error
            heOrder.String = pTextOld;
            return;
        end
        
        % now we're sure we've got a valid permutation, save the outputs
        % and update all input boxes
        pOld = p;
        p = pTemp;
        pTextOld = pText;
        
        oldVal(1) = hl1.Value;
        oldVal(2) = hl2.Value;
        oldVal(3) = hl3.Value;
        oldVal(4) = hl4.Value;
        ht1.String = ['xVar-1: ',gridVarNames{p(1)}];
        ht2.String = ['xVar-2: ',gridVarNames{p(2)}];
        ht3.String = ['xVar-3: ',gridVarNames{p(3)}];
        ht4.String = ['xVar-4: ',gridVarNames{p(4)}];
        ha1.XLabel.String = ht1.String;
        ha1.YLabel.String = ht3.String;
        ha2.XLabel.String = ht2.String;
        ha2.YLabel.String = ht3.String;
        ha3.XLabel.String = ht1.String;
        ha3.YLabel.String = ht2.String;
        ha4.XLabel.String = ht1.String;
        ha4.YLabel.String = ht2.String;
        ha4.ZLabel.String = ht3.String;
        x1 = strsplit(num2str(reshape(unique(data.(gridVarNames{p(1)})),1,[])));
        x2 = strsplit(num2str(reshape(unique(data.(gridVarNames{p(2)})),1,[])));
        x3 = strsplit(num2str(reshape(unique(data.(gridVarNames{p(3)})),1,[])));
        x4 = strsplit(num2str(reshape(unique(data.(gridVarNames{p(4)})),1,[])));
        hl1.String = x1;
        hl2.String = x2;
        hl3.String = x3;
        hl4.String = x4;
        hl1.Value = oldVal(find(pOld==p(1),1));
        hl2.Value = oldVal(find(pOld==p(2),1));
        hl3.Value = oldVal(find(pOld==p(3),1));
        hl4.Value = oldVal(find(pOld==p(4),1));
        
        % finally, we can call the main update routine
        update_yvar(hObj,[]);
        update_plots(hObj,[]);
    end

    function update_plots(~,~)
        update_plot1();
        update_plot2();
        update_plot3();
        update_plot4();
    end

    function update_plot1()
        % delete contour object and replot
        hg1.delete();
        x = permute(data.(gridVarNames{p(1)}),p);
        y = permute(data.(gridVarNames{p(3)}),p);
        z = permute(data.(yVarName),p);
        idx2 = hl2.Value;
        idx4 = hl4.Value;
        x = squeeze(x(:,idx2,:,idx4));
        y = squeeze(y(:,idx2,:,idx4));
        z = squeeze(z(:,idx2,:,idx4));
        [C,hg1] = contourf(ha1,x,y,z,levels);
        clabel(C,hg1);
    end

    function update_plot2()
        % delete contour object and replot
        hg2.delete();
        x = permute(data.(gridVarNames{p(2)}),p);
        y = permute(data.(gridVarNames{p(3)}),p);
        z = permute(data.(yVarName),p);
        idx1 = hl1.Value;
        idx4 = hl4.Value;
        x = squeeze(x(idx1,:,:,idx4));
        y = squeeze(y(idx1,:,:,idx4));
        z = squeeze(z(idx1,:,:,idx4));
        [C,hg2] = contourf(ha2,x,y,z,levels);
        clabel(C,hg2);        
    end

    function update_plot3()
        % delete contour object and replot
        hg3.delete();
        x = permute(data.(gridVarNames{p(1)}),p);
        y = permute(data.(gridVarNames{p(2)}),p);
        z = permute(data.(yVarName),p);
        idx3 = hl3.Value;
        idx4 = hl4.Value;
        x = squeeze(x(:,:,idx3,idx4));
        y = squeeze(y(:,:,idx3,idx4));
        z = squeeze(z(:,:,idx3,idx4));
        [C,hg3] = contourf(ha3,x,y,z,levels);
        clabel(C,hg3);         
    end

    function update_plot4()
        hg4(:).delete();
        x = permute(data.(gridVarNames{p(1)}),p);
        y = permute(data.(gridVarNames{p(2)}),p);
        z = permute(data.(gridVarNames{p(3)}),p);
        v = permute(data.(yVarName),p);
        idx1 = hl1.Value;
        idx2 = hl2.Value;
        idx3 = hl3.Value;
        idx4 = hl4.Value;
        x = squeeze(x(:,:,:,idx4));
        y = squeeze(y(:,:,:,idx4));
        z = squeeze(z(:,:,:,idx4));
        v = squeeze(v(:,:,:,idx4));
        xi = x(idx1,   1,   1);
        yi = y(   1,idx2,   1);
        zi = z(   1,   1,idx3);
        hg4 = nd_slice(ha4,x,y,z,v,xi,yi,zi);
%         hg4 = slice(ha4,permute(y,[2,1,3]),...
%                         permute(x,[2,1,3]),...
%                         permute(z,[2,1,3]),...
%                         permute(v,[2,1,3]),yi,xi,zi);
    end

    function update_yvar(hObj,~)
        % get new variable name, update levels, and update plot
        idx = hpy.Value;
        yVarName = hpy.String{idx};
        update_clim(hObj,[])
%         update_plots(hObj,[])
    end

    function update_clim(~,~)
        y = permute(data.(yVarName),p);
        idx4 = hl4.Value;
        y = squeeze(y(:,:,:,idx4));
        yMin = min(reshape(y,1,[]));
        yMax = max(reshape(y,1,[]));
        yMax = max(yMax,yMin+10*eps(yMin)); % ensure yMax > yMin
        levels = linspace(yMin,yMax,DEFAULT_NUM_LEVELS);
        ha1.CLim = [yMin,yMax];
        ha2.CLim = [yMin,yMax];
        ha3.CLim = [yMin,yMax];
        ha4.CLim = [yMin,yMax];
        update_plots();
    end

    function load_data_workspace(~,~)
        % get list of all structures and ask user to select
        varList = evalin('base','whos');
        idxStruct = strcmp({varList(:).class},'struct');
        structList = varList(idxStruct);
        promptstring = 'Select a structure variable:';
        [selection,ok] = listdlg('PromptString',promptstring,...
                                 'ListString',{structList(:).name},...
                                 'SelectionMode','single',...
                                 'Name','Load Data - Workspace');
        if ok
            % read data from workspace and continue load process
            tempData = evalin('base',structList(selection).name);
            load_data(tempData)
        end
    end

    function load_data_file(~,~)
        
    end

    function load_data(newData)
        % get grid var names
        newVarNames = fieldnames(newData);
        newGridVarNames = cell(1,4);
        for i1 = 1:4
            promptstring = sprintf('Select grid variable #%d',i1);
            [selection,ok] = listdlg('PromptString',promptstring,...
                                     'SelectionMode','single',...
                                     'ListString',newVarNames);
            if ~ok % cancel if requested, no data has been saved
                return;
            end
            newGridVarNames{i1} = newVarNames{selection};
        end
        % verify names are unique
        if length(unique(newGridVarNames)) < length(newGridVarNames)
            errordlg('Grid variables must be unique.  Data was not loaded.')
            return;
        end
        % 
        gridVarNames = newGridVarNames;
        data = newData;
        varNames = fieldnames(data);
        hpy.String = varNames;
        hpy.Value = 1;
        hl1.Value = 1;
        hl2.Value = 1;
        hl3.Value = 1;
        hl4.Value = 1;
    end
end
