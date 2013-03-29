function varargout = ClaudiaGUI_v2(varargin)
%CLAUDIAGUI_V2 M-file for ClaudiaGUI_v2.fig
%      CLAUDIAGUI_V2, by itself, creates a new CLAUDIAGUI_V2 or raises the existing
%      singleton*.
%
%      H = CLAUDIAGUI_V2 returns the handle to a new CLAUDIAGUI_V2 or the handle to
%      the existing singleton*.
%
%      CLAUDIAGUI_V2('Property','Value',...) creates a new CLAUDIAGUI_V2 using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to ClaudiaGUI_v2_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      CLAUDIAGUI_V2('CALLBACK') and CLAUDIAGUI_V2('CALLBACK',hObject,...) call the
%      local function named CALLBACK in CLAUDIAGUI_V2.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ClaudiaGUI_v2

% Last Modified by GUIDE v2.5 18-Aug-2011 13:52:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ClaudiaGUI_v2_OpeningFcn, ...
                   'gui_OutputFcn',  @ClaudiaGUI_v2_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ClaudiaGUI_v2 is made visible.
function ClaudiaGUI_v2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for ClaudiaGUI_v2
handles.output = hObject;
handles.x = load('ClaudiaGUI_v2.mat');

% Update handles structure
guidata(hObject, handles);



% UIWAIT makes ClaudiaGUI_v2 wait for user response (see UIRESUME)
% uiwait(gca1);

% --- Outputs from this function are returned to the command line.
function varargout = ClaudiaGUI_v2_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = handles.output;




function gene_Callback(hObject, eventdata, handles)
% hObject    handle to gene (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gene as text
%        str2double(get(hObject,'String')) returns contents of gene as a double

gene = get(hObject,'String');
if ischar(gene)
    set(handles.pushbutton5,'Enable','on')
else
    set(handles.pushbutton5,'Enable','off')
end

% --- Executes on button press in pushbutton1.
function Rad23b_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton2.
function Oct4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in pushbutton3.
function Sox2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Initializing Steps
c1 = get(handles.Rad23b,'Value');
c2 = get(handles.Oct4,'Value');
c3 = get(handles.Sox2,'Value');
c4 = get(handles.Nanog,'Value');
c5 = get(handles.PI,'Value');
c6 = get(handles.Rad23b_Repeat,'Value');
c7 = get(handles.PI_Repeat,'Value');
headers = {get(handles.Rad23b,'String') get(handles.Oct4,'String') get(handles.Sox2,'String') get(handles.Nanog,'String') ...
    get(handles.Rad23b_Repeat,'String') get(handles.PI,'String') get(handles.PI_Repeat,'String')};

Yaxis = get(handles.RawReads,'Value');


A = [c1 c2 c3 c4 c6 c5 c7];
f = sum(A);
k = find(A);
headers = headers(k);
ReadSum = handles.x.ReadSum(k);
if f == 0
    errordlg('Select an IP')
end

axisoption = get(handles.popupmenu1,'Value');


GG = get(handles.checkgene,'Value');
chrcoor = get(handles.chrsearch,'Value');
if GG && chrcoor
    errordlg('Select one of the search options')
end
if GG == 1
    gene = get(handles.gene,'String');
    g = find(strcmp(gene,handles.x.knownGene(:,2)),1,'first');
    if isempty(g)
        errordlg('Enter an appropriate gene name')
    end
    chr = handles.x.knownGene{g,4};
    x1 = cell2mat(handles.x.knownGene(g,6));
    x2 = cell2mat(handles.x.knownGene(g,7));
    set(handles.chr,'String',handles.x.knownGene{g,4})
    set(handles.rightedge,'String',num2str(x2))
    set(handles.leftedge,'String',num2str(x1))
    [~,minx1] = min(abs(x1-handles.x.(headers{1}).(chr)(1,:)));
    [~,minx2] = min(abs(x2-handles.x.(headers{1}).(chr)(1,:)));
    win = minx1:minx2+1;
elseif chrcoor == 1
    x1 = str2double(get(handles.leftedge,'String'));
    x2 = str2double(get(handles.rightedge,'String'));
    if x2 < x1
        errordlg('Right edge is less than left edge. Check coordinates')
    end
    chr = get(handles.chr,'String');
    mid = mean(x1,x2);
    y = find(strcmp(chr,handles.x.knownGene(:,4)));
    [~,zz] = min(abs(cell2mat(handles.x.knownGene(y,6))-mid));
    g = y(zz);
    gene = char(handles.x.knownGene{g,2});
    set(handles.gene,'String',gene);
    [~,minx1] = min(abs(x1-handles.x.(headers{1}).(chr)(1,:)));
    [~,minx2] = min(abs(x2-handles.x.(headers{1}).(chr)(1,:)));
    win = minx1:minx2+1;
else
    errordlg('Check one of the 2 search options')
end

%% Starts the Visualization
top = zeros(f,1);
for j = 1:length(headers)
    if ~strcmp(headers{j},'PI') && ~strcmp(headers{j},'MEF_PI')
        subplot(f,1,j),
        if Yaxis == 2
            top(j) = max(handles.x.(headers{j}).(chr)(2,win)/ReadSum(j));
        else
            top(j) = max(handles.x.(headers{j}).(chr)(2,win));
        end
        Peaks = [get(handles.Macs,'Value') get(handles.preimmune,'Value') get(handles.internal,'Value')];
        
        z = .1*top(j);
        if Peaks(1)
            mac = find(strcmp(chr,handles.x.(headers{j}).Macs_text(1:end-1,1)));
            macs = find(str2double(handles.x.(headers{j}).Macs_text(mac,2)) > x1-20000 & str2double(handles.x.(headers{j}).Macs_text(mac,3)) < x2);
            if isempty(macs) == 0
                for i = 1:length(macs)
                    macwin = str2double(handles.x.(headers{j}).Macs_text(mac(macs(i)),2)):str2double(handles.x.(headers{j}).Macs_text(mac(macs(i)),3));
                    plot([macwin(1) macwin(end)],[top(j)+z top(j)+z],'g-','LineWidth',3), hold on
                end
                z = z + .1*top(j);
            end
        end
        if Peaks(2)
            lin_c = find(strcmp(chr,handles.x.(headers{j}).Lin_Conden_text(1:end-1,5)));
            linc = find(handles.x.(headers{j}).Lin_Conden(lin_c,4) > x1-20000 & handles.x.(headers{j}).Lin_Conden(lin_c,5) < x2);
            if isempty(linc) == 0
                for i = 1:length(linc)
                    cond = handles.x.(headers{j}).Lin_Conden(lin_c(linc(i)),4):handles.x.(headers{j}).Lin_Conden(lin_c(linc(i)),5);
                    plot([cond(1) cond(end)],[top(j)+z top(j)+z],'b-','LineWidth',3), hold on
                end
            end
            z = z + .1*top(j);
        end
        if Peaks(3)
            poi_c = find(strcmp(chr,handles.x.(headers{j}).Poisson_Conden_text(1:end-1,5)));
            poic = find(handles.x.(headers{j}).Poisson_Conden(poi_c,4) > x1-20000 & handles.x.(headers{j}).Poisson_Conden(poi_c,5) < x2);
            if isempty(poic) == 0
                for i = 1:length(poic)
                    cond = handles.x.(headers{j}).Poisson_Conden(poi_c(poic(i)),4):handles.x.(headers{j}).Poisson_Conden(poi_c(poic(i)),5);
                    plot([cond(1) cond(end)],[top(j)+z top(j)+z],'k-','LineWidth',3), hold on
                end
            end
            z = z + .1*top(j);
        end
        
        tss = cell2mat(handles.x.knownGene(g,6));
        tes = cell2mat(handles.x.knownGene(g,7));
        if cell2mat(handles.x.knownGene(g,5)) == 1
            plot(tss*ones(1,2),[0 top(j)],'go-',tes*ones(1,2),[0 top(j)],'ro-','LineWidth',2),hold on
        else
            plot(tss*ones(1,2),[0 top(j)],'ro-',tes*ones(1,2),[0 top(j)],'go-','LineWidth',2),hold on
        end
        if Yaxis == 2
            area(handles.x.(headers{j}).(chr)(1,win),handles.x.(headers{j}).(chr)(2,win)./ReadSum(j))
        else
            area(handles.x.(headers{j}).(chr)(1,win),handles.x.(headers{j}).(chr)(2,win))
        end
        ylabel(sprintf('%s Reads',headers{j}))
        axis([x1 x2 0 top(j)+z])
        string = get(gca,'XTick');
        for ZZ = 1:length(string)
            labels(ZZ) = {num2str(string(ZZ))};
        end
        set(gca,'XTickLabel',labels)
        hold off
    else
        subplot(f,1,j)
        if Yaxis == 2
            top(j) = max(handles.x.(headers{j}).(chr)(2,win)/ReadSum(j));
        else
            top(j) = max(handles.x.(headers{j}).(chr)(2,win));
        end
        tss = cell2mat(handles.x.knownGene(g,6));
        tes = cell2mat(handles.x.knownGene(g,7));
        if cell2mat(handles.x.knownGene(g,5)) == 1
            plot(tss*ones(1,2),[0 top(j)],'go-',tes*ones(1,2),[0 top(j)],'ro-','LineWidth',2),hold on
        else
            plot(tss*ones(1,2),[0 top(j)],'ro-',tes*ones(1,2),[0 top(j)],'go-','LineWidth',2), hold on
        end
        if Yaxis == 2
            area(handles.x.(headers{j}).(chr)(1,win),handles.x.(headers{j}).(chr)(2,win)./ReadSum(j))
        else
            area(handles.x.(headers{j}).(chr)(1,win),handles.x.(headers{j}).(chr)(2,win))
        end
        ylabel(sprintf('%s Reads',headers{j}))
        axis([x1 x2 0 top(j)+.1*top(j)])
        string = get(gca,'XTick');
        for ZZ = 1:length(string)
            labels(ZZ) = {num2str(string(ZZ))};
        end
        set(gca,'XTickLabel',labels)
        hold off
    end
end
if axisoption == 2
    [~,j] = max(top);
    for k = 1:f
        subplot(f,1,k)
        axis([x1 x2 0 top(j)+(sum(Peaks)+3)*(.1*top(j))]),hold off
        string = get(gca,'XTick');
        for ZZ = 1:length(string)
            labels(ZZ) = {num2str(string(ZZ))};
        end
        set(gca,'XTickLabel',labels)
    end
end


% --- Executes during object creation, after setting all properties.
function gene_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gene (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function chr_Callback(hObject, eventdata, handles)
% hObject    handle to chr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of chr as text
%        str2double(get(hObject,'String')) returns contents of chr as a double


% --- Executes during object creation, after setting all properties.
function chr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function leftedge_Callback(hObject, eventdata, handles)
% hObject    handle to leftedge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of leftedge as text
%        str2double(get(hObject,'String')) returns contents of leftedge as a double


% --- Executes during object creation, after setting all properties.
function leftedge_CreateFcn(hObject, eventdata, handles)
% hObject    handle to leftedge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rightedge_Callback(hObject, eventdata, handles)
% hObject    handle to rightedge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rightedge as text
%        str2double(get(hObject,'String')) returns contents of rightedge as a double


% --- Executes during object creation, after setting all properties.
function rightedge_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rightedge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chrsearch.
function chrsearch_Callback(hObject, eventdata, handles)
% hObject    handle to chrsearch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chrsearch


% --- Executes on button press in checkgene.
function checkgene_Callback(hObject, eventdata, handles)
% hObject    handle to checkgene (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkgene


% --------------------------------------------------------------------
function Axes_Callback(hObject, eventdata, handles)
% hObject    handle to Axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function LowestIP_Callback(hObject, eventdata, handles)
% hObject    handle to LowestIP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function HighestIP_Callback(hObject, eventdata, handles)
% hObject    handle to HighestIP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function IndividualIP_Callback(hObject, eventdata, handles)
% hObject    handle to IndividualIP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in internal.
function internal_Callback(hObject, eventdata, handles)
% hObject    handle to internal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of internal


% --- Executes on button press in preimmune.
function preimmune_Callback(hObject, eventdata, handles)
% hObject    handle to preimmune (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of preimmune


% --- Executes on button press in Macs.
function Macs_Callback(hObject, eventdata, handles)
% hObject    handle to Macs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Macs


% --- Executes on button press in Nanog.
function Nanog_Callback(hObject, eventdata, handles)
% hObject    handle to Nanog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Nanog


% --- Executes on button press in PI.
function PI_Callback(hObject, eventdata, handles)
% hObject    handle to PI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PI


% --- Executes on button press in Rad23b_Repeat.
function Rad23b_Repeat_Callback(hObject, eventdata, handles)
% hObject    handle to Rad23b_Repeat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Rad23b_Repeat


% --- Executes on button press in PI_Repeat.
function PI_Repeat_Callback(hObject, eventdata, handles)
% hObject    handle to PI_Repeat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PI_Repeat


% --- Executes on selection change in RawReads.
function RawReads_Callback(hObject, eventdata, handles)
% hObject    handle to RawReads (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns RawReads contents as cell array
%        contents{get(hObject,'Value')} returns selected item from RawReads


% --- Executes during object creation, after setting all properties.
function RawReads_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RawReads (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on mouse press over axes background.
function FIGURE_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to FIGURE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

filename = get(handles.filename,'String');
saveas(gcf,filename,'psc2')


function filename_Callback(hObject, eventdata, handles)
% hObject    handle to filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filename as text
%        str2double(get(hObject,'String')) returns contents of filename as a double


% --- Executes during object creation, after setting all properties.
function filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
