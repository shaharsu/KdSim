function varargout = KdSimGUI(varargin)
% KDSIMGUI MATLAB code for KdSimGUI.fig
%      KdSimGUI is a GUI to assisst the user in generating simuated
%      cellular volume change data. For details, please see the original
%      paper.

% Last Modified by GUIDE v2.5 23-Mar-2017 20:34:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @KdSimGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @KdSimGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
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


% --- Executes just before KdSimGUI is made visible.
function KdSimGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to KdSimGUI (see VARARGIN)

% Choose default command line output for KdSimGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes KdSimGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = KdSimGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)

%---------------------------------
% Executes on pushing the "run" button
%----------------------------------

n_runs=str2double(get(handles.simNo,'string'))*6;
k_on=str2double(get(handles.kon,'string'));
k_off=str2double(get(handles.koff,'string'));
set(handles.Kd,'string',num2str(k_off/k_on));
% A_in=str2double(get(handles.A_conc,'string'));
% B_in=str2double(get(handles.B_conc,'string'));
stoiA=str2double(get(handles.stoiA,'string'));
stoiB=str2double(get(handles.stoiB,'string'));
E_C=str2double(get(handles.E_C,'string'));
R_i=zeros(n_runs,1);
G_i=R_i;
dG=R_i;
dR=R_i;
dV=R_i;ks_on=R_i;ks_off=R_i;
A_i=R_i;B_i=R_i;C_i=R_i;

for i=1:n_runs
    [A_in,B_in]=deal(-1);
    while (A_in<0)
        A_in = (str2double(get(handles.A_conc,'string'))*randn(1)+6)*1e-6; % A from experimental distribution (sigma=6)
    end
    while (B_in<0)
        B_in = (str2double(get(handles.B_conc,'string'))*randn(1)+6)*1e-6; % A from experimental distribution (sigma=6)
    end
    V_choice=[1.25 1.07 1 0.92 0.83 0.74];
    dV_i = V_choice(mod(i,6)+1)+randn()*0.01; % dV from experimental values
    start=[A_in B_in 0 dV_i k_on k_off E_C stoiA stoiB]; % run_simpleKd input
    [R_i(i,1),G_i(i,1),dR(i,1),dG(i,1),dF(i,1),dV(i,1),A_i(i,1),B_i(i,1),C_i(i,1),ks_on(i,1),ks_off(i,1)]=run_simpleKd(start);
end

% These are experimental values that can be plotted on top of the
% simulations for comparison. Just uncomment lines beginning with h1p1,
% h2p1, h3p1.

% AcGFP1-mCherry
% dG_exp=[0.15612 0.02111 0 -0.02719 -0.03314 -0.09905];
% dR_exp=[-0.08422 -0.02823 0 0.03719 0.10335 0.15899];
% dF_exp=[-0.16535 -0.04172 0 0.05331 0.12445 0.22032];

% GFP-GAPDH + PGK-mCherry
dG_exp = [0.03 0.001 0 -0.02 -0.03 -0.06];
dR_exp = [-0.04 -0.02 0 0.02 0.1 0.15];
dF_exp = [-0.17 -0.04 0 0.05 0.12 0.22];

h1=scatter(handles.axes1,dV,dG,'markeredgecolor','black','linewidth',1,'markerfacecolor','green');
hold on
% h1p1=line(handles.axes1,V_choice,dG_exp,'linestyle','none','marker','o','markeredgecolor','black','markerfacecolor','red','markersize',15);
hold off
ylabel(handles.axes1,'\chi_{green}');
set(handles.axes1,'xticklabel','');

h2=scatter(handles.axes6,dV,dR,'markeredgecolor','black','linewidth',1,'markerfacecolor','red');
hold on
% h2p1=line(handles.axes6,V_choice,dR_exp,'linestyle','none','marker','o','markeredgecolor','black','markerfacecolor','red','markersize',15);
hold off
ylabel(handles.axes6,'\chi_{red}');
set(handles.axes6,'xticklabel','');

h3=scatter(handles.axes7,dV,dF,'markeredgecolor','black','linewidth',1,'markerfacecolor',[1 0.6 0]);
hold on
% h3p1=line(handles.axes7,V_choice,dF_exp,'linestyle','none','marker','o','markeredgecolor','black','markerfacecolor','red','markersize',15);
hold off
ylabel(handles.axes7,'\chi_{FRET}');
xlabel(handles.axes7,'$\tilde{V}$','Interpreter','LaTex','FontSize',20);

set([handles.axes1,handles.axes6,handles.axes7],'userdata',[R_i,G_i,dR,dG,dF,dV,A_i,B_i,C_i,ks_on,ks_off],'colororder',[(dV-min(dV))/(max(dV)-min(dV)) dV-dV 1-(dV-min(dV))/(max(dV)-min(dV))]);
set([h1,h2,h3],'ButtonDownFcn',@plotPoint);



function plotPoint(src,evntData)
%----------------------------------------
% plot point when clicked on either green or red
%----------------------------------------
ax = get(src,'Parent');
ax2 = findobj(gcf,'tag','axes2');
axG = findobj(gcf,'tag','axes3');
axR = findobj(gcf,'tag','axes5');
axF = findobj(gcf,'tag','axes4');
table = findobj(gcf,'tag','pointTable');
E_C = str2double(get(findobj(gcf,'tag','E_C'),'string'));
kon = str2double(get(findobj(gcf,'tag','kon'),'string'));
koff = str2double(get(findobj(gcf,'tag','koff'),'string'));
stoiA=str2double(get(findobj(gcf,'tag','stoiA'),'string'));
stoiB=str2double(get(findobj(gcf,'tag','stoiB'),'string'));
Cp = get(ax,'CurrentPoint');
allData = get(ax,'userdata');
X = get(src,'xdata');
Y = get(src,'ydata');
Xp = Cp(2,1);  % X-point
Yp = Cp(2,2);  % Y-point
[dp,Ip] = min((X-Xp).^2+(Y-Yp).^2);

% Extract the closest coordinate values

Xc = X(Ip);
Yc = Y(Ip);
Xpos = find((X-Xc)==0);

% Add some text

str1 = sprintf('\\leftarrow %i',Xpos);
ht(1) = text(Xp,Yp,str1,'Clipping','off');
set(ht,'FontSize',8,'Color','red')

point = (allData(Xpos,:)); %point = [R_i,G_i,dR,dG,dV,A_i,B_i,C_i,k_on,k_off]
start = [point(7),point(8),point(9),point(6),kon,koff,E_C,stoiA,stoiB]; %start = [A_i,B_i,C_i,dV,k_on,k_off]
[R_i,G_i,dR,dG,dF,dV,A_i,B_i,C_i,k_on,k_off,t,A,B,C,R,G,F] = run_simpleKd(start); 
                                     
delete(findobj([ax2 axR axG axF],'type','line'));

t_ini=900;%9e2;
t_fin=size(t,1);%1500;
h1=line('xdata',t(t_ini:t_fin),'ydata',[A(t_ini:t_fin)],'marker','.','parent',ax2,'color','g');
h2=line('xdata',t(t_ini:t_fin),'ydata',[B(t_ini:t_fin)],'marker','.','parent',ax2,'color','r');
h3=line('xdata',t(t_ini:t_fin),'ydata',[C(t_ini:t_fin)],'marker','.','parent',ax2,'color','b');
title(ax2,'Concentrations')
legend(ax2,'donor','acceptor','complex')
ylabel(ax2,'c (\muM)')
xlabel(ax2,'t (s)')
assignin('base','simul',[t; A; B; C; G; R; F]');

h31=line('xdata',t(t_ini:t_fin),'ydata',[G(t_ini:t_fin)],'marker','.','parent',axG,'color','g');
title(axG,'Green emission')
ylabel(axG,'I (AU)')
xlabel(axG,'t (s)')

h32=line('xdata',t(t_ini:t_fin),'ydata',[R(t_ini:t_fin)],'marker','.','parent',axR,'color','r');
title(axR,'Red emission')
ylabel(axR,'I (AU)')
xlabel(axR,'t (s)')

h33=line('xdata',t(t_ini:t_fin),'ydata',F(t_ini:t_fin),'marker','.','parent',axF,'color','m');
title(axF,'FRET')
ylabel(axF,'E_F (AU)')
xlabel(axF,'t (s)')

tableData=[{point(6)} ;{A(1)} ;{B(1)} ;{point(3)} ;{point(4)}];
set(table,'data',tableData);

function kon_Callback(hObject, eventdata, handles)
% hObject    handle to kon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kon as text
%        str2double(get(hObject,'String')) returns contents of kon as a double


% --- Executes during object creation, after setting all properties.
function kon_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on kon and none of its controls.
function kon_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to kon (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)



function koff_Callback(hObject, eventdata, handles)
% hObject    handle to koff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of koff as text
%        str2double(get(hObject,'String')) returns contents of koff as a double


% --- Executes during object creation, after setting all properties.
function koff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to koff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Kd_Callback(hObject, eventdata, handles)
% hObject    handle to Kd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Kd as text
%        str2double(get(hObject,'String')) returns contents of Kd as a double


% --- Executes during object creation, after setting all properties.
function Kd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Kd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on pushbutton1 and none of its controls.
function pushbutton1_KeyPressFcn(hObject, eventdata, handles)

% hObject    handle to pushbutton1 (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)



function simNo_Callback(hObject, eventdata, handles)
% hObject    handle to simNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of simNo as text
%        str2double(get(hObject,'String')) returns contents of simNo as a double


% --- Executes during object creation, after setting all properties.
function simNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to simNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)




% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function E_C_Callback(hObject, eventdata, handles)
% hObject    handle to E_C (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of E_C as text
%        str2double(get(hObject,'String')) returns contents of E_C as a double


% --- Executes during object creation, after setting all properties.
function E_C_CreateFcn(hObject, eventdata, handles)
% hObject    handle to E_C (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function A_conc_Callback(hObject, eventdata, handles)
% hObject    handle to A_conc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of A_conc as text
%        str2double(get(hObject,'String')) returns contents of A_conc as a double


% --- Executes during object creation, after setting all properties.
function A_conc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A_conc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function B_conc_Callback(hObject, eventdata, handles)
% hObject    handle to B_conc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of B_conc as text
%        str2double(get(hObject,'String')) returns contents of B_conc as a double


% --- Executes during object creation, after setting all properties.
function B_conc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to B_conc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in copyClp.
function copyClp_Callback(hObject, eventdata, handles)

%----------------------------------------
% Executes on "Copy data" button press
% copies all points in chi_red, chi_green, and chi_FRET to clipboard.
% format is [A_i B_i dV chi_G chi_R] for all volume changes.
% First two rows are AVG and SD of the data.
%----------------------------------------

format long;
allData=get(handles.axes1,'userdata');
selData=[allData(:,7),allData(:,8),allData(:,6) allData(:,5) allData(:,4) allData(:,3)]; %dV, dF, dG, dR
A_i=circshift(reshape(selData(:,1)',6,size(selData,1)/6)',[0 1]);
B_i=circshift(reshape(selData(:,2)',6,size(selData,1)/6)',[0 1]);
dV=circshift(reshape(selData(:,3)',6,size(selData,1)/6)',[0 1])-1;
dF=circshift(reshape(selData(:,4)',6,size(selData,1)/6)',[0 1]);
dG=circshift(reshape(selData(:,5)',6,size(selData,1)/6)',[0 1]);
dR=circshift(reshape(selData(:,6)',6,size(selData,1)/6)',[0 1]);
finalDataWorkspace=[dV(:,1) dG(:,1) dR(:,1); dV(:,2) dG(:,2) dR(:,2); dV(:,3) dG(:,3) dR(:,3); dV(:,4) dG(:,4) dR(:,4); dV(:,5) dG(:,5) dR(:,5) ;dV(:,6) dG(:,6) dR(:,6)];
finalData=[dV(:,1) dG(:,1) dR(:,1) dV(:,2) dG(:,2) dR(:,2) dV(:,3) dG(:,3) dR(:,3) dV(:,4) dG(:,4) dR(:,4) dV(:,5) dG(:,5) dR(:,5) dV(:,6) dG(:,6) dR(:,6)];
finalConc=[A_i(:,1) B_i(:,1) A_i(:,2) B_i(:,2) A_i(:,3) B_i(:,3) A_i(:,4) B_i(:,4) A_i(:,5) B_i(:,5) A_i(:,6) B_i(:,6)];
meanConc=reshape(mean(finalConc,'omitnan'),2,6)';
stdConc=reshape(std(finalConc,'omitnan'),2,6)';
meanData=reshape(mean(finalData,'omitnan'),3,6)';
stdData=reshape(std(finalData,'omitnan'),3,6)';
clipData=[meanConc(:,1) stdConc(:,1) meanConc(:,2) stdConc(:,2) meanData(:,1),stdData(:,1),meanData(:,2),stdData(:,2),meanData(:,3),stdData(:,3)];
assignin('base','allData',finalDataWorkspace);
mat2clip(finalDataWorkspace);

% hObject    handle to copyClp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over copyClp.
function copyClp_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to copyClp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in FitInit.
function FitInit_Callback(hObject, eventdata, handles)

%----------------------------------------
% Executes on "Fit" button press
% Not currently implemented.
%----------------------------------------

A=get(handles.axes1,'userdata');
% dG_exp=[0.15612;0.02111;0;-0.02719;-0.03314;-0.09905];
% dR_exp=[-0.08422;-0.02823;0;0.03719;0.10335;0.15899];
% dF_exp=[-0.16535;-0.04172;0;0.05331;0.12445;0.22032];
dG_exp = [0.03;0.001;0;-0.02;-0.03;-0.06];
dR_exp = [-0.04;-0.02;0;0.02;0.1;0.15];
dF_exp = [-0.17;-0.04;0;0.05;0.12;0.22];

dG=reshape(A(:,4),6,34);
dG_ave=circshift(mean(dG,2),1);
sse_dG=sum(sqrt((dG_ave-dG_exp).^2));
dR=reshape(A(:,3),6,34);
dR_ave=circshift(mean(dR,2),1);
sse_dR=sum(sqrt((dR_ave-dR_exp).^2));
dF=reshape(A(:,5),6,34);
dF_ave=circshift(mean(dF,2),1);
sse_dF=sum(sqrt((dF_ave-dF_exp).^2));
dV=reshape(A(:,6),6,34);


function getAvg_Callback(hObject, eventdata, handles)

%----------------------------------------
% Executes on "Averages" button press
% Copy average and SD of the time trace for red, green and FRET 
% fluorescence for all points to clipboard
%----------------------------------------

allData=get(handles.axes1,'userdata');
E_C = str2double(get(findobj(gcf,'tag','E_C'),'string'));
kon = str2double(get(findobj(gcf,'tag','kon'),'string'));
koff = str2double(get(findobj(gcf,'tag','koff'),'string'));
stoiA=str2double(get(findobj(gcf,'tag','stoiA'),'string'));
stoiB=str2double(get(findobj(gcf,'tag','stoiB'),'string'));
point = (allData(1,:)); %point = [R_i,G_i,dR,dG,dV,A_i,B_i,C_i,k_on,k_off]
start = [point(7),point(8),point(9),point(6),kon,koff,E_C,stoiA,stoiB]; %start = [A_i,B_i,C_i,dV,k_on,k_off]
[R_i,G_i,dR,dG,dF,dV,A_i,B_i,C_i,k_on,k_off,t] = run_simpleKd(start);
t_size=size(t,1);
R_raw=zeros(size(allData,1)/6,t_size,6);
G_raw=R_raw;
F_raw=R_raw;
for ii=1:size(allData,1)
    point = (allData(ii,:)); %point = [R_i,G_i,dR,dG,dV,A_i,B_i,C_i,k_on,k_off]
    start = [point(7),point(8),point(9),point(6),kon,koff,E_C,stoiA,stoiB]; %start = [A_i,B_i,C_i,dV,k_on,k_off]
    [R_i,G_i,dR,dG,dF,dV,A_i,B_i,C_i,k_on,k_off,t,A,B,C,R_raw(ceil(ii/6),:,mod(ii,6)+1),G_raw(ceil(ii/6),:,mod(ii,6)+1),F_raw(ceil(ii/6),:,mod(ii,6)+1)] = run_simpleKd(start);
end
R_norm=R_raw./repmat(R_raw(:,round((t_size-2)/2),:),1,t_size,1);
G_norm=G_raw./repmat(G_raw(:,round((t_size-2)/2),:),1,t_size,1);
F_norm=F_raw./repmat(F_raw(:,round((t_size-2)/2),:),1,t_size,1);
for ii=1:size(R_norm,3)
    R_normMean(:,2*ii-1:2*ii)=[mean(squeeze(R_norm(:,:,ii)'),2,'omitnan') std(squeeze(R_norm(:,:,ii)'),0,2,'omitnan')];
    G_normMean(:,2*ii-1:2*ii)=[mean(squeeze(G_norm(:,:,ii)'),2,'omitnan') std(squeeze(G_norm(:,:,ii)'),0,2,'omitnan')];
    F_normMean(:,2*ii-1:2*ii)=[mean(squeeze(F_norm(:,:,ii)'),2,'omitnan') std(squeeze(F_norm(:,:,ii)'),0,2,'omitnan')];
end
mat2clip([t R_normMean(1:t_size,:) G_normMean(1:t_size,:) F_normMean(1:t_size,:)]);



function stoiA_Callback(hObject, eventdata, handles)
% hObject    handle to stoiA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stoiA as text
%        str2double(get(hObject,'String')) returns contents of stoiA as a double


% --- Executes during object creation, after setting all properties.
function stoiA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stoiA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stoiB_Callback(hObject, eventdata, handles)
% hObject    handle to stoiB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stoiB as text
%        str2double(get(hObject,'String')) returns contents of stoiB as a double


% --- Executes during object creation, after setting all properties.
function stoiB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stoiB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
