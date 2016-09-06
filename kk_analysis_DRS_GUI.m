% kk_analysis_DRS_GUI.m is the graphical user interface for
% kk_analysis_DRS_.m. For details about the program see
% kk_analysis_DRS_GUI.m

function varargout = kk_analysis_DRS_GUI(varargin)
% KK_ANALYSIS_DRS_GUI MATLAB code for kk_analysis_DRS_GUI.fig
%      KK_ANALYSIS_DRS_GUI, by itself, creates a new KK_ANALYSIS_DRS_GUI or raises the existing
%      singleton*.
%
%      H = KK_ANALYSIS_DRS_GUI returns the handle to a new KK_ANALYSIS_DRS_GUI or the handle to
%      the existing singleton*.
%
%      KK_ANALYSIS_DRS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in KK_ANALYSIS_DRS_GUI.M with the given input arguments.
%
%      KK_ANALYSvIS_DRS_GUI('Property','Value',...) creates a new KK_ANALYSIS_DRS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before kk_analysis_DRS_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to kk_analysis_DRS_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help kk_analysis_DRS_GUI

% Last Modified by GUIDE v2.5 27-Jun-2016 09:59:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @kk_analysis_DRS_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @kk_analysis_DRS_GUI_OutputFcn, ...
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


% --- Executes just before kk_analysis_DRS_GUI is made visible.
function kk_analysis_DRS_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to kk_analysis_DRS_GUI (see VARARGIN)
assignin('base', 'xxx', 0.5)

% Choose default command line output for kk_analysis_DRS_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes kk_analysis_DRS_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = kk_analysis_DRS_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uistack(handles.uipanel2)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uistack(handles.uipanel1)

% --- Executes on buuiptton press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
[baseName, folder] = uigetfile('drs_ml.txt');
try
fullFileName = fullfile(folder)
set(handles.text2,'String', [fullFileName(1:12),'...',...
    fullFileName(max(13,numel(fullFileName)-12):numel(fullFileName))]);
disp('button 3')
assignin('base', 'pathname_data', fullFileName);
evalin('base', 'settings.pathname_data=pathname_data');
catch
    disp('no file selected.')
end

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fontSize = 20; % Whatever you want.
% First for the first axes:
try
axes(handles.axes1);
evalin('base', 'clear subplot_gII');
catch
end

handles.pushbutton3.Enable='off';
handles.pushbutton6.Enable='off';

handles.edit_transparent_range_max.Enable='off';
handles.edit_max_small_absorption.Enable='off';
handles.edit_energy_min.Enable='off';
handles.edit_energy_max.Enable='off';
handles.edit_d_energy.Enable='off';
handles.edit_thickness_range_min.Enable='off';
handles.edit_thickness_range_max.Enable='off';
handles.popupmenu_substrate.Enable='off';
handles.doxide.Enable='off';
handles.pushbutton_load_settings.Enable='off';
handles.pushbutton_load_advanced.Enable='off';

handles.checkbox_use_initial_parameters.Enable='off';
handles.checkbox_use_previous_spectrum.Enable='off';

handles.popupmenu_select_settings.Enable='off';
set(handles.popupmenu_select_settings,'Value', 1);
popupmenu_select_settings_Callback(hObject, eventdata, handles)

assignin('base', 'run', 1);
try
settings=evalin('base', 'settings');
catch
end

initial_parms=struct('use', ...
               get(handles.checkbox_use_initial_parameters,'Value'),...
               'use_parms_of_previous_thickness',...
               get(handles.checkbox_use_previous_spectrum,'Value'));
settings.initial_parms=initial_parms;

settings.transparent_range.min=get(handles.edit_transparent_range_min,'String');
settings.transparent_range.max=get(handles.edit_transparent_range_max,'String');
settings.transparent_range.max_small_abs=get(handles.edit_max_small_absorption,'String');

settings.energy.min=get(handles.edit_energy_min,'String');
settings.energy.max=get(handles.edit_energy_max,'String');
settings.energy.d=get(handles.edit_d_energy,'String');

settings.thickness_range.min=get(handles.edit_thickness_range_min,'String');
settings.thickness_range.max=get(handles.edit_thickness_range_max,'String');

settings.substrate.substrate=get(handles.popupmenu_substrate,'Value');
settings.substrate.d_oxide=get(handles.doxide,'String');

settings.fit_parm.weights.st=get(handles.fit_parm_weights_st,'String');
settings.fit_parm.weights.st_diff_c=get(handles.fit_parm_weights_st_diff_c,'String');
settings.fit_parm.weights.st_del=get(handles.fit_parm_weights_st_del,'String');
settings.fit_parm.weights.diff_exp_mod_transparent_range=get(handles.fit_parm_weights_diff_exp_mod_transparent_range,'String');
settings.fit_parm.weights.e2_curv_transparent_range=get(handles.fit_parm_weights_e2_curv_transparent_range,'String');
settings.fit_parm.weight_e2_small_factor=get(handles.fit_parm_weight_e2_small_factor,'String');

set(handles.pathname_fit,'Visible', 'off');
set(handles.text_results,'Visible', 'off');

drawnow
try
kk_analysis_DRS_(settings);
catch error
    assignin('base', 'error', error)
    if strcmp(error.message, 'Undefined function or variable ''stopped''.')
        disp('stopped by button')
    else
        pushbutton5_Callback(hObject, eventdata, handles)
        error.rethrow
    end
end
disp('button 4')
assignin('base', 'run', 1);

pushbutton5_Callback(hObject, eventdata, handles)

% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('button 5; stop')
try
assignin('base', 'run', 0);
pathname_fit=evalin('base', 'pathname_fit');
evalin('base', 'pathname_fit_old=pathname_fit');
set(handles.pathname_fit,'String', pathname_fit);
set(handles.pathname_fit,'Visible', 'on');
set(handles.text_results,'Visible', 'on');
evalin('base', 'clear pathname_fit');
end
try
handles.pushbutton3.Enable='on';
handles.pushbutton6.Enable='on';

handles.edit_transparent_range_max.Enable='on';
handles.edit_max_small_absorption.Enable='on';
handles.edit_energy_min.Enable='on';
handles.edit_energy_max.Enable='on';
handles.edit_d_energy.Enable='on';
handles.edit_thickness_range_min.Enable='on';
handles.edit_thickness_range_max.Enable='on';
handles.popupmenu_substrate.Enable='on';
handles.doxide.Enable='on';
popupmenu_substrate_Callback(hObject, eventdata, handles);
handles.pushbutton_load_settings.Enable='on';
handles.pushbutton_load_advanced.Enable='on';

handles.checkbox_use_initial_parameters.Enable='on';
handles.checkbox_use_previous_spectrum.Enable='on';
checkbox_use_initial_parameters_Callback(hObject, eventdata, handles);

handles.popupmenu_select_settings.Enable='on';
end

drawnow


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_use_initial_parameters.
function checkbox_use_initial_parameters_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_use_initial_parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_use_initial_parameters
s=get(handles.checkbox_use_initial_parameters,'Value');
if s
handles.checkbox_use_previous_spectrum.Enable='on';
else
handles.checkbox_use_previous_spectrum.Enable='off';
end


% --- Executes on selection change in popupmenu_substrate.
function popupmenu_substrate_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_substrate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_substrate contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_substrate
substrate=handles.popupmenu_substrate.Value
if substrate ==1
    handles.text5.Enable='off';
    handles.doxide.Enable='off';
    doxide=0;
end

if substrate ==2
    handles.text5.Enable='on';
    handles.doxide.Enable='on';
    doxide=0;
end
% --- Executes during object creation, after setting all properties.
function popupmenu_substrate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_substrate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function doxide_Callback(hObject, eventdata, handles)
% hObject    handle to doxide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of doxide as text
%        str2double(get(hObject,'String')) returns contents of doxide as a double


% --- Executes during object creation, after setting all properties.
function doxide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to doxide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_thickness_range_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_thickness_range_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_thickness_range_min as text
%        str2double(get(hObject,'String')) returns contents of edit_thickness_range_min as a double
edit_thickness_range_min=str2double(get(handles.edit_thickness_range_min,'String'));
edit_thickness_range_max=str2double(get(handles.edit_thickness_range_max,'String'));
set(handles.edit_thickness_range_max, 'String', max(edit_thickness_range_min,edit_thickness_range_max));


% --- Executes during object creation, after setting all properties.
function edit_thickness_range_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_thickness_range_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_thickness_range_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_thickness_range_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_thickness_range_max as text
%        str2double(get(hObject,'String')) returns contents of edit_thickness_range_max as a double
edit_thickness_range_min=str2double(get(handles.edit_thickness_range_min,'String'));
edit_thickness_range_max=str2double(get(handles.edit_thickness_range_max,'String'));
set(handles.edit_thickness_range_min, 'String', min(edit_thickness_range_min,edit_thickness_range_max));

% --- Executes during object creation, after setting all properties.
function edit_thickness_range_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_thickness_range_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_transparent_range_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_transparent_range_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_transparent_range_min as text
%        str2double(get(hObject,'String')) returns contents of edit_transparent_range_min as a double


% --- Executes during object creation, after setting all properties.
function edit_transparent_range_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_transparent_range_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_transparent_range_max_Callback(hObject, eventdata, handles)
energy_min=str2double(get(handles.edit_energy_min,'String'));
energy_max=str2double(get(handles.edit_energy_max,'String'));
edit_transparent_range_max=str2double(get(handles.edit_transparent_range_max,'String'));
edit_max_small_absorption=str2double(get(handles.edit_max_small_absorption,'String'));

energy_min=min(energy_min, edit_transparent_range_max-0.05);
set(handles.edit_transparent_range_min, 'String', num2str(energy_min));
set(handles.edit_energy_min           , 'String', num2str(energy_min));
set(handles.edit_energy_max           , 'String', num2str(max([energy_max,energy_min+0.8,edit_transparent_range_max+0.5+0.02])));
set(handles.edit_max_small_absorption , 'String', num2str(max(edit_max_small_absorption,edit_transparent_range_max+0.02)));

% --- Executes during object creation, after setting all properties.
function edit_transparent_range_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_transparent_range_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_energy_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_energy_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_energy_min as text
%        str2double(get(hObject,'String')) returns contents of edit_energy_min as a double
energy_min=str2double(get(handles.edit_energy_min,'String'));
energy_max=str2double(get(handles.edit_energy_max,'String'));
edit_transparent_range_max=str2double(get(handles.edit_transparent_range_max,'String'));
edit_max_small_absorption=str2double(get(handles.edit_max_small_absorption,'String'));
if (energy_min+0.1)>(edit_max_small_absorption+0.0);
warningstring='The (almost) transparent range should cover at least 0.1 eV';
h = warndlg(warningstring);
end
set(handles.edit_transparent_range_min, 'String', num2str(energy_min));
set(handles.edit_energy_max           , 'String', num2str(max(energy_max,energy_min+0.8)));
set(handles.edit_transparent_range_max, 'String', num2str(max(edit_transparent_range_max,energy_min+0.05)));
edit_transparent_range_max=str2double(get(handles.edit_transparent_range_max,'String'));
set(handles.edit_max_small_absorption, 'String', num2str(max(edit_max_small_absorption,edit_transparent_range_max+0.02)));



% --- Executes during object creation, after setting all properties.
function edit_energy_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_energy_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_energy_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_energy_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_energy_max as text
%        str2double(get(hObject,'String')) returns contents of edit_energy_max as a double
energy_min=str2double(get(handles.edit_energy_min,'String'));
energy_max=str2double(get(handles.edit_energy_max,'String'));
set(handles.edit_energy_min           , 'String', num2str(min(energy_min, energy_max-0.8)));
set(handles.edit_transparent_range_min, 'String', num2str(min(energy_min, energy_max-0.8)));
edit_energy_min_Callback(hObject, eventdata, handles)
edit_max_small_absorption=min(str2double(get(handles.edit_max_small_absorption,'String')),energy_max-0.5);
set(handles.edit_max_small_absorption, 'String', num2str(edit_max_small_absorption));
edit_max_small_absorption_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit_energy_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_energy_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_d_energy_Callback(hObject, eventdata, handles)
% hObject    handle to edit_d_energy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_d_energy as text
%        str2double(get(hObject,'String')) returns contents of edit_d_energy as a double
edit_d_energy=str2double(get(handles.edit_d_energy,'String'));
set(handles.edit_d_energy, 'String', num2str(min(max(0.001,edit_d_energy),0.05)));

% --- Executes during object creation, after setting all properties.
function edit_d_energy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_d_energy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles, fullFileName)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% initial parms
if nargin==3
[baseName, folder] = uigetfile('*.*');
fullFileName = fullfile(folder, baseName);
end
try
set(handles.text22,'String', [fullFileName(1:12),'...',...
    fullFileName(max(13,numel(fullFileName)-12):numel(fullFileName))]);
[~,data]=hdrload(fullFileName);
E_e1_e2=data(:,1:3);
min_e1=max(0,(max(E_e1_e2(:,2))-1.5));
try
subplot(evalin('base', 'subplot_g'));
catch
axes(handles.axes4);
end
g=subplot(1,1,1);
assignin('base', 'subplot_g', g);
p = get(g,'position');
    p(3) = p(3)*0.7;
    p(1) = p(1)+p(3)*0.3;
    p(4) = p(4)*0.5;
    p(2) = p(2)+p(4)*0.8;
    set(g, 'position', p);
plot(E_e1_e2(:,1), [E_e1_e2(:,2)-min_e1,E_e1_e2(:,3)])
numstring=[num2str(min_e1),'0000'];
legend(['e1xy-',numstring(1:3)],'e2xy')    
axis([min(E_e1_e2(:,1)) max(E_e1_e2(:,1)) 0 1.5])
xlabel('Energy [eV]','FontSize',11)
ylabel(['\epsilon _{2},\epsilon _{1}-', numstring(1:3)],'FontSize',11)
assignin('base', 'pathname_initial', fullFileName);
evalin('base', 'settings.pathname_initial=pathname_initial');
catch
disp('no file selected')    
end

% --- Executes on button press in exit.
function exit_Callback(hObject, eventdata, handles)
% hObject    handle to exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
assignin('base', 'run', 0);
delete(gcbf)


% --- Executes on button press in checkbox_use_previous_spectrum.
function checkbox_use_previous_spectrum_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_use_previous_spectrum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_use_previous_spectrum
%pushbutton4_Callback(hObject, eventdata, handles)

% --- Executes on button press in pushbutton_load_settings.
function pushbutton_load_settings_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load_settings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[baseName, folder] = uigetfile('settings.mat');
fullFileName = fullfile(folder,baseName);
try
settings_=load(fullFileName);
settings=settings_.settings;

set(handles.edit_transparent_range_min, 'String', settings.transparent_range.min);
set(handles.edit_transparent_range_max, 'String', settings.transparent_range.max);
set(handles.edit_max_small_absorption, 'String', settings.transparent_range.max_small_abs);
set(handles.edit_energy_min, 'String', settings.energy.min);
set(handles.edit_energy_max, 'String', settings.energy.max);
set(handles.edit_d_energy, 'String',   settings.energy.d);
set(handles.edit_thickness_range_min, 'String', settings.thickness_range.min);
set(handles.edit_thickness_range_max, 'String', settings.thickness_range.max);
set(handles.popupmenu_substrate, 'Value', settings.substrate.substrate);
set(handles.doxide, 'String', settings.substrate.d_oxide);
set(handles.checkbox_use_initial_parameters, 'Value', settings.initial_parms.use);
set(handles.checkbox_use_previous_spectrum,  'Value', settings.initial_parms.use_parms_of_previous_thickness);

%set advanced settings:
set(handles.fit_parm_weights_st, 'String', settings.fit_parm.weights.st);
set(handles.fit_parm_weights_st_diff_c, 'String', settings.fit_parm.weights.st);
set(handles.fit_parm_weights_st_del, 'String', settings.fit_parm.weights.st_del);
set(handles.fit_parm_weights_diff_exp_mod_transparent_range, 'String', settings.fit_parm.weights.diff_exp_mod_transparent_range);
set(handles.fit_parm_weights_e2_curv_transparent_range, 'String', settings.fit_parm.weights.e2_curv_transparent_range);
set(handles.fit_parm_weight_e2_small_factor, 'String', settings.fit_parm.weight_e2_small_factor);

try
assignin('base', 'pathname_data', settings.pathname_data);
evalin('base', 'settings.pathname_data=pathname_data');
end
try
assignin('base', 'pathname_initial', settings.pathname_initial);
evalin('base', 'settings.pathname_initial=pathname_initial');
pushbutton6_Callback(hObject, eventdata, handles, settings.pathname_initial);
end

checkbox_use_initial_parameters_Callback(hObject, eventdata, handles)
catch
    disp('error: no initial parameters loaded.')
end

function edit_max_small_absorption_Callback(hObject, eventdata, handles)
energy_min=str2double(get(handles.edit_energy_min,'String'));
energy_max=str2double(get(handles.edit_energy_max,'String'));
edit_transparent_range_max=str2double(get(handles.edit_transparent_range_max,'String'));
edit_max_small_absorption=str2double(get(handles.edit_max_small_absorption,'String'));


if (energy_min+0.1)>(edit_max_small_absorption+0.0);
warningstring='The (almost) transparent range should cover at least 0.1 eV';
h = warndlg(warningstring);
end


energy_min=min(energy_min, edit_max_small_absorption-0.05-0.02);
set(handles.edit_transparent_range_min, 'String', num2str(energy_min));
set(handles.edit_energy_min           , 'String', num2str(energy_min));
set(handles.edit_energy_max           , 'String', num2str(max([energy_max,energy_min+0.8,edit_max_small_absorption+0.5])));
set(handles.edit_transparent_range_max, 'String', num2str(min(edit_transparent_range_max,edit_max_small_absorption-0.02)));


% --- Executes during object creation, after setting all properties.
function edit_max_small_absorption_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_max_small_absorption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes on selection change in popupmenu_select_settings.
function popupmenu_select_settings_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_select_settings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_select_settings contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_select_settings

select_=get(handles.popupmenu_select_settings,'Value');

switch select_
    case 1
        advanced_visible='off';
        string_='show advanced settings';
    case 2
        advanced_visible='on';
        string_='hide advanced settings';
end

assignin('base','advanced_visible',advanced_visible);
set(handles.uipanel_advanced_settings, 'visible', advanced_visible);

% --- Executes during object creation, after setting all properties.
function popupmenu_select_settings_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_select_settings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fit_parm_weights_st_Callback(hObject, eventdata, handles)
% hObject    handle to fit_parm_weights_st (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fit_parm_weights_st as text
%        str2double(get(hObject,'String')) returns contents of fit_parm_weights_st as a double


% --- Executes during object creation, after setting all properties.
function fit_parm_weights_st_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fit_parm_weights_st (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fit_parm_weights_st_diff_c_Callback(hObject, eventdata, handles)
% hObject    handle to fit_parm_weights_st_diff_c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fit_parm_weights_st_diff_c as text
%        str2double(get(hObject,'String')) returns contents of fit_parm_weights_st_diff_c as a double


% --- Executes during object creation, after setting all properties.
function fit_parm_weights_st_diff_c_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fit_parm_weights_st_diff_c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fit_parm_weights_st_del_Callback(hObject, eventdata, handles)
% hObject    handle to fit_parm_weights_st_del (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fit_parm_weights_st_del as text
%        str2double(get(hObject,'String')) returns contents of fit_parm_weights_st_del as a double


% --- Executes during object creation, after setting all properties.
function fit_parm_weights_st_del_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fit_parm_weights_st_del (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fit_parm_weight_e2_small_factor_Callback(hObject, eventdata, handles)
% hObject    handle to fit_parm_weight_e2_small_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fit_parm_weight_e2_small_factor as text
%        str2double(get(hObject,'String')) returns contents of fit_parm_weight_e2_small_factor as a double


% --- Executes during object creation, after setting all properties.
function fit_parm_weight_e2_small_factor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fit_parm_weight_e2_small_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fit_parm_weights_e2_curv_transparent_range_Callback(hObject, eventdata, handles)
% hObject    handle to fit_parm_weights_e2_curv_transparent_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fit_parm_weights_e2_curv_transparent_range as text
%        str2double(get(hObject,'String')) returns contents of fit_parm_weights_e2_curv_transparent_range as a double


% --- Executes during object creation, after setting all properties.
function fit_parm_weights_e2_curv_transparent_range_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fit_parm_weights_e2_curv_transparent_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fit_parm_weights_diff_exp_mod_transparent_range_Callback(hObject, eventdata, handles)
% hObject    handle to fit_parm_weights_diff_exp_mod_transparent_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fit_parm_weights_diff_exp_mod_transparent_range as text
%        str2double(get(hObject,'String')) returns contents of fit_parm_weights_diff_exp_mod_transparent_range as a double


% --- Executes during object creation, after setting all properties.
function fit_parm_weights_diff_exp_mod_transparent_range_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fit_parm_weights_diff_exp_mod_transparent_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pathname_fit_Callback(hObject, eventdata, handles)
% hObject    handle to pathname_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pathname_fit as text
%        str2double(get(hObject,'String')) returns contents of pathname_fit as a double
pathname_fit=evalin('base','pathname_fit_old');
set(handles.pathname_fit,'String', [pathname_fit]);
