function varargout = demo(varargin)
% DEMO MATLAB code for demo.fig
%      DEMO, by itself, creates a new DEMO or raises the existing
%      singleton*.
%
%      H = DEMO returns the handle to a new DEMO or the handle to
%      the existing singleton*.
%
%      DEMO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DEMO.M with the given input arguments.
%
%      DEMO('Property','Value',...) creates a new DEMO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before demo_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to demo_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help demo

% Last Modified by GUIDE v2.5 22-Feb-2018 21:03:52

% Begin initialization code - DO NOT EDIT
global ud;

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @demo_OpeningFcn, ...
                   'gui_OutputFcn',  @demo_OutputFcn, ...
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


% --- Executes just before demo is made visible.
function demo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to demo (see VARARGIN)

% Choose default command line output for demo
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
axes(handles.axes8);
imshow(imread('colorbar.png'));
% UIWAIT makes demo wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = demo_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
global ud;
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile({'*.jpg;*.tif;*.png;*.gif;*.bmp','All Image Files';...
          '*.*','All Files' },'Image');
 I = imread([pathname,filename]);
 ud.imgl = I;
 axes(handles.axes1);
 imshow(I);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
global ud;
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile({'*.jpg;*.tif;*.png;*.gif;*.bmp','All Image Files';...
          '*.*','All Files' },'Image');
 J = imread([pathname,filename]);
 ud.imgr = J;
 axes(handles.axes5);
 imshow(J);


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
global ud;
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

topdown = get(get(handles.mg, 'SelectedObject'),'String');
btupm = get(get(handles.btup, 'SelectedObject'),'String');

% value = get(handles.ckbhdr, 'Value');
min = str2num(get(handles.edtmin, 'String'));
max = str2num(get(handles.edtmax, 'String'));
s = str2num(get(handles.edtWNum, 'String'));
ppd = str2num(get(handles.editPPD, 'String'));
threshold_sliency = get(handles.slider2,'Value');

map=[];
mx=0;
mi=0;
jnd=0;
background = rgb2gray(ud.imgl);

%     ppd = hdrvdp_pix_per_deg( 21, [size(0,1) size(0,2)], 1 );

tic;
if strcmp(btupm, 'HDRVDP') == 1
    if strcmp(topdown, 'SSIM') == 1
        res1 = hdrvdp_ssim( ud.imgl, ud.imgr, 'sRGB-display', ppd, [], s, min, max, 'normal' );

        if res1.same == 0
            map = hdrvdp_visualize_ssim(background, res1.P_map, res1.patchPoints);
        else
            map = cat(3, background, background, background);  
        end
        mx=res1.maxval;
        mi=res1.minval;
    elseif strcmp(topdown, 'CW_SSIM') == 1
        res1 = hdrvdp_ssim( ud.imgl, ud.imgr, 'sRGB-display', ppd, [], s, min, max, 'cw' );
        if res1.same == 0
            map = hdrvdp_visualize_ssim(background, res1.P_map, res1.patchPoints);
        else
            map = cat(3, background, background, background);  
        end
        mx=res1.maxval;
        mi=res1.minval;
    elseif strcmp(topdown, 'SaliencyMap') == 1
        res1 = hdrvdp_saliencymap( ud.imgl, ud.imgr, 'sRGB-display', ppd, [], s, min, max, threshold_sliency);
        if res1.same == 0
            map = hdrvdp_visualize_ssim(background, res1.P_map, res1.patchPoints);
        else
            map = cat(3, background, background, background);  
        end
        mx=res1.maxval;
        mi=res1.minval;
    else
        T = double(ud.imgl)/2^8;
        R = double(ud.imgr)/2^8;
        res = hdrvdp( T, R, 'sRGB-display', ppd );   
        map= hdrvdp_visualize( 'pmap', res.P_map, {'context_image', T});  
    end
elseif strcmp(btupm, 'VDM') == 1
    if strcmp(topdown, 'SSIM') == 1
        res1 = VDM_IDX_SSIM( ud.imgl, ud.imgr, ppd, min, max, 'normal', s );
        if res1.same == 0
            map = vdm_visualize_ssim(background, res1.dmap, res1.patchPoints);
        else
            map = cat(3, background, background, background);  
        end
        mx=res1.maxval;
        mi=res1.minval;
    elseif strcmp(topdown, 'CW_SSIM') == 1
        res1 = VDM_IDX_SSIM( ud.imgl, ud.imgr, ppd, min, max, 'cw', s );
        if res1.same == 0
            map = vdm_visualize_ssim(background, res1.dmap, res1.patchPoints);
        else
            map =cat(3, background, background, background);  
        end
        mx=res1.maxval;
        mi=res1.minval;
    elseif strcmp(topdown, 'SaliencyMap') == 1
        res1 = VDM_IDX_Saliency( ud.imgl, ud.imgr, ppd, min, max, 'normal', s, threshold_sliency);
        if res1.same == 0
            map = vdm_visualize_ssim(background, res1.dmap, res1.patchPoints);
        else
            map = cat(3, background, background, background);  
        end
        mx=res1.maxval;
        mi=res1.minval;
    else
        res = VDM_IDX( ud.imgl, ud.imgr, ppd );  
        T = double(ud.imgl)/2^8;
        map= hdrvdp_visualize( 'pmap', res.dmap, {'context_image', T});  
        jnd = res.dprime;
    end
else
    if strcmp(topdown, 'SSIM') == 1
        T = rgb2gray(ud.imgl);
        R = rgb2gray(ud.imgr);
        %[map,mx,mi,jnd] = ssim_map(T,R, s, 'normal');
        [jnd,map] = ssim(T,R);
        map = double(cat(3, map, map, map));
    elseif strcmp(topdown, 'CW_SSIM') == 1
        T = rgb2gray(ud.imgl);
        R = rgb2gray(ud.imgr);
%         [map,mx,mi,jnd] = ssim_map(T,R, s, 'cw');
        [jnd,map]=cwssim_index(T,R,1,4,0,0);
%         map(isnan(map))=1;
        map = mat2gray(map);
        map = double(cat(3, map, map, map));
    else
        map= cat(3, background, background, background);  
    end
end
time = toc;
textLabel = sprintf('Min: %0.3f \nMax: %0.3f \nTime:%0.3f s\nJND:  %0.3f', mi,mx,time,jnd);
set(handles.txtime, 'String', textLabel);

ud.result = map;

axes(handles.axes6);
imshow(map,[])
axis off
axis image
% imshow(map);
% imshow(hdrvdp_visualize( 'pmap', res1.P_map, { 'context_image', R } ) );


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3


% --- Executes on button press in ckbhdr.
function ckbhdr_Callback(hObject, eventdata, handles)
% hObject    handle to ckbhdr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ckbhdr



function edtmin_Callback(hObject, eventdata, handles)
% hObject    handle to edtmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtmin as text
%        str2double(get(hObject,'String')) returns contents of edtmin as a double


% --- Executes during object creation, after setting all properties.
function edtmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtmax_Callback(hObject, eventdata, handles)
% hObject    handle to edtmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtmax as text
%        str2double(get(hObject,'String')) returns contents of edtmax as a double


% --- Executes during object creation, after setting all properties.
function edtmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtWNum_Callback(hObject, eventdata, handles)
% hObject    handle to edtWNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtWNum as text
%        str2double(get(hObject,'String')) returns contents of edtWNum as a double


% --- Executes during object creation, after setting all properties.
function edtWNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtWNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ppmlevel.
function ppmlevel_Callback(hObject, eventdata, handles)
% hObject    handle to ppmlevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ppmlevel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ppmlevel


% --- Executes during object creation, after setting all properties.
function ppmlevel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ppmlevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ppmangle.
function ppmangle_Callback(hObject, eventdata, handles)
% hObject    handle to ppmangle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ppmangle contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ppmangle


% --- Executes during object creation, after setting all properties.
function ppmangle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ppmangle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 global ud;
 levellist = get(handles.ppmlevel, 'String');
 levelval = get(handles.ppmlevel,'Value');
 level = str2num(levellist{levelval});
 anglelist = get(handles.ppmangle, 'String');
 angleval = get(handles.ppmangle, 'Value');
 angle = anglelist{angleval};
 I = ud.img;
 mask = cortexMask(level, angle, I);
 gray = rgb2gray(I);
 F = fftshift(fft2(gray));
%  F2 = log(abs(F));
 F=F.*mask;
 F2 = log(abs(F));
 
 filteredImage = ifft2(fftshift(F));
 amplitudeImage3 = abs(filteredImage);
  
 ud.result = mat2gray(amplitudeImage3);

 axes(handles.axes5);
 imshow(F2,[]);
 axes(handles.axes6);
 imshow(amplitudeImage3,[]);


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ud;
[filename, pathname] = uigetfile({'*.jpg;*.tif;*.png;*.gif;*.bmp','All Image Files';...
          '*.*','All Files' },'Image');
 I = imread([pathname,filename]);
 ud.img = I;
 axes(handles.axes1);
 imshow(I);


% --- Executes on button press in btnSave.
function btnSave_Callback(hObject, eventdata, handles)
% hObject    handle to btnSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ud;
[filename, pathname] = uiputfile({'*.png'},'Image');
imwrite(ud.result,[pathname, filename], 'png');



function editPPD_Callback(hObject, eventdata, handles)
% hObject    handle to editPPD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPPD as text
%        str2double(get(hObject,'String')) returns contents of editPPD as a double


% --- Executes during object creation, after setting all properties.
function editPPD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPPD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ud;
in = load('parameters.mat');
params = in.params;

% make Saliencymap
img_ref = initializeImage('', ud.imgl);

if (params.foaSize < 0)
  p = defaultSaliencyParams(img_ref.size);
  params.foaSize = p.foaSize;
end

[salMap,salData] = makeSaliencyMap(img_ref,params);
wta = initializeWTA(salMap,params);
[wta,thisWinner] = evolveWTA(wta);
wtaMap = emptyMap(img_ref.size(1:2),'Winner Take All');
wtaMap.data = imresize(wta.sm.V, img_ref.size(1:2),'bilinear');
bw = wtaMap.data/max(abs(wtaMap.data(:)));
ud.sliencyMap = bw;
val = get(handles.slider2,'Value');
threshold = val;
input = bw>threshold;
[B,L] = bwboundaries(input);
axes(handles.axes6);
imshow(ud.imgl,[]);
hold on
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1)
end

% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global ud;
val = get(handles.slider2,'Value');

bw = ud.sliencyMap;
threshold = val;
input = bw>threshold;
[B,L] = bwboundaries(input);
axes(handles.axes6);
imshow(ud.imgl);
hold on
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1)
end


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
