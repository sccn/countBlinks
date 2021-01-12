% countBlinks() - A tool to support to annotate blinks manually using VEOG and HEOG 
%                 IC/channel time-series. Developed for a collabolation with UCLA Sandra Loo.

% History
% 06/28/2018 ver 0.40 by Makoto. VEOG and HEOG can be plotted. Function structure optimized by eliminating redundant code for plotting (by separating updateTimeSeriesPlots() as a subfunction)
% 06/27/2018 ver 0.31 by Makoto. Scalp topo not mean but peak latency plotted.
% 06/20/2018 ver 0.30 by Makoto. Scalp topo added.
% 06/15/2018 ver 0.20 by Makoto. Updated. 
% 05/25/2018 ver 0.10 by Makoto. Created.
%
% Copyright (C) 2018, Makoto Miyakoshi. SCCN, INC, UCSD.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1.07  USA

function varargout = countBlinks(varargin)
% COUNTBLINKS MATLAB code for countBlinks.fig
%      COUNTBLINKS, by itself, creates a new COUNTBLINKS or raises the existing
%      singleton*.
%
%      H = COUNTBLINKS returns the handle to a new COUNTBLINKS or the handle to
%      the existing singleton*.
%
%      COUNTBLINKS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COUNTBLINKS.M with the given input arguments.
%
%      COUNTBLINKS('Property','Value',...) creates a new COUNTBLINKS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before countBlinks_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to countBlinks_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help countBlinks

% Last Modified by GUIDE v2.5 28-Jun-2018 19:10:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @countBlinks_OpeningFcn, ...
                   'gui_OutputFcn',  @countBlinks_OutputFcn, ...
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


% --- Executes just before countBlinks is made visible.
function countBlinks_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to countBlinks (see VARARGIN)

% Choose default command line output for countBlinks
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes countBlinks wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = countBlinks_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on button press in markAsBlinkButton.
function markAsBlinkButton_Callback(hObject, eventdata, handles)
% hObject    handle to markAsBlinkButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load data.
EEG                       = evalin('base', 'EEG');
extendedVeogTimeSeries    = evalin('base', 'extendedVeogTimeSeries');
extendedHeogTimeSeries    = evalin('base', 'extendedHeogTimeSeries');
extendedTimes             = evalin('base', 'extendedTimes');
blinkPeakLatencyFrameList = evalin('base', 'blinkPeakLatencyFrameList');
visMaskWidthSec           = evalin('base', 'visMaskWidthSec');

% Detect the global peak.
[~, peakFrameIdx] = max(extendedVeogTimeSeries(1,:));

% Calculat the indices for the current window.
maskWidthSec   = str2num(get(handles.replaceWindowLengthEdit, 'String'));
maskWidthFrame = EEG.srate*maskWidthSec;
currentWindowFrameIdx = peakFrameIdx-ceil(maskWidthFrame/2):peakFrameIdx+ceil(maskWidthFrame/2);

% Accept the current peak as blink and store the current peak latency.
blinkPeakLatencyFrameList = [blinkPeakLatencyFrameList peakFrameIdx];

% Mask the current time series.
extendedVeogTimeSeries(1,currentWindowFrameIdx) = NaN;

% Report the number of blinks.
disp(sprintf('Total number of blinks: %.0f', length(blinkPeakLatencyFrameList)))

% Save data.
assignin('base', 'extendedVeogTimeSeries',    extendedVeogTimeSeries)
assignin('base', 'extendedTimes',             extendedTimes)
assignin('base', 'blinkPeakLatencyFrameList', blinkPeakLatencyFrameList)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Update the time-series plots. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
updateTimeSeriesPlots(hObject, eventdata, handles)



% --- Executes on button press in skipButton.
function skipButton_Callback(hObject, eventdata, handles)
% hObject    handle to skipButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load data.
EEG                       = evalin('base', 'EEG');
extendedVeogTimeSeries    = evalin('base', 'extendedVeogTimeSeries');
extendedHeogTimeSeries    = evalin('base', 'extendedHeogTimeSeries');
extendedTimes             = evalin('base', 'extendedTimes');
blinkPeakLatencyFrameList = evalin('base', 'blinkPeakLatencyFrameList');
visMaskWidthSec           = evalin('base', 'visMaskWidthSec');

% Detect the global peak.
[~, peakFrameIdx] = max(extendedVeogTimeSeries(1,:));

% Calculat the indices for the current window.
maskWidthSec   = str2num(get(handles.replaceWindowLengthEdit, 'String'));
maskWidthFrame = EEG.srate*maskWidthSec;
currentWindowFrameIdx = peakFrameIdx-ceil(maskWidthFrame/2):peakFrameIdx+ceil(maskWidthFrame/2);

% Accept the current peak as blink and store the current peak latency.
%blinkPeakLatencyFrameList = [blinkPeakLatencyFrameList peakFrameIdx]; % Because skipping.

% Mask the current time series.
extendedVeogTimeSeries(1,currentWindowFrameIdx) = NaN;

% Report the number of blinks.
%disp(sprintf('Total number of blinks: %.0f',%length(blinkPeakLatencyFrameList))) % Because skipping.

% Save data.
assignin('base', 'extendedVeogTimeSeries',    extendedVeogTimeSeries)
assignin('base', 'extendedTimes',             extendedTimes)
%assignin('base', 'blinkPeakLatencyFrameList', blinkPeakLatencyFrameList) % Because skipping.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Update the time-series plots. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
updateTimeSeriesPlots(hObject, eventdata, handles)



function replaceWindowLengthEdit_Callback(hObject, eventdata, handles)
% hObject    handle to replaceWindowLengthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of replaceWindowLengthEdit as text
%        str2double(get(hObject,'String')) returns contents of replaceWindowLengthEdit as a double

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Update the time-series plots. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
updateTimeSeriesPlots(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function replaceWindowLengthEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to replaceWindowLengthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in veogPopupmenu.
function veogPopupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to veogPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns veogPopupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from veogPopupmenu


% --- Executes during object creation, after setting all properties.
function veogPopupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to veogPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function veogEdit_Callback(hObject, eventdata, handles)
% hObject    handle to veogEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of veogEdit as text
%        str2double(get(hObject,'String')) returns contents of veogEdit as a double


% --- Executes during object creation, after setting all properties.
function veogEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to veogEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in startPushbutton.
function startPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to startPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load EEG.
EEG = evalin('base', 'EEG');

% Obtain channel index.
veogIdx = str2num(get(handles.veogEdit, 'string'));
heogIdx = str2num(get(handles.heogEdit, 'string'));

% Obtain the current time-series data for eye blinks.
if get(handles.veogPopupmenu, 'value') == 1
    veogTimeSeries = EEG.icaact(veogIdx,:);
else
    veogTimeSeries = EEG.data(veogIdx,:);
end

% Obtain the current time-series data for horizontal eye movements.
if get(handles.heogEdit, 'value') == 1
    heogTimeSeries = EEG.icaact(heogIdx,:);
else
    heogTimeSeries = EEG.data(heogIdx,:);
end

% Crop the visualization window.
visMaskWidthSec   = 5;
visMaskWidthFrame = EEG.srate*visMaskWidthSec;

% Zero-pad the both ends of the time series.
extendedVeogTimeSeries = [nan(1,visMaskWidthFrame) veogTimeSeries nan(1,visMaskWidthFrame)];
extendedHeogTimeSeries = [nan(1,visMaskWidthFrame) heogTimeSeries nan(1,visMaskWidthFrame)];
timeInterval       = mean(diff(EEG.times));
extendedTimeLeft   = visMaskWidthFrame*-timeInterval:timeInterval:-timeInterval;
extendedTimeRight  = EEG.times(end)+timeInterval:timeInterval:EEG.times(end)+visMaskWidthFrame*timeInterval;
extendedTimes      = [extendedTimeLeft EEG.times extendedTimeRight];

% Duplicate extendedVeogTimeSeries
extendedVeogTimeSeries = repmat(extendedVeogTimeSeries, [2 1]);

% Save the current peak frame and window frames.
blinkPeakLatencyFrameList = []; 
assignin('base', 'extendedVeogTimeSeries',    extendedVeogTimeSeries)
assignin('base', 'extendedHeogTimeSeries',    extendedHeogTimeSeries)
assignin('base', 'extendedTimes',             extendedTimes)
assignin('base', 'blinkPeakLatencyFrameList', blinkPeakLatencyFrameList)
assignin('base', 'visMaskWidthSec',           visMaskWidthSec)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Update the time-series plots. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
updateTimeSeriesPlots(hObject, eventdata, handles)


  
% --- Executes on button press in replotPushbutton.
function replotPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to replotPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Update the time-series plots. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
updateTimeSeriesPlots(hObject, eventdata, handles)



% --- Executes on button press in finishButton.
function finishButton_Callback(hObject, eventdata, handles)
% hObject    handle to finishButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load data.
EEG                       = evalin('base', 'EEG');
blinkPeakLatencyFrameList = evalin('base', 'blinkPeakLatencyFrameList');
blinkPeakLatencyFrameList = sort(blinkPeakLatencyFrameList);
visMaskWidthSec           = evalin('base', 'visMaskWidthSec');

% Apply the -5 sec here.
blinkPeakLatencyFrameList = blinkPeakLatencyFrameList-EEG.srate*visMaskWidthSec;

% Populate blink markers.
tmpEvent = EEG.event;
for blinkIdx = 1:length(blinkPeakLatencyFrameList)
    tmpEvent(end+1).latency = blinkPeakLatencyFrameList(blinkIdx);
    tmpEvent(end).type      = 'blink';
end
tmpEventLatency = cell2mat({tmpEvent.latency}');
[~, sortIdx]    = sort(tmpEventLatency);
tmpEvent        = tmpEvent(sortIdx);
EEG.event       = tmpEvent;
EEG             = eeg_checkset(EEG, 'makeur');

% Update EEG.
assignin('base', 'EEG', EEG)
eeglab redraw

disp('EEG structure updated with blink events.')


% --- Executes on button press in halfPushbutton.
function halfPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to halfPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Obtain the current window width.
currentValue = str2num(get(handles.replaceWindowLengthEdit, 'String'));

% Apply the calculation.
updatedValue = round((currentValue/2)*1000)/1000;

% Update handles structure.
set(handles.replaceWindowLengthEdit, 'String', num2str(updatedValue));
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Update the time-series plots. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
updateTimeSeriesPlots(hObject, eventdata, handles)



% --- Executes on button press in doublePushbutton.
function doublePushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to doublePushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Obtain the current window width.
currentValue = str2num(get(handles.replaceWindowLengthEdit, 'String'));

% Apply the calculation.
updatedValue = round((currentValue*2)*1000)/1000;

% Update handles structure.
set(handles.replaceWindowLengthEdit, 'String', num2str(updatedValue));
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Update the time-series plots. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
updateTimeSeriesPlots(hObject, eventdata, handles)


% --- Executes on selection change in heogPopupmenu.
function heogPopupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to heogPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns heogPopupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from heogPopupmenu


% --- Executes during object creation, after setting all properties.
function heogPopupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to heogPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function heogEdit_Callback(hObject, eventdata, handles)
% hObject    handle to heogEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of heogEdit as text
%        str2double(get(hObject,'String')) returns contents of heogEdit as a double


% --- Executes during object creation, after setting all properties.
function heogEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to heogEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function updateTimeSeriesPlots(hObject, eventdata, handles)

% Load data.
EEG                       = evalin('base', 'EEG');
extendedVeogTimeSeries    = evalin('base', 'extendedVeogTimeSeries');
extendedHeogTimeSeries    = evalin('base', 'extendedHeogTimeSeries');
extendedTimes             = evalin('base', 'extendedTimes');
blinkPeakLatencyFrameList = evalin('base', 'blinkPeakLatencyFrameList');
visMaskWidthSec           = evalin('base', 'visMaskWidthSec');

% Detect the global peak.
[~, peakFrameIdx] = max(extendedVeogTimeSeries(1,:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Update the full-length time series plot. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cla(handles.wholeTimeSeriesAxes)
maskedVeogTimeSeries   = extendedVeogTimeSeries(1,EEG.srate*visMaskWidthSec+1:end-EEG.srate*visMaskWidthSec);
originalVeogTimeSeries = extendedVeogTimeSeries(2,EEG.srate*visMaskWidthSec+1:end-EEG.srate*visMaskWidthSec);
plot(handles.wholeTimeSeriesAxes, EEG.times/1000, originalVeogTimeSeries, 'color', [1 0.7 0.7]);
hold(handles.wholeTimeSeriesAxes, 'on')
plot(handles.wholeTimeSeriesAxes, EEG.times/1000, maskedVeogTimeSeries)
hold(handles.wholeTimeSeriesAxes, 'off')
xlim(handles.wholeTimeSeriesAxes, [EEG.times(1)/1000 EEG.times(end)/1000])

% If this is the first time to plot the whole VEOG time series.
if any(blinkPeakLatencyFrameList)
    peakValues = extendedVeogTimeSeries(2,blinkPeakLatencyFrameList);
    maxValue   = median(peakValues);
else
    maxValue = max(maskedVeogTimeSeries);
end
ylim(handles.wholeTimeSeriesAxes, [-maxValue maxValue])
xlabel(handles.wholeTimeSeriesAxes, 'Time (s)')
ylabel(handles.wholeTimeSeriesAxes, '(RMS) Amplitude (uV)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Update the zoom-window time series plot. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cla(handles.currentWindowAxes)

% Crop the data for visualization.
visMaskWidthFrame   = visMaskWidthSec*EEG.srate;
visWindowFrameIdx   = peakFrameIdx-ceil(visMaskWidthFrame/2):peakFrameIdx+ceil(visMaskWidthFrame/2);
visVeog             = extendedVeogTimeSeries(2, visWindowFrameIdx);
visHeog             = extendedHeogTimeSeries(1, visWindowFrameIdx);
visDataNoRejMaskIdx = find(~isnan(extendedVeogTimeSeries(1, visWindowFrameIdx)));
visDataRejMaskIdx   = find(isnan(extendedVeogTimeSeries(1, visWindowFrameIdx)));

% Crop the masking window.
maskWidthSec   = str2num(get(handles.replaceWindowLengthEdit, 'String'));
maskWidthFrame = EEG.srate*maskWidthSec;
currentWindowFrameIdx = peakFrameIdx-ceil(maskWidthFrame/2):peakFrameIdx+ceil(maskWidthFrame/2);
maskData       = extendedVeogTimeSeries(2, currentWindowFrameIdx);
maskTime       = extendedTimes(currentWindowFrameIdx)/1000;

% Prepare NaN-masked waveforms to show.
visVeogNoRej = visVeog; 
visVeogRej   = visVeog;
visVeogNoRej(visDataRejMaskIdx) = NaN;
visVeogRej(visDataNoRejMaskIdx) = NaN;

% Plot the window.
plotTimes = extendedTimes(visWindowFrameIdx)/1000;
hold(handles.currentWindowAxes, 'on')
if get(handles.plotHeogCheckbox, 'value') == 1    
    plot(handles.currentWindowAxes, plotTimes, visHeog,    'color', [0.76 0.76 0.76])
end
plot(handles.currentWindowAxes, plotTimes, visVeogRej, 'r')
plot(handles.currentWindowAxes, plotTimes, visVeogNoRej)
xlim(handles.currentWindowAxes, [plotTimes(1) plotTimes(end)])
ylim(handles.currentWindowAxes, [-maxValue maxValue])
xlabel('Time (s)')
ylabel('(RMS) Amplitude (uV)')

% Plot the correction mask.
axes(handles.currentWindowAxes)
fillHarmony = fill([maskTime maskTime(end:-1:1)], [ones(1,length(maskTime))*maxValue -ones(1,length(maskTime))*maxValue], [1 0.85 0.85]);
set(fillHarmony, 'lineStyle', 'none')
uistack(fillHarmony, 'bottom')
hold(handles.currentWindowAxes, 'off')

% Plot the mask on the whole time-series
axes(handles.wholeTimeSeriesAxes)
hold(handles.wholeTimeSeriesAxes, 'on')
maxValue = get(handles.wholeTimeSeriesAxes, 'YLim');
fillHarmony = fill([maskTime maskTime(end:-1:1)], [ones(1,length(maskTime))*maxValue(1) ones(1,length(maskTime))*maxValue(2)], [0.8 1 0.8]);  
set(fillHarmony, 'lineStyle', 'none')
uistack(fillHarmony, 'bottom')
hold(handles.wholeTimeSeriesAxes, 'off')

% Show scalp topo of the peak frame.
axes(handles.scalpTopoAxes)
cla
topoplotCustom(EEG.data(:, peakFrameIdx-visMaskWidthFrame), EEG.chanlocs);
title('Scalp topography at peak latency');
set(gcf, 'color', [0.66 0.76 1]);


% --- Executes on button press in plotHeogCheckbox.
function plotHeogCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to plotHeogCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotHeogCheckbox

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Update the time-series plots. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
updateTimeSeriesPlots(hObject, eventdata, handles)
