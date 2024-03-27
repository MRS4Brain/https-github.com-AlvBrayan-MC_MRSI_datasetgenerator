function hThresh = interactiveImageThresholder(handle)

% Do sanity checking on handles and take care of the zero-argument case.
if (nargin == 0)
    handle = get(0, 'CurrentFigure');
    if isempty(handle)
        error(message('images:common:notAFigureHandle', upper( mfilename )))
    end
end

iptcheckhandle(handle, {'figure', 'axes', 'image', 'uipanel'},...
    mfilename, 'H', 1);

[imageHandle, ~, figHandle] = imhandles(handle);

if (isempty(imageHandle))
    error(message('images:common:noImageInFigure'))
end

% Find and validate target image/axes.
imageHandle = imageHandle;
axHandle    = ancestor(imageHandle,'axes');
set(axHandle,'NextPlot','add');
set(figHandle,'Units','Normalized');
figpos      = get(figHandle,'Position');
hThresh     = figure('Name','INTERACTIVE THRESHOLDER','NumberTitle','off',...
                    'Units','Normalized','Position',[figpos(1) 0.2 figpos(3) 0.15],...
                    'DefaultuicontrolUnits','Normalized','Menubar','none'); 

im          = get(imageHandle,'CData'); 
maxval      = max(im(:));
minval      = min(im(:));
imrange     = maxval - minval;
set(hThresh,'DefaultuicontrolUnits','Normalized','DefaultuicontrolFontUnits','Points','DefaultuicontrolFontSize',14);
set(hThresh,'DefaultuipanelUnits','Normalized','DefaultuipanelFontSize',14);
iptwindowalign(figHandle,'left',hThresh,'left');
iptwindowalign(figHandle,'bottom',hThresh,'top');

h_pprimary = uipanel(hThresh,'Position',[0.01 0.01 .98 .98],'Title','Threshold Parameters');
h_sthresh  = uicontrol(h_pprimary,'Style','Slider','Position',[.05 .55 .6 .1]);
h_txtthreh = uicontrol(h_pprimary,'Style','Text','Position',[0.1 .7 .8 .15],'String','Threshold');
h_ethresh  = uicontrol(h_pprimary,'Style','Edit','Position',[0.7 .50 .2 .2],'Enable','on','BackgroundColor',[1 1 1]);
h_salpha   = uicontrol(h_pprimary,'Style','Slider','Position',[.05 .05 .6 .1]);
h_txtalpha = uicontrol(h_pprimary,'Style','Text','Position',[.1 .2 .8 .15],'String','Alpha');
h_bsavebin = uicontrol(h_pprimary,'Style','Pushbutton','Position',[0.7 0.025, .2 .2],'String','Save Binary');

h_overlay  = imshow(cat(3,zeros(size(im)),.7*ones(size(im)),zeros(size(im))),'Parent',axHandle);
set(h_overlay,'AlphaData',0.5);

initval    = min(im(:))-1;
set(h_sthresh,'Min',min(im(:))-1,'Max',max(im(:))+1,'Value',initval,'SliderStep',[1/imrange 10/imrange]);
set(h_ethresh,'String',num2str(initval));
set(h_salpha,'Min',0,'Max',1,'SliderStep',[1/100 1/10],'Value',0.5);

f_setthresh = @(hobject,eventdata)setthresh(hobject,eventdata);
lh_sthresh  = addlistener(h_sthresh,'ContinuousValueChange',f_setthresh);
setappdata(hThresh,'ThreshSliderListener',lh_sthresh);

f_setalpha  = @(hobject,eventdata)setalpha(hobject,eventdata);
lh_salpha   = addlistener(h_salpha,'ContinuousValueChange',f_setalpha);
setappdata(hThresh,'AlphaSliderListener',lh_salpha);

f_manthresh = @(hobject,eventdata)manthresh(hobject,eventdata);
set(h_ethresh,'Callback',f_manthresh);

f_quitThresholder = @(hobject,eventdata)quitThresholder(hobject,eventdata);
set(hThresh,'DeleteFcn',f_quitThresholder);

f_savebinary = @(hobject,eventdata)savebinary(hobject,eventdata);
set(h_bsavebin,'Callback',f_savebinary);

    function setthresh(varargin)
        h   = varargin{1};
        val = get(h,'Value');
        set(h_ethresh,'String',num2str(val));
        updatethreshimg;
    end

    function setalpha(varargin)
        h   = varargin{1};
        val = get(h,'Value');
        thresh = get(h_sthresh,'Value');
        updatethreshimg;
    end

    function manthresh(varargin)
        h    = varargin{1};
        val  = get(h,'String');
        nval = str2num(val);
        if isnumeric(nval);
            if nval > maxval
                newval = maxval;
            elseif nval < minval
                newval = minval;
            else
                newval = nval;
            end
            set(h_sthresh,'Value',newval);
            set(h_ethresh,'String',num2str(newval));
            updatethreshimg;
        end
    end

    function updatethreshimg()
        alphaval  = get(h_salpha,'Value');
        threshval = get(h_sthresh,'Value');
        set(h_overlay,'AlphaData',alphaval*(im > threshval));
    end
        
    function savebinary(varargin)
        threshval = get(h_sthresh,'Value');
        binimg    = im > threshval;
        assignin('base','imThresh',binimg);
    end

    function quitThresholder(varargin)
        h = varargin{1};
        if ishandle(h_overlay)
            delete(h_overlay);
        end
        delete(h);
    end
end





