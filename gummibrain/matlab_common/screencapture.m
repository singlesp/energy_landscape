function imageData = screencapture(varargin)
% screencapture - get a screen-capture of a figure frame, component handle, or screen area rectangle
%
% ScreenCapture gets a screen-capture of any Matlab GUI handle (including desktop, 
% figure, axes or uicontrol), or a specified area rectangle located relative to the
% specified handle. Screen area capture is possible by specifying the root (desktop)
% handle (=0). The output can be either to an image file or to a Matlab matrix (useful
% for displaying via imshow() or for further processing). This utility also enables
% adding a toolbar button for easy interactive screen-capture.
%
% Syntax:
%    imageData = screencapture(handle, position, filename,
%    'PropName',PropValue, ...)
%
% Input Parameters:
%    handle      - optional handle to be used for screen-capture origin.
%                    If empty/unsupplied then current figure (gcf) will be used.
%    position    - optional position array in pixels: [x,y,width,height].
%                    If empty/unsupplied then the handle's position vector will be used.
%                    If both handle and position are empty/unsupplied then the position
%                      will be retrieved via interactive mouse-selection.
%    filename    - optional filename for storing the screen-capture.
%                    If empty/unsupplied then no output to file will be done.
%                    The file format will be determined from the filename (JPG/PNG/...).
%                    Supported formats are those supported by the imwrite function.
%                    If neither filename nor imageData were specified, the user will be
%                      asked to interactively specify the output filename.
%    'PropName',PropValue - 
%                  optional list of property pairs (e.g., 'filename','myImage.png','pos',[10,20,30,40],'handle',gca)
%                  PropNames may be abbreviated and are case-insensitive.
%                  PropNames may also be given in whichever order.
%                  Supported PropNames are:
%                    - 'handle'    (default: gcf handle)
%                    - 'position'  (default: gcf position array)
%                    - 'filename'  (default: '')
%                    - 'toolbar'   (figure handle; default: gcf)
%                         this adds a screen-capture button to the figure's toolbar
%                         If this parameter is specified, then no screen-capture
%                           will take place and the returned imageData will be [].
%
% Output parameters:
%    imageData   - image data in a format acceptable by the imshow function
%                    If neither filename nor imageData were specified, the user will be
%                      asked to interactively specify the output filename.
%
% Examples:
%    imageData = screencapture;  % interactively select screen-capture rectangle
%    imageData = screencapture(hListbox);  % capture image of a uicontrol
%    imageData = screencapture(0,  [20,30,40,50]);  % select a small desktop region
%    imageData = screencapture(gcf,[20,30,40,50]);  % select a small figure region
%    imageData = screencapture(gca,[10,20,30,40]);  % select a small axes region
%      imshow(imageData);  % display the captured image in a matlab figure
%      imwrite(imageData,'myImage.png');  % save the captured image to file
%    screencapture(gcf,[],'myFigure.jpg');  % capture the entire figure into file
%    screencapture('handle',gcf,'filename','myFigure.jpg');  % same as previous
%    screencapture('toolbar',gcf);  % adds a screen-capture button to gcf's toolbar
%    screencapture('toolbar',[],'file','sc.bmp'); % same with default output filename
%
% Bugs and suggestions:
%    Please send to Yair Altman (altmany at gmail dot com)
%
% See also:
%    imshow, imwrite, print
%
% Release history:
%    1.1 2009-06-03: Handle missing output format; performance boost (thanks to Urs); fix minor root-handle bug; added toolbar button option
%    1.0 2009-06-02: First version posted on <a href="http://www.mathworks.com/matlabcentral/fileexchange/authors/27420">MathWorks File Exchange</a>

% License to use and modify this code is granted freely to all interested, as long as the original author is
% referenced and attributed as such. The original author maintains the right to be solely associated with this work.

% Programmed and Copyright by Yair M. Altman: altmany(at)gmail.com
% $Revision: 1.1 $  $Date: 2009/06/03 19:25:29 $

  % Ensure that java awt is enabled...
  if ~usejava('awt')
      error('YMA:screencapture:NeedAwt','ScreenCapture requires Java to run.');
  end

  % Process optional arguments
  paramsStruct = processArgs(varargin{:});

  % If toolbar button requested, add it and exit
  if ~isempty(paramsStruct.toolbar)
      addToolbarButton(paramsStruct);
      % Exit immediately (do NOT take a screen-capture)
      if nargout,  imageData = [];  end
      return;
  end

  % Capture the requested screen rectangle using java.awt.Robot
  imgData = getScreenCaptureImageData(paramsStruct.position);

  % Save image data in filename, if supplied
  if ~isempty(paramsStruct.filename)
      imwrite(imgData,paramsStruct.filename);
  end

  % Return image raster data to user, if requested
  if nargout
      imageData = imgData;

  % If neither output formats was specified (neighter filename nor output data)
  elseif isempty(paramsStruct.filename)
      % Ask the user to specify a filename
      %error('YMA:screencapture:noOutput','No output specified for ScreenCapture: specify the output filename and/or output data');
      [filename,pathname] = uiputfile('*.*','Save screen capture as');
      if ~isequal(filename,0) & ~isequal(pathname,0)  %#ok Matlab6 compatibility
          imwrite(imgData,fullfile(pathname,filename));
      end
  end
  return;  % debug breakpoint

%% Process optional arguments
function paramsStruct = processArgs(varargin)

    % Get the properties in either direct or P-V format
    [regParams, pvPairs] = parseparams(varargin);

    % Now process the optional P-V params
    try
        % Initialize
        paramName = [];
        paramsStruct = [];
        paramsStruct.handle = [];
        paramsStruct.position = [];
        paramsStruct.filename = '';
        paramsStruct.toolbar = [];

        % Parse the regular (non-named) params in recption order
        if ~isempty(regParams) && (isempty(regParams{1}) || ishandle(regParams{1}(1)))
            paramsStruct.handle = regParams{1};
            regParams(1) = [];
        end
        if ~isempty(regParams) && isnumeric(regParams{1}) && (length(regParams{1}) == 4)
            paramsStruct.position = regParams{1};
            regParams(1) = [];
        end
        if ~isempty(regParams)
            paramsStruct.filename = regParams{1};
        end

        % Parse the optional param PV pairs
        supportedArgs = {'handle','position','filename','toolbar'};
        while ~isempty(pvPairs)

            % Ensure basic format is valid
            paramName = '';
            if ~ischar(pvPairs{1})
                error('YMA:screencapture:invalidProperty','Invalid property passed to ScreenCapture');
            elseif length(pvPairs) == 1
                if isempty(paramsStruct.filename)
                    paramsStruct.filename = pvPairs{1};
                    break;
                else
                    error('YMA:screencapture:noPropertyValue',['No value specified for property ''' pvPairs{1} '''']);
                end
            end

            % Process parameter values
            paramName  = pvPairs{1};
            paramValue = pvPairs{2};
            pvPairs(1:2) = [];
            idx = find(strncmpi(paramName,supportedArgs,length(paramName)));
            if ~isempty(idx)
                paramsStruct.(lower(supportedArgs{idx(1)})) = paramValue;

                % If 'toolbar' param specified, then it cannot be left empty - use gcf
                if strncmpi(paramName,'toolbar',length(paramName)) & isempty(paramsStruct.toolbar)  %#ok ML6
                    paramsStruct.toolbar = getCurrentFig;
                end

            elseif isempty(paramsStruct.filename)
                paramsStruct.filename = paramName;
                pvPairs = {paramValue, pvPairs{:}};

            else
                supportedArgsStr = sprintf('''%s'',',supportedArgs{:});
                error('YMA:screencapture:invalidProperty','%s \n%s', ...
                      'Invalid property passed to ScreenCapture', ...
                      ['Supported property names are: ' supportedArgsStr(1:end-1)]);
            end
        end  % loop pvPairs

    catch
        if ~isempty(paramName),  paramName = [' ''' paramName ''''];  end
        error('YMA:screencapture:invalidProperty','Error setting ScreenCapture property %s:\n%s',paramName,lasterr);
    end

    % No handle specified
    if isempty(paramsStruct.handle)

        % Set default handle, if not supplied
        if ~isempty(paramsStruct.position)
            paramsStruct.handle = getCurrentFig;
        else
            [paramsStruct.handle, paramsStruct.position] = getInteractivePosition;
        end

    % Handle was supplied - ensure it is a valid handle
    elseif ~ishandle(paramsStruct.handle)
        error('YMA:screencapture:invalidHandle','Invalid handle passed to ScreenCapture');

    % Handle was supplied but position was not, so use the handle's position
    elseif isempty(paramsStruct.position)
        paramsStruct.position = getPixelPos(paramsStruct.handle);
        paramsStruct.position(1:2) = 0;

    % Both handle & position were supplied - ensure a valid pixel position vector
    elseif ~isnumeric(paramsStruct.position) || (length(paramsStruct.position) ~= 4)
        error('YMA:screencapture:invalidPosition','Invalid position vector passed to ScreenCapture: \nMust be a [x,y,w,h] numeric pixel array');
    end

    % Convert position from handle-relative to desktop Java-based pixels
    paramsStruct.position = convertPos(paramsStruct.handle, paramsStruct.position);
%end  % processArgs

%% Convert position from handle-relative to desktop Java-based pixels
function position = convertPos(hParent, position)

    try
        % Get the screen-size for later use
        screenSize = get(0,'ScreenSize');

        % First get the parent handle's desktop-based Matlab pixel position
        if ~isa(handle(hParent),'figure')
            parentPos = getPixelPos(hParent, true);
            hParent = ancestor(hParent,'figure');
        else
            % Compensate 4px figure boundaries = difference betweeen OuterPosition, Position
            parentPos = [1-4,1-4,0,0];
        end
        figurePos = getPixelPos(hParent);
        desktopPos = figurePos + parentPos;

        % Now convert to Java-based pixels based on screen size
        % TODO: handle multiple monitors
        javaY = screenSize(4) - desktopPos(2) - position(2) - position(4) - 2;
        position = [desktopPos(1)+position(1)+2 javaY position(3:4)];

        % Ensure the figure is at the front so it can be screen-captured
        figure(hParent);
        drawnow;
        pause(0.02);
    catch
        % Maybe root/desktop handle (root does not have a 'Position' prop so getPixelPos croaks
        if double(hParent)==0  % =root/desktop handle
            javaY = screenSize(4) - position(2) - position(4) - 2;
            position = [position(1)-1 javaY position(3:4)];
        end
    end
%end  % convertPos

%% Interactively get the requested capture rectangle
function [hFig, positionRect] = getInteractivePosition
    uiwait(msgbox('Mouse-click in any figure and drag a bounding rectangle for screen-capture','ScreenCapture'));
    k = waitforbuttonpress;  %#ok k is unused
    hFig = getCurrentFig;
    %p1 = get(hFig,'CurrentPoint');
    positionRect = rbbox; 
    %p2 = get(hFig,'CurrentPoint');
%end  % getInteractivePosition

%% Get current figure (even if its handle is hidden)
function hFig = getCurrentFig
    oldState = get(0,'showHiddenHandles');
    set(0,'showHiddenHandles','on');
    hFig = get(0,'CurrentFigure');
    set(0,'showHiddenHandles',oldState);
%end  % getCurrentFig

%% Get ancestor figure - used for old Matlab versions that don't have a built-in ancestor()
function hObj = ancestor(hObj,type)
    if ~isempty(hObj) & ishandle(hObj)  %#ok for Matlab 6 compatibility
        try
            hObj = get(hObj,'Ancestor');
        catch
            % never mind...
        end
        try
            %if ~isa(handle(hObj),type)  % this is best but always returns 0 in Matlab 6!
            %if ~isprop(hObj,'type') | ~strcmpi(get(hObj,'type'),type)  % no isprop() in ML6!
            objType=''; try objType=get(hObj,'type'); catch, end  %#ok
            if ~strcmpi(objType,type)
                try
                    parent = get(handle(hObj),'parent');
                catch
                    parent = hObj.getParent;  % some objs have no 'Parent' prop, just this method...
                end
                if ~isempty(parent)  % empty parent means root ancestor, so exit
                    hObj = ancestor(parent,type);
                end
            end
        catch
            % never mind...
        end
    end
%end  % ancestor

%% Get position of an HG object in specified units
function pos = getPos(hObj,field,units)
    % Matlab 6 did not have hgconvertunits so use the old way...
    oldUnits = get(hObj,'units');
    if strcmpi(oldUnits,units)  % don't modify units unless we must!
        pos = get(hObj,field);
    else
        set(hObj,'units',units);
        pos = get(hObj,field);
        set(hObj,'units',oldUnits);
    end
%end  % getPos

%% Get pixel position of an HG object - for Matlab 6 compatibility
function pos = getPixelPos(hObj,varargin)
    try
        if ~isa(handle(hObj),'figure')
            % getpixelposition is unvectorized unfortunately!
            pos = getpixelposition(hObj,varargin{:});
        else
            pos = getPos(hObj,'OuterPosition','pixels');
        end
    catch
        try
            % Matlab 6 did not have getpixelposition nor hgconvertunits so use the old way...
            pos = getPos(hObj,'Position','pixels');
        catch
            % Maybe the handle does not have a 'Position' prop (e.g., text/line/plot) - use its parent
            pos = getPixelPos(get(hObj,'parent'),varargin{:});
        end
    end
%end  % getPixelPos

%% Adds a ScreenCapture toolbar button
function addToolbarButton(paramsStruct)
    % Ensure we have a valid toolbar handle
    hToolbar = findall(paramsStruct.toolbar,'type','uitoolbar');
    if isempty(hToolbar)
        error('YMA:screencapture:noToolbar','the ''Toolbar'' parameter must contain a figure handle possessing a valid toolbar');
    end
    hToolbar = hToolbar(1);  % just in case there are several toolbars... - use only the first

    % Prepare the camera icon
    icon = ['3333333333333333'; ...
            '3333333333333333'; ...
            '3333300000333333'; ...
            '3333065556033333'; ...
            '3000000000000033'; ...
            '3022222222222033'; ...
            '3022220002222033'; ...
            '3022203110222033'; ...
            '3022201110222033'; ...
            '3022204440222033'; ...
            '3022220002222033'; ...
            '3022222222222033'; ...
            '3000000000000033'; ...
            '3333333333333333'; ...
            '3333333333333333'; ...
            '3333333333333333'];
    cm = [   0      0      0; ...  % black
             0   0.60      1; ...  % light blue
          0.53   0.53   0.53; ...  % light gray
           NaN    NaN    NaN; ...  % transparent
             0   0.73      0; ...  % light green
          0.27   0.27   0.27; ...  % gray
          0.13   0.13   0.13];     % dark gray
    cdata = ind2rgb(uint8(icon-'0'),cm);

    % If the button does not already exit
    hButton = findall(hToolbar,'Tag','ScreenCaptureButton');
    tooltip = 'Screen capture';
    if ~isempty(paramsStruct.filename)
        tooltip = [tooltip ' to ' paramsStruct.filename];
    end
    if isempty(hButton)
        % Add the button with the icon to the figure's toolbar
        hButton = uipushtool(hToolbar, 'CData',cdata, 'Tag','ScreenCaptureButton', 'TooltipString',tooltip, 'ClickedCallback',['screencapture(''' paramsStruct.filename ''')']);  %#ok unused
    else
        % Otherwise, simply update the existing button
        set(hButton, 'CData',cdata, 'Tag','ScreenCaptureButton', 'TooltipString',tooltip, 'ClickedCallback',['screencapture(''' paramsStruct.filename ''')']);
    end
%end  % addToolbarButton

%% Java-get the actual screen-capture image data
function imgData = getScreenCaptureImageData(positionRect)
    rect = java.awt.Rectangle(positionRect(1), positionRect(2), positionRect(3), positionRect(4));
    robot = java.awt.Robot;
    jImage = robot.createScreenCapture(rect);

    % Convert the resulting Java image to a Matlab image
    % Adapted for a much-improved performance from:
    % http://www.mathworks.com/support/solutions/data/1-2WPAYR.html
    h = jImage.getHeight;
    w = jImage.getWidth;
    imgData = zeros([h,w,3],'uint8');
    %pixelsData = uint8(jImage.getData.getPixels(0,0,w,h,[]));
    %for i = 1 : h
    %    base = (i-1)*w*3+1;
    %    imgData(i,1:w,:) = deal(reshape(pixelsData(base:(base+3*w-1)),3,w)');
    %end

    % Performance further improved based on feedback from Urs:
    pixelsData = reshape(typecast(jImage.getData.getDataStorage,'uint32'),w,h).';
    imgData(:,:,3) = bitshift(bitand(pixelsData,256^1-1),-8*0);
    imgData(:,:,2) = bitshift(bitand(pixelsData,256^2-1),-8*1);
    imgData(:,:,1) = bitshift(bitand(pixelsData,256^3-1),-8*2);
%end  % getInteractivePosition



%%%%%%%%%%%%%%%%%%%%%%%%%% TODO %%%%%%%%%%%%%%%%%%%%%%%%%
% - Fix: handle multiple monitors
% - Enh: enable interactive desktop rbbox that does not necesarily start within a figure
% - Enh: capture current object (uicontrol/axes/figure) if w=h=0
