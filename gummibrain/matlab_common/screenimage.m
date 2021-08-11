function img = screenimage(frame)

if isscalar(frame) && ishandle(frame) && (frame > 0)
    inputType = get(frame,'type');
else
    error('Invalid input argument.  frame must be a handle to a figure or axis.');
end

if(strcmpi(inputType,'axes'))
    fig = get(frame,'parent');
else
    fig = frame;
end

set(fig,'PaperPositionMode','auto');

renderer = get(fig,'renderer');

if strcmp(renderer,'painters')
    renderer = 'opengl';
end

pixelsperinch = get(0,'screenpixelsperInch');

img = hardcopy(frame, ['-d' renderer], ['-r' num2str(round(pixelsperinch))]);