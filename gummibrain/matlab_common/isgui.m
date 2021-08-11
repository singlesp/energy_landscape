function y = isgui

y = ~(usejava('jvm') && ~feature('ShowFigureWindows'));