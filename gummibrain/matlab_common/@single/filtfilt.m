function y = filtfilt(b,a,x)
%w = warning('off','MATLAB:dispatcher:nameConflict');

y = cast(filtfilt(b,a,double(x)),class(x));

%warning(w.state,'MATLAB:dispatcher:nameConflict');
