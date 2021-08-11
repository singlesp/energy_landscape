function F = dft_fftreplace(y, fidx_replace, dftfreq, srate,dim)
%F = dft_fftreplace(y, fidx_replace, dftfreq, srate,dim)
%
%F = fft(y,[],dim);
%F(fidx_replace,:...) = dft(y,dftfreq,srate,dim);
% unless fidx_relpace is empty

if(~exist('dim','var') || isempty(dim))
    dim = 1;
end

if(~exist('srate','var'))
    srate = [];
end

F = fft(y,[],dim);
if(~isempty(fidx_replace))
    Fc = dft(y,dftfreq,srate,dim);
    str = repmat(',:',1,ndims(F)-1);
    eval(['F(fidx_replace' str ') = Fc;']);
end
