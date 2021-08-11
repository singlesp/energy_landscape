
function out = convert_int8(in)
%from zmliu code

in(in>intmax('int8'))=intmax('int8');
in(in<-intmax('int8'))=-intmax('int8');
out=int8(in);