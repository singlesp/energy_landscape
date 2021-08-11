%function stackfillplot(x,Y,do_norm)
% x = 1xN x-values
% Y = MxN y-values (each row is a trace)
% do_norm = true: stack sums up to 1. false(default): just sum them up.

function stackfillplot(x,Y,do_norm)

if(nargin < 3)
    do_norm = false;
end

if(size(x,1) > size(x,2))
    x = x'; %make sure the x's are row-vector, not column
end

c = get(gca,'ColorOrder');
sz_c = size(c,1);

[num_traces trace_length] = size(Y);

Y = cumsum(abs(Y),1); %add up magnitudes of previous traces so they stack up

if(do_norm)
    Y=Y./(ones(num_traces,1)*Y(end,:)); %make the last trace=1 if do_norm
end

%add a line of zeros below the first trace
Y = [zeros(1,trace_length); Y];

%concat reversed copy of NEXT trace to end of each trace so that the
%polygon fills in all of the space between trace i and trace i+1
Y = [Y Y([2:end 1],end:-1:1)]; 

%concat reversed copy of x so polygon fills 1->N for values of trace i,
%then N->1 for values of trace i+1
x = [x x(end:-1:1)]; 

hold on;
for i = 1:num_traces
    patch(x,Y(i,:),c(mod(i-1,sz_c)+1,:),'EdgeColor','none');
end
hold off;