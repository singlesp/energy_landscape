function result = mricron(structural_image, functional_image, clut)

if(nargin < 1)
    structural_image = '';
end

if(nargin < 2)
    functional_image = '';
end

if(nargin < 3)
    if(isempty(functional_image))
        clut = '';
    else
        clut = 'actc';
    end
end

structural_image = regexprep(structural_image,'\.img$','.hdr');
functional_image = regexprep(functional_image,'\.img$','.hdr');

if(~isempty(functional_image))
    functional_image = [' -o ' functional_image];
end

if(~isempty(clut))
    clut = [' -c ' clut];
end

cmdstr = sprintf('mricron %s%s%s &',structural_image, functional_image, clut);
result = system(cmdstr);