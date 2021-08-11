function d = rivalrydir

if(isunix)
    d = [researchdir '/rivalry/data'];
else
    d = [researchdir '/rivalry/'];
end