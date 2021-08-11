function d = learningdir

if(isunix)
    d = '~/v1learning/data';
elseif(ispc)
    d = 'c:/users/keith/research/v1learning2';
end