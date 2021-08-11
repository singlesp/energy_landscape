function idx = FindListInList(list_to_search,list_to_find)

if(size(list_to_search,1) < size(list_to_search,2)) % n x 1
    list_to_search = list_to_search';
end
n = size(list_to_search,1);

if(size(list_to_find,1) > size(list_to_find,2)) % 1 x m
    list_to_find = list_to_find';
end
m = size(list_to_find,2);

idx = any(repmat(list_to_search,1,m) == repmat(list_to_find,n,1),2);