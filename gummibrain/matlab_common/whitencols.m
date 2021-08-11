function n = whitencols(v)

n = (v-repmat(mean(v,1),size(v,1),1))./repmat(std(v,1),size(v,1),1);