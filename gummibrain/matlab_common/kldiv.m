function D = kldiv(P,Q)
Q = Q./sum(Q);
P = P./sum(P);

s = P.*log2(P./Q);
s(isnan(s) | isinf(s)) = 0;
D = sum(s);
