function D = jsdiv(P,Q)
Q = Q./sum(Q);
P = P./sum(P);
M = (P+Q)/2;

D = .5*kldiv(P,M)+.5*kldiv(Q,M);