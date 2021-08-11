function W = cosinewindow(N)

W = cos(pi*[0:N-1]/(N-1) - pi/2);
W = W/sum(W);