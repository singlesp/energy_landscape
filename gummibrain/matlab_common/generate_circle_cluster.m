function circ = circle_cluster(numpoints,r,c,sigma)

r_noise = r+randn(numpoints,1)*sigma;
theta = rand(numpoints,1)*2*pi;
circ = [cos(theta) sin(theta)].*repmat(r_noise,1,2) + repmat(c,numpoints,1);
