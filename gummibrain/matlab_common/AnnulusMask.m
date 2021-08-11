function pict = AnnulusMask(radius, inner_frac)

if(nargin < 2)
    inner_frac = .5;
end

inner_radius = radius*inner_frac;

[x y] = meshgrid(linspace(-radius,radius,radius*2));
r = sqrt(x.^2+y.^2);

pict = r >= inner_radius & r <= radius;