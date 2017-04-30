function [rval, r_x, r_xx] = fuselage_geom(xx, Rmax, length)

xval = xx./length;

rval = Rmax.*((4*xval.*(1-xval)).^0.75);
r_x = atand(3*Rmax.*(((4.*xval.*(1-xval)).^(-0.25)).*(1-2.*xval))/length);
r_xx = -3*Rmax.*((4.*xval.*(1-xval)).^(-5/4) .* (1 - 2.*xval)^2 + 2*(4.*xval.*(1-xval)).^(-0.25))/(length^2);

end