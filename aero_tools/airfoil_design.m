function airfoil_design(fname, ftype, dx, c_r, t_c, wt, sim_type)

%% Define domain

% dx = 0.01;
% t_c = 0.7;
% wt = 0.5;

if exist([fname '.' ftype], 'file')
    delete([fname '.' ftype]);
    fprintf('Warning: File Deleted\n');
end

fid = fopen([fname '.' ftype], 'wt');

if strcmp(sim_type, 'openvsp')
    bott = (0:dx:1);
    top = (0:dx:1);
    coords = [top', 2.*((1-wt)*t_c).*((1-top).*top)'; % top of airfoil
           bott', 2.*t_c.*(-(1-bott).*bott)']; % bottom of airfoil
    
    fprintf(fid, 'DEMO GEOM AIRFOIL FILE\n');
    fprintf(fid, [fname ' Supersonic Airfoil\n']);
    fprintf(fid, '0\tSym Flag\n');
    fprintf(fid, '%i\t Num Pnts Upper\n', length(coords(:,1))*0.5);
    fprintf(fid, '%i\t Num Pnts Lower\n', length(coords(:,1))*0.5);
    for i = 1:length(coords)
        fprintf(fid, '%3.7f\t%4.7f\n', coords(i,1), coords(i,2));
        if i == 0.5*length(coords(:,1))
            fprintf(fid, '\n');
        end
    end
    fclose(fid);
    
elseif strcmp(sim_type, 'ansys')
%     coords = [(0:dx:1)', 2.*((1-wt)*t_c).*((1-(0:dx:1)).*(0:dx:1))'; % top of airfoil
%            ((1-dx):-dx:0)', 2.*t_c.*(-(1-((1-dx):-dx:0)).*((1-dx):-dx:0))']; % bottom of airfoil
    top = (0:dx:1);
    bott = ((1-dx):-dx:0);
    coords = [top', 2.*(t_c).*((1-top).*top)'; % top of airfoil
       bott', 2.*t_c.*(-(1-bott).*bott)']; % bottom of airfoil
       
    fprintf(fid, '\t\t#Supersonic Airfoil\n\n');
    fprintf(fid, '#group\t#point\t#x_cord\t#y_cord\t#z_cord\n');
    for i = 1:length(coords)-1
    %     fprintf(fid, '%s Supersonic Airfoil\n\n', WING.supersonic.name);
        fprintf(fid, '%i\t%i\t%12.7f\t%12.7f\t%12.7f\n', 1, i, coords(i,1), coords(i,2), 0.0000);
    end
    fprintf(fid, '%i\t%i\t%12.7f\t%12.7f\t%12.7f\n', 1, 0, coords(1,1), coords(1,2), 0.0000);
    fclose(fid);
    
elseif strcmp(sim_type, 'Selig')
    top = 1:-dx:0;
    bott = dx:dx:1;
    coords = [top', 2.*((1-wt)*t_c).*((1-top).*top)'; % top of airfoil
       bott', 2.*t_c.*(-(1-bott).*bott)']; % bottom of airfoil
       
    fprintf(fid, [fname ' Supersonic Airfoil\n']);
    for i = 1:length(coords)
        fprintf(fid, '%3.7f\t%4.7f\n', coords(i,1), coords(i,2));
        if i == 0.5*length(coords(:,1))
            fprintf(fid, '\n');
        end
    end
    fclose(fid);
    
elseif strcmp(sim_type, 'solidworks')
    top = 1:-dx:0;
    bott = dx:dx:1;
    coords = [top', 2.*((1-wt)*t_c).*((1-top).*top)'; % top of airfoil
       bott', 2.*t_c.*(-(1-bott).*bott)']; % bottom of airfoil
       
    fprintf(fid, [fname ' Supersonic Airfoil\n']);
    for i = 1:length(coords)
        fprintf(fid, '%3.7f\t%4.7f\t%i\n', c_r*coords(i,1), c_r*coords(i,2), 0.0);
        if i == 0.5*length(coords(:,1))
            fprintf(fid, '\n');
        end
    end
    
end

figure();
plot(coords(:,1), coords(:,2));
title('Biconvex Airfoil');

end