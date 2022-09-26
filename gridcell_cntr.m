function [var_avgcntr,uncrtn,delta_h_hat_nought_squared,a,b,c,Lxx,Lxy,Lyy,Rxh,Ryh] = ...
    gridcell_cntr(var,lon,lat,h_uncrtn,var_avg9,lon_avg9,lat_avg9,gridcntr_lon,gridcntr_lat)
% var is vector of ATL12 geophysical variable to be interpolated to gridcntr_lon, gridcntr_lat
% lon is ATL12 longitude vector for var, in the nine cells
% lat is ATL12 latitude vector for var, in the nine cells
% h_uncrtn is the uncertainty of the heights vector. Input ones([length(df),1]) if want simple average
% var_avg9 is simple average of var or degrees-of-freedom-uncertainty-weighted average of ATL12 var
%    ie. = dot product of (variable,(1./h_uncrtn).^2) / sum((1./h_uncrtn).^2);  eqn(B15)
% lon_avg9 is simple average of var or degrees-of-freedom-uncertainty-weighted average of ATL12 longitude
% lat_avg9 is simple average of var or degrees-of-freedom-uncertainty-weighted average of ATL12 latitude
% gridcntr_lon, gridcntr_lat is center of gridcell


% Methodology follows ATL19 ATBD Appendix B

            hp = (var - var_avg9);  % eqn(B16)
            xp = (lon - lon_avg9);  
            yp = (lat - lat_avg9);  
            
            wi = (1./h_uncrtn).^2;
            Lxx = sum(wi.*xp.*xp);   % eqn(B17)
            Lyy = sum(wi.*yp.*yp);
            Lxy = sum(wi.*xp.*yp);
            Rxh = sum(wi.*xp.*hp);
            Ryh = sum(wi.*yp.*hp);
            
            a = (Rxh*Lyy - Ryh*Lxy) / (Lxx*Lyy - Lxy*Lxy);  % eqn(B18)
            b = (Ryh*Lxx - Rxh*Lxy) / (Lxx*Lyy - Lxy*Lxy);
            c = var_avg9 - (a*lon_avg9 + b*lat_avg9);
            
            var_avgcntr = a*gridcntr_lon + b*gridcntr_lat + c;  
            
            ab = a*xp + b*yp;  % eqn(B19)
            delta_h_hat_nought_squared = length(var)/(length(var)-2) * 1/sum(wi) * sum(wi.*(ab - hp).^2); % eqn(B8)

            La = Lxx*Lyy - Lxy^2;   % denominator eqn(B20)
            xparen = Lyy^2*Lxx - Lyy*Lxy^2;  % parts of eqn(B21)
            yparen = Lyy*Lxx^2 - Lxx*Lxy^2;
            xhatp = gridcntr_lon - lon_avg9;
            yhatp = gridcntr_lat - lat_avg9;

            delta_h_hat_squared = delta_h_hat_nought_squared * ( (1/La)^2 * (xparen*xhatp^2 + yparen*yhatp^2) + 1);  %eqn(B21)

            uncrtn = sqrt(delta_h_hat_squared);

end

