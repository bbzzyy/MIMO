function [link,Pr] = link3GPP( scenario,d_2d,h_UT)
%Link transition from 3GPP table 7.4.2-1
% scenario:'RMa','UMi','UMa'
% h_UT: the height of user terminal in meters
% d_2d: the 2d distance in meters

switch (scenario)
    case 'RMa'
        if d_2d <= 10
            Pr = 1;
        else
            Pr = exp(-(d_2d-10)/1000);
        end
    case 'UMi'
        if d_2d <= 18
            Pr = 1;
        else
            Pr = 18/d_2d+exp(-d_2d/36)*(1-18/d_2d);
        end
    case 'UMa'
        if d_2d <= 18
            Pr = 1;
        else
            if h_UT <= 13
                C_p = 0;
            else
                C_p = ((h_UT-13)/10)^1.5;
            end
            Pr = (18/d_2d+exp(-d_2d/36)*(1-18/d_2d))*(1+C_p*(5/4)*(d_2d/100)^3*exp(-d_2d/150));
        end
end

x = rand;
if x <= Pr
    link = 'LOS';
else
    link = 'NLOS';
end

end

