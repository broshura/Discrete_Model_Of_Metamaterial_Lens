clc;
syms x y
disp(culc_vpo(0, 1, 1, 1, 0.2, 0.1));
a = int(culc_vpo(x, 1, 1, 1, 0.2, 0.1), 0, 2*pi);
disp(a)
b = int(y(x), x, 0, 1)

function y = v(x)
    v = x
end

function VPO=culc_vpo(q, x_rel, y_rel, z_rel, r1, r2)

rho_rel=sqrt((x_rel - r2*sin(q)).^2 + y_rel.^2);     %horizontal distance between the loop 1 and element dl of loop 2
kap_sq=4*r1*rho_rel./((r1+rho_rel).^2+(z_rel - r2*cos(q)).^2); %kappa-parameter needed for elliptic functions
[myK,myE]=ellipke(kap_sq);                  %elliptic functions K,E
PROTIRKA = 1e-7;    % very small offset to prevent division by zero for unlucky mutual positions
kap=sqrt(kap_sq) + PROTIRKA; 
Ahelp=1/(2*pi).*sqrt(r1./(rho_rel+PROTIRKA)).*((2./kap-kap).*myK-2./kap.*myE);  %Aphi vector potential, only phi-component with respect to loop1
VPO=r2*Ahelp.*(y_rel*cos(q))./(rho_rel+PROTIRKA);     %this is element (Aphi*r2) which will be integrated over phi=0..2*pi;
end