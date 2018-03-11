close all

xlen = max(x);
ylen = max(y);

Pmat = zeros(ylen,xlen);
utemp = zeros(ylen,xlen);
vtemp = zeros(ylen,xlen);
ttemp = zeros(ylen,xlen);

for i = 1:1:length(p)
    Pmat(y(i),x(i)) = p(i);
    utemp(y(i),x(i)) = u(i);
    vtemp(y(i),x(i)) = v(i);
    ttemp(y(i),x(i)) = T(i);
end

P_fin = zeros(ylen-1,xlen-1);
U_fin = zeros(ylen-1,xlen);
V_fin = zeros(ylen,xlen);
T_fin = zeros(ylen-1,xlen-1);

for i = 1:1:xlen
    for j = 1:1:ylen-1
        U_fin(j,i) = utemp(j,i);
    end
end

for i = 1:1:xlen-1
    for j = 1:1:ylen-1
        P_fin(j,i) = Pmat(j,i);
    end
end

for i = 1:1:xlen-1
    for j = 1:1:ylen-1
        T_fin(j,i) = ttemp(j,i);
    end
end

% Convert to P_fin to atm
P_fin = P_fin/101325;

% Non-dimensionalize U-Velocity
U_fin = U_fin;

figure('Name','u-velocity')
contour(U_fin,10)
colorbar
grid minor
title 'U-Velocity Field (m/s)'
xlabel 'i node'
ylabel 'j node'

figure('Name','v-velocity')
contour(vtemp,10)
colorbar
grid minor
title 'V-Velocity Field (m/s)'
xlabel 'i node'
ylabel 'j node'

figure('Name','Pressure')
contourf(P_fin,10)
colorbar
grid minor
title 'Pressure Field (atm)'
xlabel 'i node'
ylabel 'j node'

figure('Name','Temperature')
contourf(T_fin,10)
colorbar
xlabel 'i node'
ylabel 'j node'

figure('Name','Velocity Plot')
quiver(x,y,u,v)
xlabel 'i node'
ylabel 'j node'
