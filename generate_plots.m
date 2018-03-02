close all

xlen = max(x);
ylen = max(y);

Pmat = zeros(ylen,xlen);
umat = zeros(ylen,xlen);
vmat = zeros(ylen,xlen);
tmat = zeros(ylen,xlen);

for i = 1:1:length(p)
    Pmat(y(i),x(i)) = p(i);
    umat(y(i),x(i)) = u(i);
    vmat(y(i),x(i)) = v(i);
    tmat(y(i),x(i)) = T(i);
end

figure('Name','u-velocity')
contour(umat,20)
colorbar

figure('Name','v-velocity')
contour(vmat,20)
colorbar

figure('Name','Pressure')
contourf(Pmat,20)
colorbar

figure('Name','Temperature')
contourf(tmat,50)
colorbar