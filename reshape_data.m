
% u
for j = 1:1:14
    for i = 1:1:15
        u_mat(i,j) = u(i+(j-1)*15);
    end
end

for j = 1:1:15
    for i = 1:1:14
        v_mat(i,j) = v(i+(j-1)*14);
    end
end