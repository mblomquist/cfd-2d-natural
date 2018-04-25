
for j = 1:1:14
    for i = 1:1:14
        P(i,j) = fielddatapressure(i+(j-1)*14);
        T(i,j) = temperature(i+(j-1)*14);
    end
end
