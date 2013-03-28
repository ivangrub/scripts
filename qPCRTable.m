[m,~] = size(A);
nump = 2;
Err = zeros(m,1);
Ct = zeros(m,1);
Delta = zeros(m-nump,3);

j = 1;
for i = 1:m
   
    Ct(i) = -(A(i,1)-A(nump*j));
    Err(i) = sqrt(A(i,2)^2+A(nump*j,2)^2);
    if mod(i,nump) == 0
        if (m-i) < nump
            continue
        end
        j = j + 1;
    end
end

k = 1;
for i = 1:m-nump 
    Delta(i,1) = 2^((Ct(i+nump)-Ct(k))-Err(i+nump));
    Delta(i,2) = 2^((Ct(i+nump)-Ct(k))+Err(i+nump));
    Delta(i,3) = mean([Delta(i,1),Delta(i,2)]);
    if mod(i+nump,nump) == 0
        k = 1;
    else
        k = k + 1;
    end
    
end


