n = input("Matrix dimension n = ");
A = input("elements of Matrix A: ");
method = input("method choice, D for Doolittle, C for Crout: ", 's');

% n = 4;
% A = [1,1,0,3;2,1,-1,1;3,-1,-1,2;-1,2,3,-1];

LU(A, n, method)

function LU(a,n,method)
    if method == "D"
        m = zeros(n,n);             % initial the n*n matrix m zeros 
        for i = 1:n;     m(i,i) = 1;  end   % let the diagonal elements be 1 
        for j = 1 : n-1
            if abs(a(j,j))<eps
                error('factorization impossible');    % when the zero pivot happens,end the process
            end
            for i = j+1 : n
                mult = a(i,j)/a(j,j);
                m(i,j) = mult;  
                for k = j:n
                    a(i,k) = a(i,k) - mult*a(j,k);
                end
            end
        end
        disp('L=');  disp(m);
        disp('U=');  disp(a);
        % disp("L*U"); disp(m*a);
    elseif method == "C"
        L=zeros(n,n);
        U=zeros(n,n);
        for i=1:n
            U(i,i)=1;
        end
        for k=1:n
            for i=k:n 
                temp_sum=0; 
                for t=1:k-1
                    temp_sum=temp_sum+L(i,t)*U(t,k);
                end
                L(i,k)=a(i,k)-temp_sum; 
            end
            for j=k+1:n 
                temp_sum=0;
                for t=1:k-1
                    temp_sum=temp_sum+L(k,t)*U(t,j); 
                end
                U(k,j)=(a(k,j)-temp_sum)/L(k,k);
            end
        end
        disp('L='); disp(L);
        disp('U='); disp(U);
        % disp("L*U"); disp(L*U);
    end
end