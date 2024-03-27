function [ Q,R ] = GramSchmidt(A )

Q=zeros(size(A));
R=zeros(size(A,2),size(A,2));
for j=1:size(A,2)
    v=A(:,j);
    for i=1:(j-1)
        R(i,j)=Q(:,i)'*A(:,j);
        v=v-R(i,j)*Q(:,i);
    end
    R(j,j)=norm(v);
    Q(:,j)=v/R(j,j);
end

end

