function Ts = diagmat (Nmat,T)

N3 = 3*Nmat;
oN = ones(Nmat,1);

i3x = [1:3:N3-2].';
i3y = [2:3:N3-1].';
i3z = [3:3:N3].';

ii = [i3x;i3x;i3x;i3y;i3y;i3y;i3z;i3z;i3z];
jj = [i3x;i3y;i3z;i3x;i3y;i3z;i3x;i3y;i3z];
ss = [T(1,1)*oN;T(1,2)*oN;T(1,3)*oN; ...
      T(2,1)*oN;T(2,2)*oN;T(2,3)*oN; ...
      T(3,1)*oN;T(3,2)*oN;T(3,3)*oN];
Ts = sparse(ii,jj,ss,N3,N3);