%enforcing rank-2 constraint for fundamental matrix and scale s.t F(3,2) =
%1
function FPrime = enforceFundamentalConstraint(F)
[uf,sf,vf] = svd(F);
FPrime = uf*diag([sf(1) sf(5) 0])*(vf');
FPrime = FPrime/FPrime(3,2);
end