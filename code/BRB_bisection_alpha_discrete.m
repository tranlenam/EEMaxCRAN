function [ myalpha ] = BRB_bisection_alpha_discrete( varupbound,varlowbound,e_i,cbv, nUsers, nBS)
%

tempx = varupbound - e_i;
fconnect = zeros(nUsers,1);
for iUser=1:nUsers
    fconnect(iUser)= -sum(tempx(nBS+iUser:nUsers:nBS+nUsers*nBS))+1;
end
fsx = zeros(nBS,1);
for iBS = 1:nBS
    fsx(iBS) = (nnz(varlowbound(nBS+(iBS-1)*nUsers+1:nBS+iBS*nUsers) <= tempx(iBS)) == nUsers);
end
fBS = zeros(nBS,1);
for iBS=1:nBS
    fBS(iBS) = varlowbound(iBS) - sum(tempx(nBS+(iBS-1)*nUsers+1:nBS+iBS*nUsers));
end

if (nnz((fconnect)<=0)== nUsers) && (nnz(fBS <= 0)==nBS) && (nnz(fsx) == nBS)
    myalpha = 1;
else
    myalpha = 0;
end

end

