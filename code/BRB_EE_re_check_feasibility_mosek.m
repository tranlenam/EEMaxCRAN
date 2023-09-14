function [ stt , W, t_b,beamformer,time_solving] = BRB_EE_re_check_feasibility_mosek( channel,nUsers,nBS,nTx,sigma, z_lb,z_ub,x_i,s_ub,t_lb,t_ub ,P_ant,  power, aeff)

x_ub =  x_i;


prob = [];
% varlength = 2*nBS*nUsers*nTx+nBS*nUsers+(nUsers*nBS-1)*nUsers*nBS*2+nBS+nUsers*nBS+(nUsers*nBS-1)*nUsers*nBS*2+2*nUsers*nBS+nUsers*nBS+nBS+nBS+1;
varlength = 3*2*nBS*nUsers*nTx+nUsers+(nUsers-1)*nUsers*2+nUsers+2*nBS*nTx+nBS*nUsers+nBS+nBS*nTx;
% Im ==0
a1=zeros(nUsers*nBS,varlength);
for iUser=1:nUsers
    a1(iUser,(iUser-1)*2*nTx*nBS+1:(iUser)*2*nTx*nBS ) =[ myvec([imag(channel(iUser,:)); real(channel(iUser,:))])'];
end

nlength1=2*nBS*nUsers*nTx;
% Real() >= a2
a2=zeros(nUsers,varlength);

for iUser=1:nUsers
    a2(iUser,(iUser-1)*2*nTx*nBS+1:(iUser)*2*nTx*nBS ) =[ myvec([ real(channel(iUser,:)); -imag(channel(iUser,:))])'];
    a2(iUser,2*nTx*nBS*nUsers+iUser)=-sqrt(exp(z_lb(iUser))-1); %
end

%interference
nlength2=nlength1+nUsers;
a3=[];
for iUser = 1:nUsers
    for jUser = 1:nUsers
        if jUser~=iUser
            temp = zeros(2,varlength);
            temp(1,(jUser-1)*2*nTx*nBS+1:(jUser)*2*nTx*nBS )...
                =[ myvec([ real(channel(iUser,:)); - imag(channel(iUser,:))])'];
            %                     temp(1,2*nTx*nBS*nUsers+nBS*nUsers+(iBS-1)*(nUsers*nBS-1)+(jj-1)*nUsers+jUser)=-1;
            temp(2,(jUser-1)*2*nTx*nBS+1:(jUser)*2*nTx*nBS )...
                =[ myvec([ imag(channel(iUser,:));  real(channel(iUser,:))])'];
            a3= [a3;temp];
        end
    end
end


for i=1:size(a3,1)
    a3(i,nlength2+i)=-1;
end

nlength3=nlength2+2*(nUsers-1)*nUsers + nUsers ;% nUser for sigma

a4 = []; %sum(u_bk) <=s*P_b
for iBS=1:nBS
    
    temp = zeros(1,varlength);
    temp(1,nlength3+(iBS-1)*nTx+1:nlength3+(iBS)*nTx)= 1;
    a4 = [a4; temp];
end

nlength4=nlength3 +  nBS*nTx ; %nBS*nUsers*nTx u_bki
% nlength5= nlength4+ nBS*nUsers; %x_bk
% a5=[];
% for iBS=1:nBS
%     for iUser=1:nUsers
%         temp = zeros(1,varlength);
%         temp(1,nlength3+(iBS-1)*nUsers*nTx+(iUser-1)*nTx+1:nlength3+(iBS-1)*nUsers*nTx+iUser*nTx)=1;
%         boundlb_a5((iBS-1)*nUsers+iUser) = x_ub((iBS-1)*nUsers+iUser)*1e-7;
%         boundub_a5((iBS-1)*nUsers+iUser) = x_ub((iBS-1)*nUsers+iUser)*t_ub(iBS);
%         a5=[a5;temp];
%     end
% end
nlength6=nlength4;

% a7 add nUsers z



for iUser=1:nUsers
    for iBS=1:nBS
        if x_ub((iBS-1)*nUsers+iUser) == 0
            bm_ub(2*(iUser-1)*nBS*nTx+2*(iBS-1)*nTx+1:2*(iUser-1)*nBS*nTx+2*(iBS)*nTx) = 0;
            bm_lb(2*(iUser-1)*nBS*nTx+2*(iBS-1)*nTx+1:2*(iUser-1)*nBS*nTx+2*(iBS)*nTx) = 0;
            %             u_ub((iBS-1)*nUsers+iUser) = 0;
            %             u_lb((iBS-1)*nUsers+iUser) = 0;
        else bm_ub(2*(iUser-1)*nBS*nTx+2*(iBS-1)*nTx+1:2*(iUser-1)*nBS*nTx+2*(iBS)*nTx) = inf;
            bm_lb(2*(iUser-1)*nBS*nTx+2*(iBS-1)*nTx+1:2*(iUser-1)*nBS*nTx+2*(iBS)*nTx) = -inf;
            %             u_ub((iBS-1)*nUsers+iUser) = s_ub(iBS)*t_ub(iBS)/2;
            %             u_lb((iBS-1)*nUsers+iUser) = s_lb(iBS)*t_lb(iBS)/2;
        end
        
    end
end
% u_ub = kron(t_ub,ones(1,nUsers));
% u_lb = kron(t_lb,ones(1,nUsers));
% prob.c = a8;

% prob.cones = cell(nUsers*nBS+nUsers,1);
%% real() >= sqrt(interference)
for i=1:nUsers
    prob.cones{i}.type = 'MSK_CT_QUAD';
    prob.cones{i}.sub = [nlength1+i,nlength2+2*(i-1)*(nUsers-1)+1:nlength2+2*i*(nUsers-1),nlength3-nUsers+i];
end
cone = nUsers;
%% norm(w_bk)^2< x_ub*power => norm(w_bk)< sqrt(x_ub*power)
sqr_powxbk = sqrt(x_ub*power) ;
nlength7 = nlength6+nUsers*nBS;
for i=1:nUsers
    for k=1:nBS
        prob.cones{cone+(i-1)*nBS+k}.type = 'MSK_CT_QUAD';
        prob.cones{cone+(i-1)*nBS+k}.sub = [nlength6+(k-1)*nUsers+i,(i-1)*nBS*2*nTx+(k-1)*2*nTx+1:(i-1)*nBS*2*nTx+(k)*2*nTx];
    end
end

% sqr_powxbk = power/2*(x_ub*power) ;
% nlength7 = nlength6+nUsers*nBS;
% for i=1:nUsers
%     for k=1:nBS
%     prob.cones{cone+(i-1)*nBS+k}.type = 'MSK_CT_RQUAD';
%     prob.cones{cone+(i-1)*nBS+k}.sub = [nlength6+(k-1)*nUsers+i,nlength4+(k-1)*nUsers+i,(i-1)*nBS*2*nTx+(k-1)*2*nTx+1:(i-1)*nBS*2*nTx+(k)*2*nTx];
%     end
% end
cone = cone + nUsers*nBS;


%% sum(||w_bk||^2) <= s_ub*P_b => sqrt(sum(||w_bk||^2)) <= sqrt(s_ub*P_b)
a9=zeros(2*nUsers*nBS*nTx,varlength);
a9(1:2*nUsers*nBS*nTx,1:2*nUsers*nBS*nTx) = diag(ones(1,2*nUsers*nBS*nTx));
a9(1:2*nUsers*nBS*nTx,nlength7+1:nlength7+2*nUsers*nBS*nTx) =  -diag(ones(1,2*nUsers*nBS*nTx));
nlength8 = nlength7+2*nUsers*nBS*nTx;
beam_B_index = [];
for k=1:nBS
    temp =[];
    for i=1:nUsers
        temp = [temp nlength7+(i-1)*nBS*2*nTx+(k-1)*2*nTx+1:nlength7+(i-1)*nBS*2*nTx+(k)*2*nTx];
    end
    beam_B_index(k,:)= temp;
end

sqr_powsb = sqrt(s_ub*power) ;
nlength9 = nlength8 + nBS;
for i=1:nBS
    prob.cones{cone+i}.type = 'MSK_CT_QUAD';
    prob.cones{cone+i}.sub = [nlength8+i,beam_B_index(i,:)];
end
cone= cone + nBS;
%% |w_bki| <= u_bki
a10=zeros(2*nUsers*nBS*nTx,varlength);
a10(1:2*nUsers*nBS*nTx,1:2*nUsers*nBS*nTx) = diag(ones(1,2*nUsers*nBS*nTx));
a10(1:2*nUsers*nBS*nTx,nlength9+1:nlength9+2*nUsers*nBS*nTx) =  -diag(ones(1,2*nUsers*nBS*nTx));
nlength10 = nlength9+2*nUsers*nBS*nTx;
% for i=1:nUsers
for j=1:nBS
    for k = 1:nTx
        prob.cones{cone+ (j-1)*nTx + k}.type = 'MSK_CT_QUAD';
        prob.cones{cone+ (j-1)*nTx + k}.sub = [nlength3 +  (j-1)*nTx + k, nlength9 + (j-1)*2*nTx+(k-1)*2+1:nBS*nTx*2:nlength9+2*nTx*nBS*nUsers,nlength9 + (j-1)*2*nTx+(k-1)*2+2:nBS*nTx*2:nlength9+2*nTx*nBS*nUsers];
    end
end
% end
cone = cone+nBS*nTx;
%% sum(u_bki^2) <= s_b*P_ant => sqrt(sum(u_bki^2)) <= sqrt(s_b*P_ant)
a11=zeros(nUsers*nBS*nTx,varlength);
a11(1:nBS*nTx,nlength3+1:nlength3+nBS*nTx) = diag(ones(1,nBS*nTx));
a11(1:nBS*nTx,nlength10+1:nlength10+nBS*nTx) =  -diag(ones(1,nBS*nTx));
nlength11 = nlength10+nBS*nTx;
sqr_powsba = sqrt(kron(s_ub,ones(1,nTx)).*(P_ant*ones(1,nBS*nTx))) ;
nlength12 = nlength11+ nBS*nTx;
for j=1:nBS
    for k=1:nTx
        %         for i=1:nUsers
        prob.cones{cone+(j-1)*nTx+k}.type = 'MSK_CT_QUAD';
        prob.cones{cone+(j-1)*nTx+k}.sub = [nlength11+(j-1)*nTx+k,nlength10 + (j-1)*nTx+k ];
        %         end
    end
end
cone = cone + nTx*nBS;
%% sum(u_bki^2) <=obj
% a12=zeros(nUsers*nBS*nTx,varlength);
% a12(1:nUsers*nBS*nTx,nlength3+1:nlength3+nUsers*nBS*nTx) = diag(ones(1,nUsers*nBS*nTx));
% a12(1:nUsers*nBS*nTx,nlength12+1:nlength12+nUsers*nBS*nTx) =  -diag(ones(1,nUsers*nBS*nTx));
% nlength13 = nlength12+nUsers*nBS*nTx;
%     prob.cones{cone+1}.type = 'MSK_CT_RQUAD';
%     prob.cones{cone+1}.sub = [nlength13+1,nlength13+2,nlength12+1:nlength12+nUsers*nBS*nTx];
%%

prob.c = [zeros(1,nlength3) ones(1,nBS*nTx) zeros(1,nUsers*nBS+nBS+nBS*nTx+ 2*2*nTx*nUsers*nBS + nTx*nBS) ];
% prob.c = [zeros(1,varlength)];
prob.a = sparse([a1;a2;a3;a4;a9;a10;a11]);

prob.blc = [zeros(1,size(a1,1)) zeros(1,size(a2,1)) 0*ones(1,size(a3,1)) t_lb.*s_ub    zeros(1,size(a9,1))  zeros(1,size(a10,1))  zeros(1,size(a11,1)) ];
prob.buc = [zeros(1,size(a1,1)) 0*ones(1,size(a2,1)) zeros(1,size(a3,1)) t_ub.*s_ub    zeros(1,size(a9,1))  zeros(1,size(a10,1))  zeros(1,size(a11,1)) ];
prob.blx = [bm_lb 0*ones(1,nUsers), -inf*ones(1,(nUsers-1)*nUsers*2) sigma*ones(1,nUsers) zeros(1,nBS*nTx)         sqr_powxbk bm_lb  sqr_powsb  bm_lb zeros(1,nBS*nTx)      sqr_powsba       ];
prob.bux = [bm_ub inf*ones(1,nUsers) inf*ones(1,(nUsers-1)*nUsers*2) sigma*ones(1,nUsers) sqrt(P_ant)*ones(1,nBS*nTx)    sqr_powxbk bm_ub  sqr_powsb  bm_ub sqrt(P_ant)*ones(1,nBS*nTx) sqr_powsba       ];


% prob.a = sparse([a1;a2;a3;a4;a9;a10;a11;a12]);
% prob.blc = [zeros(1,size(a1,1)) zeros(1,size(a2,1)) 0*ones(1,size(a3,1)) t_lb.*ones(1,size(a4,1))   sum(z_lb)  zeros(1,size(a9,1))  zeros(1,size(a10,1))  zeros(1,size(a11,1))  zeros(1,size(a12,1))];
% prob.buc = [zeros(1,size(a1,1)) 0*ones(1,size(a2,1)) zeros(1,size(a3,1)) t_ub.*ones(1,size(a4,1))   sum(z_ub)  zeros(1,size(a9,1))  zeros(1,size(a10,1))  zeros(1,size(a11,1))  zeros(1,size(a12,1))];
% prob.blx = [bm_lb 0*ones(1,nUsers), -inf*ones(1,(nUsers-1)*nUsers*2) sigma*ones(1,nUsers) zeros(1,nUsers*nBS*nTx)      x_lb z_lb   sqr_powxbk   sqr_powsb  sqr_powsba   bm_lb bm_lb  zeros(1,nUsers*nBS*nTx) zeros(1,nUsers*nBS*nTx)             0         1/2];
% prob.bux = [bm_ub inf*ones(1,nUsers) inf*ones(1,(nUsers-1)*nUsers*2) sigma*ones(1,nUsers) P_ant*ones(1,nUsers*nBS*nTx) x_ub z_ub   sqr_powxbk   sqr_powsb  sqr_powsba   bm_ub bm_ub  P_ant*ones(1,nUsers*nBS*nTx) P_ant*ones(1,nUsers*nBS*nTx) inf 1/2];

%%

[r,res]       = mosekopt('minimize echo(0) info',prob);
time_solving = res.info.MSK_DINF_OPTIMIZER_TIME;
if ~r
    W = res.sol.itr.xx;
    % res.sol.itr.solsta
    % if strncmp(res.rcodestr, 'MSK_RES_TRM_STALL',17)
    %     aadas
    % end
    % res.sol.itr.solsta
    
    % res.sol.itr.pobjval
    % clear j
    if strcmp(res.sol.itr.solsta, 'OPTIMAL') ||  strcmp(res.sol.itr.solsta, 'NEAR_OPTIMAL')
        stt = 0;
        beamformer=zeros(nBS*nTx,nUsers);
        for iUser=1:nUsers
            %
            reBeamformer = W((iUser-1)*2*nTx*nBS+[1:2:2*nTx*nBS-1],1);
            imBeamformer = W((iUser-1)*2*nTx*nBS+[2:2:2*nTx*nBS],1);
            beamformer(:,iUser) = reBeamformer + 1i*imBeamformer;

        end
        u_bki = (W(nlength3+1:nlength4))';
        u_bki = reshape(u_bki,nTx,nBS);
        t_b = sum(u_bki);
        %     t_bl = t_blxmin;
    else
        stt = 1;
        beamformer = -inf;
        %     imBeamformer = -inf;
        t_b = 0;
    end
else
    stt = 1;
    beamformer = -inf;
    %     imBeamformer = -inf;
    t_b = 0;
    W = 0;
end