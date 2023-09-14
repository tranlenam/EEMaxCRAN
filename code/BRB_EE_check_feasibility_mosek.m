function [ sta , cobj,cbup, x_bk,x_bkl, t_bl,z_k,z_u,s_b, fsb,time_solving] = BRB_EE_check_feasibility_mosek( channel,nUsers,nBS,nTx,sigma,scalefactor, lowerboundx,upperboundx,C_b,P_ac,P_sl,P_const,P_SBa, P_ant,power, aeff,m,cbv,cup)

t_lb =  lowerboundx(nBS+nUsers*nBS+nUsers+1:nBS+nUsers*nBS+nUsers+nBS)*scalefactor^2;
t_ub =  upperboundx(nBS+nUsers*nBS+nUsers+1:nBS+nUsers*nBS+nUsers+nBS)*scalefactor^2;
z_lb =  lowerboundx(nBS+nUsers*nBS+1:nBS+nUsers*nBS+nUsers);
z_ub =  upperboundx(nBS+nUsers*nBS+1:nBS+nUsers*nBS+nUsers);
x_ub =  upperboundx(nBS+1:nBS+nUsers*nBS);
x_lb =  lowerboundx(nBS+1:nBS+nUsers*nBS);
s_ub =  upperboundx(1:nBS);
s_lb =  lowerboundx(1:nBS);



%% compute upper and lower bound of constraints w.r.t input variable set.
fbhu = zeros(nBS,1);
fbhl = zeros(nBS,1);
for iBS=1:nBS % compute upper and lower bound of FH capacity w.r.t input variable set.
    fbhu(iBS)= x_ub((iBS-1)*nUsers+1:iBS*nUsers)*z_ub';
    fbhl(iBS)= x_lb((iBS-1)*nUsers+1:iBS*nUsers)*z_lb';
end


cbup = sum(z_ub)/(aeff*sum(t_lb)/scalefactor^2+ (P_ac-P_sl)*max([sum(s_lb(1:nBS))])+P_SBa*max([sum(fbhl.^m) ])+P_const);
cobj = sum(z_lb)/(aeff*sum(t_ub)/scalefactor^2+ (P_ac-P_sl)*max([sum(s_ub(1:nBS))])+P_SBa*min([C_b]) +P_const);
fsb = 1;
sta = 1;
z_k = z_lb;
z_u = z_ub;
x_bk = x_ub;
x_bkl=x_lb;
t_bl = t_lb/scalefactor^2;
s_b = s_lb;
beamformer = inf;
time_solving = 0;
if sum(z_lb) > C_b*sum(s_ub)
    return
end
if nnz((fbhl>C_b))
    return
else
    %%
    if nnz(x_ub - x_lb) % if binary variable not converged
        prob = [];
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
        
        %% real() >= sqrt(interference)
        for iUser=1:nUsers
            prob.cones{iUser}.type = 'MSK_CT_QUAD';
            prob.cones{iUser}.sub = [nlength1+iUser,nlength2+2*(iUser-1)*(nUsers-1)+1:nlength2+2*iUser*(nUsers-1),nlength3-nUsers+iUser];
        end
        cone = nUsers;
        %% norm(w_bk)^2< x_ub*power => norm(w_bk)< sqrt(x_ub*power)
        sqr_powxbk = sqrt(x_ub*power) ;
        nlength7 = nlength6+nUsers*nBS;
        for iUser=1:nUsers
            for iBS=1:nBS
                prob.cones{cone+(iUser-1)*nBS+iBS}.type = 'MSK_CT_QUAD';
                prob.cones{cone+(iUser-1)*nBS+iBS}.sub = [nlength6+(iBS-1)*nUsers+iUser,(iUser-1)*nBS*2*nTx+(iBS-1)*2*nTx+1:(iUser-1)*nBS*2*nTx+(iBS)*2*nTx];
            end
        end
        
        
        cone = cone + nUsers*nBS;
        
        
        %% sum(||w_bk||^2) <= s_ub*P_b => sqrt(sum(||w_bk||^2)) <= sqrt(s_ub*P_b)
        a9=zeros(2*nUsers*nBS*nTx,varlength);
        a9(1:2*nUsers*nBS*nTx,1:2*nUsers*nBS*nTx) = diag(ones(1,2*nUsers*nBS*nTx));
        a9(1:2*nUsers*nBS*nTx,nlength7+1:nlength7+2*nUsers*nBS*nTx) =  -diag(ones(1,2*nUsers*nBS*nTx));
        nlength8 = nlength7+2*nUsers*nBS*nTx;
        beam_B_index = [];
        for k=1:nBS
            temp =[];
            for iUser=1:nUsers
                temp = [temp nlength7+(iUser-1)*nBS*2*nTx+(k-1)*2*nTx+1:nlength7+(iUser-1)*nBS*2*nTx+(k)*2*nTx];
            end
            beam_B_index(k,:)= temp;
        end
        
        sqr_powsb = sqrt(s_ub*power) ;
        nlength9 = nlength8 + nBS;
        for iBS=1:nBS
            prob.cones{cone+iBS}.type = 'MSK_CT_QUAD';
            prob.cones{cone+iBS}.sub = [nlength8+iBS,beam_B_index(iBS,:)];
        end
        cone= cone + nBS;
        %% sum(|w_bki| <= u_bki
        a10=zeros(2*nUsers*nBS*nTx,varlength);
        a10(1:2*nUsers*nBS*nTx,1:2*nUsers*nBS*nTx) = diag(ones(1,2*nUsers*nBS*nTx));
        a10(1:2*nUsers*nBS*nTx,nlength9+1:nlength9+2*nUsers*nBS*nTx) =  -diag(ones(1,2*nUsers*nBS*nTx));
        nlength10 = nlength9+2*nUsers*nBS*nTx;
        % for iUser=1:nUsers
        for jBS=1:nBS
            for k = 1:nTx
                prob.cones{cone+ (jBS-1)*nTx + k}.type = 'MSK_CT_QUAD';
                prob.cones{cone+ (jBS-1)*nTx + k}.sub = [nlength3 +  (jBS-1)*nTx + k, nlength9 + (jBS-1)*2*nTx+(k-1)*2+1:nBS*nTx*2:nlength9+2*nTx*nBS*nUsers,nlength9 + (jBS-1)*2*nTx+(k-1)*2+2:nBS*nTx*2:nlength9+2*nTx*nBS*nUsers];
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
        for jBS=1:nBS
            for iTx=1:nTx
                %         for iUser=1:nUsers
                prob.cones{cone+(jBS-1)*nTx+iTx}.type = 'MSK_CT_QUAD';
                prob.cones{cone+(jBS-1)*nTx+iTx}.sub = [nlength11+(jBS-1)*nTx+iTx,nlength10 + (jBS-1)*nTx+iTx ];
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
        
        prob.c = [zeros(1,nlength3) ones(1,nBS*nTx) zeros(1,nUsers*nBS+nBS+nBS*nTx+ 2*2*nUsers*nTx*nBS + nTx*nBS)  ];
        prob.a = sparse([a1;a2;a3;a4;a9;a10;a11]);
        
        prob.blc = [zeros(1,size(a1,1)) zeros(1,size(a2,1)) 0*ones(1,size(a3,1)) (s_lb.*t_lb).*ones(1,size(a4,1))    zeros(1,size(a9,1))  zeros(1,size(a10,1))  zeros(1,size(a11,1))  ];
        prob.buc = [zeros(1,size(a1,1)) 0*ones(1,size(a2,1)) zeros(1,size(a3,1)) (s_ub.*t_ub).*ones(1,size(a4,1))   zeros(1,size(a9,1))  zeros(1,size(a10,1))  zeros(1,size(a11,1))  ];
        prob.blx = [bm_lb 0*ones(1,nUsers), -inf*ones(1,(nUsers-1)*nUsers*2) sigma*ones(1,nUsers) zeros(1,nBS*nTx)             sqr_powxbk bm_lb  sqr_powsb  bm_lb zeros(1,nBS*nTx)      sqr_powsba      ];
        prob.bux = [bm_ub inf*ones(1,nUsers) inf*ones(1,(nUsers-1)*nUsers*2) sigma*ones(1,nUsers) sqrt(P_ant)*ones(1,nBS*nTx)    sqr_powxbk bm_ub  sqr_powsb  bm_ub sqrt(P_ant)*ones(1,nBS*nTx) sqr_powsba      ];
        
        
        %%
        
        [~,res]       = mosekopt('minimize echo(0) info',prob);
        time_solving =  res.info.MSK_DINF_OPTIMIZER_TIME;
        
        if strcmp(res.rcodestr, 'MSK_RES_OK') || strcmp(res.rcodestr, 'MSK_RES_TRM_STALL')
            W = res.sol.itr.xx;
            if (strcmp(res.sol.itr.solsta, 'OPTIMAL') || strcmp(res.sol.itr.solsta, 'NEAR_OPTIMAL'))
                
                sta = 0;
                obj = res.sol.itr.pobjval ;
                
                x = (W(nlength4+1:nlength4+nUsers*nBS))';
                u_bk = (W(nlength3+1:nlength4))';
                u_bk = reshape(u_bk,nTx,nBS);
                t_blxmin = sum(u_bk);
                t_bl = t_blxmin;
                
                beamformer=zeros(nBS*nTx,nUsers);
                for iUser=1:nUsers
                    %
                    reBeamformer = W((iUser-1)*2*nTx*nBS+[1:2:2*nTx*nBS-1],1);
                    imBeamformer = W((iUser-1)*2*nTx*nBS+[2:2:2*nTx*nBS],1);
                    beamformer(:,iUser) = reBeamformer + 1i*imBeamformer;
                    for iBS=1:nBS
                        beampow(iBS,iUser) = norm(beamformer((iBS-1)*nTx+1:(iBS)*nTx,iUser));
                    end
                    %
                end
                x_b = double(beampow >= 1e-8);
                
                bhcxmin = x_b*z_lb';
                s_b = max(x_b');
                
                if ~nnz(bhcxmin > C_b)
                    fsb = 0;
                else
                end
                %
            end
        end
        % else a=1;
        
    else % if binary variables converged, check if continous variables are feasible
        [ stt , ~, t_bl, ~,time_solving ] = BRB_EE_re_check_feasibility_mosek( channel,nUsers,nBS,nTx,sigma, z_lb,z_ub,x_lb,s_ub,t_lb,t_ub ,P_ant,  power, aeff);
        if stt == 0
            sta = 0;
            fsb = 0;
        end
        
        
    end
    if ~(nnz(t_bl >= t_lb)==nBS)
        t_bl = t_lb;
    end
    
end


end
