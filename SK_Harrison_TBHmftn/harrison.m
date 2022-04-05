%%% author : Yueh-Ting Yao, 2021/04/15

%%% Harrison's distance dependent tight-binding model
%%% sp,sp interaction : https://doi.org/10.1103/PhysRevB.20.2420
%%% sp,d interaction : https://doi.org/10.1103/PhysRevB.21.3214
%%% d,d interaction : Walter A. Harrison, Surface Science 299/300 (1994) 29X-3 10

%%
function HAR = harrison(orb1,orb2,t_sk,dij,por1,por2,rd1,rd2,ro1,ro2,ro3,ro4)

% orb1 : orbital type of first orbital
% orb2 : orbital type of second orbital
% t_sk : t from slaterkoster function
% dij : distance between 2 sites
% por1 : order of distance of sp-d interaction
% por2 : order of distance of d-d interaction
% rd1 : radius of orb1 on first site
% rd2 : radius of orb2 on second site
% ro1 : order of orb1 radius when sp-d interaction
% ro2 : order of orb2 radius when sp-d interaction
% ro3 : order of orb1 radius when d-d interaction
% ro4 : order of orb2 radius when d-d interaction

hbarm = 2*13.6058;      % --- binding energy of hydrogen

if (0<orb1) && (orb1<5)
    if (0<orb2) && (orb2<5)
        HAR  = t_sk*hbarm/dij^2;                                            % sp,sp interaction , formula(1)
    end
	if (4<orb2) && (orb2<10)
        hm   = hbarm*(rd1^ro1)*(rd2^ro2);
        HAR  = t_sk*hm/dij^por1;                                            % sp,d interaction , formula(18)
    end
end

if (4<orb1) && (orb1<10) 
	if (0<orb2) && (orb2<5)
        hm   = hbarm*(rd1^ro1)*(rd2^ro2);
        HAR  = t_sk*hm/dij^por1;                                            % sp,d interaction , formula(18)
    end
	if (4<orb2) && (orb2<10)
        hm   = hbarm*(rd1^ro3)*(rd2^ro4);
        HAR  = t_sk*hm/dij^por2;                                            % d,d interaction , formula(2)
	end
end
