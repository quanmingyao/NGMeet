function  [Init_Index]  =  Block_matching3D(V_pat, par, Neighbor_arr,Num_arr,SelfIndex_arr)
% Block matching for full band patches;
L            =   length(Num_arr);
Init_Index   =  zeros(par.patnum,L);
sizeV        =  size(V_pat);
for  i  =  1 : L
    Patch           = V_pat(:,:,SelfIndex_arr(i));
    Patch           = Patch(:);
    Neighbors       = Unfold(V_pat(:,:,Neighbor_arr(1:Num_arr(i),i)),[sizeV(1:2),Num_arr(i)],3)';    
    Dist            = sum((repmat(Patch,1,Num_arr(i))-Neighbors).^2);   
    [~, index]      = sort(Dist);
    Init_Index(:,i) =Neighbor_arr(index(1:par.patnum),i);
end
