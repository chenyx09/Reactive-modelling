clear
global app
app=1;
params
load LC_data
phi_p=phi_positive;
iter = 1;
sol=[];
phi_negative=phi_negative_anchor;
while iter<20
    if iter==1
        new_phi_negative = lane_change_staliro(0);
    else
        new_phi_negative = lane_change_staliro(1);
    end
    phi_negative= [phi_negative;new_phi_negative];
    disp(['Number of negative samples from falsification: ',num2str(size(phi_negative,1))])
    if isempty(new_phi_negative)
        disp('Success!')
        break
    end
%     CEGAR_SVM_v2_1region;
    SVM_reactive_bound;
    if flag
%         threshold_x = sol.threshold_x;
        
%         disp(['threshold_x = ', num2str(threshold_x)])

%         neg_mis=[];
%         for i =1:size(phi_n,1)
%             if all(w_val{1}*phi_n(i,:)'>=EPS*ones(n_plane,1))
%                     neg_mis=[neg_mis,i];
%                 end
% %             if phi_n(i,x1_idx)<=threshold_x
% %                 if all(w_val{1}*phi_n(i,:)'>=EPS*ones(n_plane,1))
% %                     neg_mis=[neg_mis,i];
% %                 end
% %             elseif phi_n(i,x1_idx)>threshold_x
% %                 if all(w_val{2}*phi_n(i,:)'>=EPS*ones(n_plane,1))
% %                     neg_mis=[neg_mis,i];
% %                 end
% %             end
%         end
                
%         if ~isempty(neg_mis)
%             disp('failed')
%             break
%         end
    else
        disp('SVM failed')
        break
%         phi_n = phi_negative;
    end
    
    iter=iter+1;
end
