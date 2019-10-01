bezier_regressor
addpath FirstBatchData
listing = dir('FirstBatchData');
listing = listing(3:end);
types ={'Jaywalking','Crosswalking'};
trace_table = [];
Table_name = {'Device','Trip','t','t0','tf','t_clock','GPS_heading','Latitude','Longitude','Altitude','Distance','vx','Throttle','ax','ay',...
    'yaw','Brake','ACC_on','TurnSignal','NumTarget','TargetID','ObstacleID','Range','RangeRate','Transversal','TargetType','Status','CIPV',...
    'Warning_left','Warning_right','Boundary_type_left','Boundary_type_right','Dis_left','Dis_right','lane_heading','Lane_quality_left','Lane_quality_right'};
for i=1:length(listing)
    data = readtable(listing(i).name,'ReadVariableNames',false);
    data.Properties.VariableNames=Table_name;
    TripID = unique(data.Trip);
    Obstacle_ID_range = unique(data.ObstacleID);
    Target_ID_range = unique(data.TargetID);
    if contains(listing(i).name,types{1})
        trace_type = types{1};
    elseif contains(listing(i).name,types{2})
        trace_type = types{2};
    end
    for n=1:length(TripID)
        for j=1:length(Obstacle_ID_range)
            for k=1:length(Target_ID_range)
                data1 = data(data.Trip==TripID(n)&data.ObstacleID==Obstacle_ID_range(j)&data.TargetID==Target_ID_range(k)&data.Range<=20,:);
                t = data1.t;
                if length(t)<10
                    break
                end
                t_parse=1;
                for l=1:length(t)-1
                    if t(l+1)-t(l)>10000
                        t_parse = [t_parse l+1];
                    end
                end
                t_parse = [t_parse length(t)];
                for l=1:length(t_parse)-1
                    trace=[];
                    trace.data = data1(t_parse(l):t_parse(l+1)-1,:);
                    trace.trace_type = trace_type;
                    trace.target_type = trace.data.TargetType(1);
                    if length(trace.data.t)>10 && min(trace.data.Range)<15
                        trace_table = [trace_table,trace];
                    end
                end
            end
        end
    end
end
data_jw = [];
data_cw = [];

for i=1:length(trace_table)
    if trace_table(i).target_type==3
        n = length(trace_table(i).data.t);
        if n>=9
            t = (trace_table(i).data.t-trace_table(i).data.t(1)*ones(length(trace_table(i).data.t),1))/1000;
            T = (trace_table(i).data.t(end)-trace_table(i).data.t(1))/1000;
            ped_y0 = trace_table(i).data.Dis_right-trace_table(i).data.Transversal;
            veh_y = trace_table(i).data.Dis_right;
            ped_rel_x0 = sqrt(trace_table(i).data.Range.^2-trace_table(i).data.Transversal.^2);
            [ped_rel_x,ped_rel_vx] = bezier_regression(ped_rel_x0,t);
            [ped_y,ped_vy] = bezier_regression(ped_y0,t);
            ped_vx = trace_table(i).data.vx+ped_rel_vx;
            %         ped_vx = trace_table(i).data.vx+trace_table(i).data.RangeRate;
            idx=find(abs(ped_vx)<=3&abs(ped_vy)<=3);
            %         trace_table(i).data = trace_table(i).data(idx,:);
            
            
            
            if max(abs(ped_vx))>5
                disp('')
            end
            if max(abs(ped_vy))>5
                disp('')
            end
            switch trace_table(i).trace_type
                case 'Jaywalking'
                    data_jw = [data_jw;ped_rel_x(idx),ped_y(idx),veh_y(idx),trace_table(i).data.vx(idx),ped_vx(idx),ped_vy(idx),trace_table(i).data.ax(idx),trace_table(i).data.TurnSignal(idx)];
                case 'Crosswalking'
                    data_cw = [data_cw;ped_rel_x(idx),ped_y(idx),veh_y(idx),trace_table(i).data.vx(idx),ped_vx(idx),ped_vy(idx),trace_table(i).data.ax(idx),trace_table(i).data.TurnSignal(idx)];
            end
        end
        
    end
end