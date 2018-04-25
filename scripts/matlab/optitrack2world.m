function[output] = optitrack2world(input, transform, offset)
%transform data from Optitrack to the Simulation World Reference

%get simulation time 
%workaround: work out problem with the ros clock!!!
ros_time = input(:,1) - input(1,1);
dt = 1 / 230; %!!!
time = 0:dt:length(ros_time) * dt - dt;
ROS_CORRECT = 1;
TIME_TRANSFORM = 1;
if(~TIME_TRANSFORM)
    if(~ROS_CORRECT)
        output = time';
    else
        %or :
        dummy_vec = [ros_time(2:end);ros_time(end)];
        ros_dt = dummy_vec - ros_time;
        
        remove_ind = find(ros_dt < 5e-3);
        input(remove_ind,:) = [];
        ros_time = input(:,1) - input(1,1);
        output = ros_time;
    end
else
    %remove redundant points
    time = input(:,1);
    time_shifted = circshift(time,1);
    time_shifted(1) = time(1);
    
    time_dt = time - time_shifted;
    idx_to_remove = find(time_dt < 1e-3);
    input(idx_to_remove,:) = [];
    output = input(:,1);
end
%transform coordinates
world_pos = [];
world_att = [];
for i = 1:size(input, 1)
    pos = input(i,2:4);
    att = input(i,5:8);
  
    offset_b = quatmul(quatmul(att, [0 offset']), quatinv(att));
    pos = pos + offset_b(2:4);
    
    %transform to the world frame
    x_w = quatmul(quatmul(transform, [0 pos]), quatinv(transform));
    att_w = quatmul(quatmul(transform, att), quatinv(transform));
    
    world_pos = [world_pos;x_w(2:4)];
    world_att = [world_att;att_w];
end

output = [output, world_pos, world_att];