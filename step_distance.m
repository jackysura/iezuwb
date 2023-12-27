%estimate base pos

global step step_pos step_dis

step = 0;
step_start = 1;
step_end = 0;
step_pos = zeros(3,100);
step_dis = zeros(1,100);
for k=1:len-1
    if(stance(k) == 1 && stance(k+1) == 0)
        step_end = k;
        step = step+1;
        step_pos(:,step) = mean(pnk(step_start:step_end,:));
        step_dis(step) = mean(distance_all(step_start:step_end));
    end
    if(stance(k) == 0 && stance(k+1) == 1)
        step_start = k;        
    end    
end

pa_xy = [1,0];
% dis_b = -0.25;
[x,feval] = fminunc('object_fun',[pa_xy])