function residual = object_fun(x)
    global step step_pos step_dis
    pa_pos = [x(1:2),0.9]';
    % pa_pos = x';
    % dis_b = x(3);
    residual = 0;
    for i = 1:step
        residual = residual + (norm(step_pos(:,i)-pa_pos)-(step_dis(i)))^2;
    end
    residual = residual/step;
end