function ball_rs = get_ball_rs(ball_r, num_affected, num_group)
    ball_rs = zeros(num_group, 1);
    ball_rs(1:num_affected) = ball_r;
    ball_rs(num_affected+1:end) = 2500*1e-6 + 250e-6*(rand(num_group-num_affected, 1)-0.5);
end