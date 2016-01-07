% Generate possible set of q_l
q_range = round(q_l_max/partial_q)+1;
for iVal = 0:1:q_range^N -1
    remain = iVal;
    for i=1:N   
        temp = floor(remain/q_range^(N-i));
        remain = remain - temp*q_range^(N-i);
        q_l(i) = temp*partial_q;
    end   
end