module G
function f(q::Vector{T}, q_dot::Vector{T}, xi::T) where T
    l1 = q[1]
    l2 = q[2]
    l3 = q[3]
    l1_dot = q_dot[1]
    l2_dot = q_dot[2]
    l3_dot = q_dot[3]
    
    [1700    0    0;
0 1700    0;
0    0 1700]
end
end