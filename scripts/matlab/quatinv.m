function [inv_q] = quatinv(q)
%inverse quaternion 

inv_q = [q(1) -q(2) -q(3) -q(4)] / norm(q);