function [q] = quatmul(q1,q2)
%quaternion multiplication
s1 = q1(1);
v1 = q1(2:4);

s2 = q2(1);
v2 = q2(2:4);

q = [s1*s2 - dot(v1,v2), cross(v1,v2) + s1 * v2 + s2 * v1];