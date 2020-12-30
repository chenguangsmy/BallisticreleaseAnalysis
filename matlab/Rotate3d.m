function [Data_R] = Rotate3d(Data, axis, degree)
%Rotate from original dimemsion to targeted dimension, 
% Data: original 3-by-n matrix;
% axis: could be 'x', 'y' or 'z'
% degree: the degree to rotate (90degree here 90)

%Transmation matrix:
% 'x': [1            0            0
%       0            cos(pitch)   -sin(pitch)
%       0            sin(pitch)   cos(pitch)]
%
% 'y': [cos(roll)    0            -sin(roll)
%       0            1            0
%       sin(coll)    0            cos(roll)]
%
% 'z': [cos(yaw)     sin(yaw)     0
%       -sin(yaw)    cos(yaw)     0
%       0            0            1]
%
%   Since force sensor has a rotate in setting up: x-backup, y-frontup,
%   z-left, we want it rotate to position we used to be: eg. x-right, y-front,
%   z-up. Thus, we need a linear transformation T1=('z', 45); T2=('x',-90)

theta = degree/180*pi; % the rotate degree;
% rotate direction: z-axis, clockwise;
switch (axis)
    case 'x',
        rot_M = ...
            [1,0,0;...
            0, cos(theta), -sin(theta);...
            0, sin(theta), cos(theta)];
    case 'y',
        rot_M = ...
            [cos(theta), 0, -sin(theta);...
             0, 1, 0; ...
             sin(theta), 0, cos(theta)];
    case 'z',
        rot_M = ...
            [cos(theta), -sin(theta), 0;...
             sin(theta), cos(theta), 0;...
             0, 0, 1];
end

Data_R = inv(rot_M)*Data;
end

