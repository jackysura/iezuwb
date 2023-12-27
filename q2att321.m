function att = q2att321(qnb)
% Convert attitude quaternion to Euler attitude angles.
%
% Prototype: att = q2att(qnb)
% Input: qnb - attitude quaternion
% Output: att - Euler angles att=[pitch; roll; yaw] in radians
%
% See also  a2mat, a2qua, m2att, m2qua, q2mat, q2att1, attsyn, q2rv, incline.

% Copyright(c) 2009-2014, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 21/02/2008, 28/01/2013

     [att1,att] = m2att(q2mat(qnb));

    
    