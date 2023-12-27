%Cbn，本体坐标系到地理坐标系的转换矩阵，地理坐标系指东-磁北-天
%acc，静止状态下加速度计输出
%mag，静止状态下磁力计输出
function Cnb = init_mat(acc, mag)
mag = mag / norm(mag);
acc = acc / norm(acc);
e = cross(mag, acc);
e = e/norm(e);
n = cross(acc, e);
n = n/norm(n);

Cbn = [e(1) n(1) acc(1);
       e(2) n(2) acc(2);
       e(3) n(3) acc(3)];
Cnb = Cbn';
end