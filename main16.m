%16 states: 
% att(3);vel(3);pos(3);gyro bias(3);acc bias(3)
% distance bias(1)

clear
glvs

%data prepare
file='21131016-193240';
fid = fopen([file '.bin'], 'rb');
d = fread(fid, [13, inf], 'float32')';
fclose(fid);
plot(d(:,2:4))
be_used = 900:21000;
gyro_all = d(be_used,5:7)*pi/180;
gyro_all = gyro_all - mean(gyro_all(1:100,:));
acc_all = d(be_used,2:4)*9.8;
g0 = norm(mean(acc_all(1:100,:)));
ts = 0.005;
stance = detect_stance([gyro_all,acc_all]);
still = detect_still(stance);
len = length(gyro_all);
subplot(3,1,1); plot(0:ts:ts*(len-1),gyro_all);xlabel('time/s');ylabel('\omega /rad/s');legend('x','y','z');
subplot(3,1,2); plot(0:ts:ts*(len-1),acc_all);xlabel('time/s');ylabel('a /m/s^2');legend('x','y','z');
subplot(3,1,3); plot(0:ts:ts*(len-1),stance);ylim([0 2]);xlabel('time/s');legend('1 for stance');
distance_all = d(be_used,8)/100 ;

%imu init
acc0 = mean(acc_all(1:500,:));
acc0 = acc0/norm(acc0);
att(1) = atan2(acc0(2),acc0(3));    att(2) = -asin(acc0(1));    att(3) = 0;
Cnb0 = a2mat321(att);
qnb = m2qua(Cnb0);
        
vn = zeros(3,1); pn = zeros(3,1); gyro_b = zeros(3,1); acc_b = zeros(3,1); dis_b = 0;
pa = [-1.6317   -2.7582,0.9]';

%kf init
phi = [1;1;0]*glv.deg;
imuerr = imuerrset(500, 5000, 5, 500);
kf = avnkfinit(ts,  phi, imuerr);

[attk, xkpk, vnk, pnk, gbk, abk, disbk] = prealloc(len, 3, 2*kf.n, 3, 3, 3, 3, 1);
stance_yaw = 0;stance_yaw_num = 0;last_stance_yaw = 0;
%step by step
for k=1:len
    %imu update
    wvm = [gyro_all(k,:),acc_all(k,:)]*ts;
    wvm = wvm - [gyro_b',acc_b']*ts;
    Cnb = q2mat(qnb);
    dvn = Cnb*wvm(4:6)';
    vn = vn + dvn + [0;0;-g0]*ts;
    pn = pn + vn*ts;
    qnb = qupdt(qnb, wvm(1:3)');
        
    %kf update
    Cnbts = Cnb*ts;
    kf.Phikk_1(4:6,1:3) = askew(dvn);
    kf.Phikk_1(1:3,10:12) = -Cnbts; kf.Phikk_1(4:6,13:15) = Cnbts;
    kf.Phikk_1(7:9,4:6) = eye(3)*ts;
    kf = kfupdate(kf);
    % zupt
    if(stance(k))
        kf.Hk = [zeros(3),eye(3),zeros(3,10)];
        wvn = [0.01;0.01;0.01];
        kf.Rk = diag(wvn)^2;
        kf.adaptive = 0;
        kf = kfupdate(kf, vn, 'M');
        [kf,qnb,vn,pn,gyro_b,acc_b,dis_b] = kffeedback(kf,qnb,vn,pn,gyro_b,acc_b,dis_b);
    end
    %ZARU
    if(still(k))
        kf.Hk = [zeros(3,9),eye(3),zeros(3,4)];
        kf.Rk = diag([0.1;0.1;0.1]).^2;   
        kf.adaptive = 0;
        kf = kfupdate(kf, wvm(1:3)'/ts, 'M');
        [kf,qnb,vn,pn,gyro_b,acc_b,dis_b] = kffeedback(kf,qnb,vn,pn,gyro_b,acc_b,dis_b);
    end
    %pos_z = 0
    if(stance(k))
        kf.Hk = [zeros(1,8),1,zeros(1,7)];
        kf.Rk = 0.01^2;
        kf.adaptive = 0;
        kf = kfupdate(kf, pn(3), 'M');
        [kf,qnb,vn,pn,gyro_b,acc_b,dis_b] = kffeedback(kf,qnb,vn,pn,gyro_b,acc_b,dis_b);
    end
    %HDR 
    if(stance(k))
    % if(0)
        att = q2att321(qnb);
        d_psi = att(3) - last_stance_yaw;
        if(abs(d_psi) < 4*glv.deg)        
            kf.Hk = [0,0,1,zeros(1,13)];
            kf.Rk = 0.1^2;
            kf = kfupdate(kf, -d_psi, 'M');
            [kf,qnb,vn,pn,gyro_b,acc_b,dis_b] = kffeedback(kf,qnb,vn,pn,gyro_b,acc_b,dis_b);
        end
    end
    %distance measure
    if(stance(k))
    % if(0)
        pr = pn-pa;        
        distance_ins = norm(pr);
        kf.Hk = [zeros(1,6),pr'/distance_ins,zeros(1,6),-1];
        kf.Rk = 0.1^2;
        res = distance_ins - (distance_all(k) - dis_b);
        kf.adaptive = 0;kf.b = 0.9;kf.beta = 1;kf.Rmin = 0.01*kf.Rk;kf.Rmax = 100*kf.Rk;
        kf.m = 1;
        kf = kfupdate(kf, res, 'M');
        [kf,qnb,vn,pn,gyro_b,acc_b,dis_b] = kffeedback(kf,qnb,vn,pn,gyro_b,acc_b,dis_b);
    end

    attk(k,:) = q2att321(qnb)';    
    if(stance(k))
        stance_yaw = (stance_yaw*stance_yaw_num + attk(k,3))/(stance_yaw_num + 1);
        stance_yaw_num = stance_yaw_num + 1;
    else
        last_stance_yaw = stance_yaw;
        stance_yaw_num = 0;
    end

    xkpk(k,:) = [kf.xk; diag(kf.Pxk)];
    vnk(k,:) = vn;    pnk(k,:) = pn;    
    gbk(k,:) = gyro_b;  abk(k,:) = acc_b;
    disbk(k,:) = dis_b;    
end
avnplot(ts, attk, xkpk, vnk, pnk, gbk, abk, distance_all, disbk,pa);

function kf = avnkfinit(nts,  phi0, imuerr)
    kf = []; kf.s = 1; kf.nts = nts;
	kf.Qk = diag([imuerr.web; imuerr.wdb; zeros(10,1)])^2*nts;
    %kf.Qk = diag([imuerr.web; imuerr.wdb; zeros(10,1)])^2;
    kf.Gammak = 1;
	kf.Pxk = diag([phi0; 0*[1;1;1]; [0;0;0];imuerr.eb; imuerr.db; 0.1])^2;
	Ft = zeros(16); kf.Phikk_1 = eye(16)+Ft*nts;
	kf.n = 16;
    kf.I = eye(kf.n);
    kf.xk = zeros(kf.n, 1); 
    kf.adaptive = 0;
    kf.xconstrain = 0; kf.pconstrain = 0;
    kf.fading = 1;
end

function avnplot(ts, attk, xkpk, vnk, pnk, gbk, abk, distance_all, disbk,pa)
    global glv
    t = (1:length(attk))'*ts;
    myfigure('a');
	subplot(421); plot(t, attk(:,1:2)/glv.deg); xygo('pr')
	subplot(423); plot(t, attk(:,3)/glv.deg); xygo('y');
	subplot(425); plot(t, gbk/glv.dph); xygo('eb'); 
	subplot(427); plot(t, abk/glv.ug); xygo('db'); 
	subplot(422); plot(t, sqrt(xkpk(:,16:18))/glv.min); xygo('phi');
	subplot(424); plot(t, sqrt(xkpk(:,22:24))); xygo('dP');
	subplot(426); plot(t, sqrt(xkpk(:,25:27))/glv.dph); xygo('eb');
 	subplot(428); plot(t, sqrt(xkpk(:,28:30))/glv.ug); xygo('db');  
    myfigure('b');
    subplot(311); plot(t,vnk);
    subplot(312); plot(t,pnk);
    subplot(313); plot(t,distance_all);
    myfigure('c');    
    plot(pnk(:,1),pnk(:,2));axis equal
    hold on;
    plot(pnk(1,1),pnk(1,2),'<');plot(pnk(end,1),pnk(end,2),'s');
    plot(pa(1),pa(2),'p');
    xlabel('x/m');  ylabel('y/m');  
    myfigure('d');
    plot(t,disbk);
    xlabel('time/s');  ylabel('ranging bias/m');  
end

function stance = detect_stance(imu)
    global glv
    len = length(imu);
    acc_mag = zeros(len,1);
    gyro_mag = zeros(len,1);
    for i = 1:len
        acc_mag(i) = norm(imu(i,4:6));
        gyro_mag(i) = norm(imu(i,1:3));
    end
    acc_stance_threshold_H = 11;
    acc_stance_threshold_L = 9;
    gyro_stance_threshold = 50*glv.deg;

    stance_acc_H = (acc_mag < acc_stance_threshold_H);
    stance_acc_L = (acc_mag > acc_stance_threshold_L);
    stance_acc = stance_acc_H & stance_acc_L; %C1
    stance_gyro = (gyro_mag < gyro_stance_threshold); %C2

    stance = stance_acc & stance_gyro;

    % this window is necessary to clean stance array from false stance detection
    W = 20;
    for k = 1:len-W+1
        if (stance(k) == true) && (stance(k+W-1) == true)
            stance(k:k+W-1) = ones(W,1);
        end
    end
    for k = 1:len-W+1
        if (stance(k) == false) && (stance(k+W-1) == false)
            stance(k:k+W-1) = zeros(W,1);
        end
    end
end

function still = detect_still(stance)
len = length(stance);
still = false(len,1);
for i = 400:len
    if(stance(i-399:i) == true(400,1))
        still(i) = 1;
    end    
end
end

function [kf,qnb,vn,pn,gyro_b,acc_b,dis_b] = kffeedback(kf,qnb,vn,pn,gyro_b,acc_b,dis_b)
    qnb = qdelphi(qnb, 1*kf.xk(1:3)); kf.xk(1:3) = 0*kf.xk(1:3);
    vn = vn-1*kf.xk(4:6);  kf.xk(4:6) = 0*kf.xk(4:6);
    pn = pn-1*kf.xk(7:9);  kf.xk(7:9) = 0*kf.xk(7:9);
    gyro_b = gyro_b+1*kf.xk(10:12);  kf.xk(10:12) = 0*kf.xk(10:12);
    acc_b = acc_b+1*kf.xk(13:15);  kf.xk(13:15) = 0*kf.xk(13:15);
    dis_b = dis_b+1*kf.xk(16);  kf.xk(16) = 0*kf.xk(16);    
end
