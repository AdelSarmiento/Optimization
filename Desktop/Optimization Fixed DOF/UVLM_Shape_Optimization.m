close all
clear all
clc

M = 6;  % # of chordwise panels
N = 10;  % # of spanwise panels

c = 1;  % chord
b = 6;  % span

V_ref = 10;  % U velocity, and reference for some terms

omega = 2*7.5;  % frequency, rad/s

t_step = 2*pi/omega/40;
N_steps = 50;  % number of steps
N_w = 10;  % number of deforming wake rows (only update the newest N_w wake rings)

t = 0:t_step:t_step*N_steps;  % time vector

N_Q_begin = 11; % Q (the cyclic forces) computed between N_Q_begin and N_steps+1

cutoff = 1E-5;  % cutoff radius where the Biot-Savart does not act

rho = 1.225;  % fluid density

%% Morphing parameters:

TWamp = 0;      % twist amplitude in degrees  (45 max)
TWphase = 10/180*pi;  % twist phase in radians (-pi to pi)

Bamp = 0;  % bending amplitude in meters (0.7 max)
Bphase = -88.888889/180*pi;  % bending phase in radians

%% Specify kinematics:

for i = 1:length(t)
    path_X(i) = -V_ref*t(i);   % forward velocity
    path_Y(i) = 0;
    path_Z(i) = 0*sin(omega*t(i));  % plunge
    roll(i) = 45*sin(omega*t(i));     % flap
    pitch(i) = 5;               % pitch
    yaw(i) = 0*sin(omega*t(i));
end

xo=load('GridxWR_P3.dat');
yo=load('GridyWR_P3.dat');
zo=load('GridzWR_P3.dat');

%pres = zeros(M,N);

%% Compute un-deformed wing grid in body-attached system:

[xo,yo] = meshgrid(0:c/M:c, 0:b/2/N:b/2);
xo = xo';
yo = yo';
zo = c*(1.8585*(xo/c).^6 - 5.4041*(xo/c).^5 + 5.336*(xo/c).^4 - 1.5454*(xo/c).^3 - 0.7805*(xo/c).^2 + 0.5355*(xo/c));

%% Compute morphing modal amplitudes:

twist = TWamp*sin(omega*t + TWphase);
bend = Bamp*sin(omega*t+Bphase);
fore = (-1/b)*bend.^2;

%% Compute morphed wing shape in body-attached system:

for i = 1:length(t)
    x(:,:,i) = xo - xo.*(1-cos((twist(i)*pi/180)*2*yo/b));
    y(:,:,i) = yo + fore(i)*4*(yo/b).^2;
    z(:,:,i) = zo - xo.*sin((twist(i)*pi/180)*2*yo/b) + bend(i)*4*(yo/b).^2;
end

for i = 1:length(t)
    
    %% Compute ring coordinates in body-attached system:
    for n = 1:N+1
        for m = 1:M
            xr(m,n,i) = (x(m+1,n,i)-x(m,n,i))/4 + x(m,n,i);
            yr(m,n,i) = y(m,n,i);
            zr(m,n,i) = (z(m+1,n,i)-z(m,n,i))/4 + z(m,n,i);
        end
        xr(m+1,n,i) = (x(m+1,n,i)-x(m,n,i))/4 + x(m+1,n,i);
        yr(m+1,n,i) = y(m,n,i);
        zr(m+1,n,i) = z(m+1,n,i);
    end
    
    %% Compute collocation points in body-attached system:
    for n = 1:N
        for m = 1:M
            xc(m,n,i) = (xr(m,n,i)+xr(m+1,n,i)+xr(m,n+1,i)+xr(m+1,n+1,i))/4;
            yc(m,n,i) = (yr(m,n,i)+yr(m+1,n,i)+yr(m,n+1,i)+yr(m+1,n+1,i))/4;
            zc(m,n,i) = (zr(m,n,i)+zr(m+1,n,i)+zr(m,n+1,i)+zr(m+1,n+1,i))/4;
        end
    end
    
    %% Transformation matrix between inertial and body-attached
    
    T_roll = [1,0,0;0,cos(roll(i)*pi/180),-sin(roll(i)*pi/180);0,sin(roll(i)*pi/180),cos(roll(i)*pi/180)];
    T_pitch = [cos(pitch(i)*pi/180),0,sin(pitch(i)*pi/180);0,1,0;-sin(pitch(i)*pi/180),0,cos(pitch(i)*pi/180)];
    T_yaw = [cos(yaw(i)*pi/180),-sin(yaw(i)*pi/180),0;sin(yaw(i)*pi/180),cos(yaw(i)*pi/180),0;0,0,1];
    T = T_yaw*T_pitch*T_roll;
    
    for m = 1:M+1
        for n = 1:N+1
            
            %% Compute ring coordinates in inertial system:
            foo = T*[xr(m,n,i);yr(m,n,i);zr(m,n,i)];
            if yr(m,n,i) == 0
                foo(2) = 0;  % prevents root of cambered wing from crossing to negative Y..
            end
            Xr(m,n,i) = foo(1)+path_X(i);
            Yr(m,n,i) = foo(2)+path_Y(i);
            Zr(m,n,i) = foo(3)+path_Z(i);
            
            %% Compute collocation points in inertial system:
            if m < M+1 && n < N+1
                foo = T*[xc(m,n,i);yc(m,n,i);zc(m,n,i)];
                Xc(m,n,i) = foo(1)+path_X(i);
                Yc(m,n,i) = foo(2)+path_Y(i);
                Zc(m,n,i) = foo(3)+path_Z(i);
            end
            
        end
    end
    
end

%% Time marching:

for i = 1:length(t)
    
    for m = 1:M
        for n = 1:N
            
            %% Compute velocities due to wing movement (forward finite difference):
            if i == 1
                U_kin(m,n,i) = 0;
                V_kin(m,n,i) = 0;
                W_kin(m,n,i) = 0;
            else
                U_kin(m,n,i) = -(Xc(m,n,i)-Xc(m,n,i-1))/t_step;
                V_kin(m,n,i) = -(Yc(m,n,i)-Yc(m,n,i-1))/t_step;
                W_kin(m,n,i) = -(Zc(m,n,i)-Zc(m,n,i-1))/t_step;
            end
            
            %% Outward normal vector:
            v1 = [Xr(m,n,i)-Xr(m+1,n+1,i),Yr(m,n,i)-Yr(m+1,n+1,i),Zr(m,n,i)-Zr(m+1,n+1,i)];
            v2 = [Xr(m+1,n,i)-Xr(m,n+1,i),Yr(m+1,n,i)-Yr(m,n+1,i),Zr(m+1,n,i)-Zr(m,n+1,i)];
            outward = [v1(2)*v2(3)-v2(2)*v1(3),v2(1)*v1(3)-v1(1)*v2(3),v1(1)*v2(2)-v2(1)*v1(2)];
            outward = outward / norm(outward);
            outward_X(m,n,i) = outward(1);
            outward_Y(m,n,i) = outward(2);
            outward_Z(m,n,i) = outward(3);
            normal_flow(m,n,i) = dot([outward(1),outward(2),outward(3)],[U_kin(m,n,i),V_kin(m,n,i),W_kin(m,n,i)]);
            
            %% tau vector (orthogonal to outward)
            tau_X(m,n,i) = .5*(Xr(m+1,n,i)+Xr(m+1,n+1,i))-.5*(Xr(m,n,i)+Xr(m,n+1,i));
            tau_Y(m,n,i) = .5*(Yr(m+1,n,i)+Yr(m+1,n+1,i))-.5*(Yr(m,n,i)+Yr(m,n+1,i));
            tau_Z(m,n,i) = .5*(Zr(m+1,n,i)+Zr(m+1,n+1,i))-.5*(Zr(m,n,i)+Zr(m,n+1,i));
            foo = sqrt(tau_X(m,n,i)^2+tau_Y(m,n,i)^2+tau_Z(m,n,i)^2);
            tau_X(m,n,i) = tau_X(m,n,i)/foo;
            tau_Y(m,n,i) = tau_Y(m,n,i)/foo;
            tau_Z(m,n,i) = tau_Z(m,n,i)/foo;
            
            %% Angle of attack:
            V_h = dot([U_kin(m,n,i),V_kin(m,n,i),W_kin(m,n,i)],[tau_X(m,n,i),tau_Y(m,n,i),tau_Z(m,n,i)]);
            V_v = dot([U_kin(m,n,i),V_kin(m,n,i),W_kin(m,n,i)],[outward_X(m,n,i),outward_Y(m,n,i),outward_Z(m,n,i)]);
            alpha(m,n,i) = atan2(V_v,V_h)*180/pi;
            
            %% Lift and drag vectors:
            third = cross([tau_X(m,n,i),tau_Y(m,n,i),tau_Z(m,n,i)],[outward_X(m,n,i),outward_Y(m,n,i),outward_Z(m,n,i)]);
            
            T_foo = [[tau_X(m,n,i);tau_Y(m,n,i);tau_Z(m,n,i)],[outward_X(m,n,i);outward_Y(m,n,i);outward_Z(m,n,i)],third'];
            T_alpha = [cos(alpha(m,n,i)*pi/180),sin(alpha(m,n,i)*pi/180),0;-sin(alpha(m,n,i)*pi/180),cos(alpha(m,n,i)*pi/180),0;0,0,1];
            
            j1 = T_foo*T_alpha'*T_foo'*[tau_X(m,n,i);tau_Y(m,n,i);tau_Z(m,n,i)];
            j2 = T_foo*T_alpha'*T_foo'*[outward_X(m,n,i);outward_Y(m,n,i);outward_Z(m,n,i)];
            
            n_lift_x(m,n,i) = j2(1);
            n_lift_y(m,n,i) = j2(2);
            n_lift_z(m,n,i) = j2(3);
            n_drag_x(m,n,i) = j1(1);
            n_drag_y(m,n,i) = j1(2);
            n_drag_z(m,n,i) = j1(3);
            
        end
    end
    
    
    %% Add on to wake geometry:
    
    Xwake(i,:) = Xr(M+1,:,i);
    Ywake(i,:) = Yr(M+1,:,i);
    Zwake(i,:) = Zr(M+1,:,i);
    
    %% Set circulation of newest wake rings to be the old TE circulation:
    
    if i > 1
        circ_wake(i-1,:) = circ(M,:,i-1);
    end
    
    %% Wing on wing influence matrix and RHS:
    
    K = 0;
    for m = 1:M
        for n = 1:N
            K = K+1;
            L = 0;
            for m1 = 1:M
                for n1 = 1:N
                    L = L+1;
                    
                    ro = [Xr(m1,n1+1,i)-Xr(m1,n1,i),Yr(m1,n1+1,i)-Yr(m1,n1,i),Zr(m1,n1+1,i)-Zr(m1,n1,i)];
                    r1 = [Xc(m,n,i)-Xr(m1,n1,i),Yc(m,n,i)-Yr(m1,n1,i),Zc(m,n,i)-Zr(m1,n1,i)];
                    r2 = [Xc(m,n,i)-Xr(m1,n1+1,i),Yc(m,n,i)-Yr(m1,n1+1,i),Zc(m,n,i)-Zr(m1,n1+1,i)];
                    aa = [r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                    bb = (r1(2)*r2(3)-r1(3)*r2(2))^2+(r1(3)*r2(1)-r1(1)*r2(3))^2+(r1(1)*r2(2)-r1(2)*r2(1))^2;
                    cc = [r1/sqrt(r1(1)^2+r1(2)^2+r1(3)^2)-r2/sqrt(r2(1)^2+r2(2)^2+r2(3)^2)];
                    UVW1 = (aa/bb)*(ro(1)*cc(1)+ro(2)*cc(2)+ro(3)*cc(3))/4/pi;
                    
                    ro = [Xr(m1+1,n1+1,i)-Xr(m1,n1+1,i),Yr(m1+1,n1+1,i)-Yr(m1,n1+1,i),Zr(m1+1,n1+1,i)-Zr(m1,n1+1,i)];
                    r1 = [Xc(m,n,i)-Xr(m1,n1+1,i),Yc(m,n,i)-Yr(m1,n1+1,i),Zc(m,n,i)-Zr(m1,n1+1,i)];
                    r2 = [Xc(m,n,i)-Xr(m1+1,n1+1,i),Yc(m,n,i)-Yr(m1+1,n1+1,i),Zc(m,n,i)-Zr(m1+1,n1+1,i)];
                    aa = [r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                    bb = (r1(2)*r2(3)-r1(3)*r2(2))^2+(r1(3)*r2(1)-r1(1)*r2(3))^2+(r1(1)*r2(2)-r1(2)*r2(1))^2;
                    cc = [r1/sqrt(r1(1)^2+r1(2)^2+r1(3)^2)-r2/sqrt(r2(1)^2+r2(2)^2+r2(3)^2)];
                    UVW2 = (aa/bb)*(ro(1)*cc(1)+ro(2)*cc(2)+ro(3)*cc(3))/4/pi;
                    
                    ro = [Xr(m1+1,n1,i)-Xr(m1+1,n1+1,i),Yr(m1+1,n1,i)-Yr(m1+1,n1+1,i),Zr(m1+1,n1,i)-Zr(m1+1,n1+1,i)];
                    r1 = [Xc(m,n,i)-Xr(m1+1,n1+1,i),Yc(m,n,i)-Yr(m1+1,n1+1,i),Zc(m,n,i)-Zr(m1+1,n1+1,i)];
                    r2 = [Xc(m,n,i)-Xr(m1+1,n1,i),Yc(m,n,i)-Yr(m1+1,n1,i),Zc(m,n,i)-Zr(m1+1,n1,i)];
                    aa = [r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                    bb = (r1(2)*r2(3)-r1(3)*r2(2))^2+(r1(3)*r2(1)-r1(1)*r2(3))^2+(r1(1)*r2(2)-r1(2)*r2(1))^2;
                    cc = [r1/sqrt(r1(1)^2+r1(2)^2+r1(3)^2)-r2/sqrt(r2(1)^2+r2(2)^2+r2(3)^2)];
                    UVW3 = (aa/bb)*(ro(1)*cc(1)+ro(2)*cc(2)+ro(3)*cc(3))/4/pi;
                    
                    ro = [Xr(m1,n1,i)-Xr(m1+1,n1,i),Yr(m1,n1,i)-Yr(m1+1,n1,i),Zr(m1,n1,i)-Zr(m1+1,n1,i)];
                    r1 = [Xc(m,n,i)-Xr(m1+1,n1,i),Yc(m,n,i)-Yr(m1+1,n1,i),Zc(m,n,i)-Zr(m1+1,n1,i)];
                    r2 = [Xc(m,n,i)-Xr(m1,n1,i),Yc(m,n,i)-Yr(m1,n1,i),Zc(m,n,i)-Zr(m1,n1,i)];
                    aa = [r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                    bb = (r1(2)*r2(3)-r1(3)*r2(2))^2+(r1(3)*r2(1)-r1(1)*r2(3))^2+(r1(1)*r2(2)-r1(2)*r2(1))^2;
                    cc = [r1/sqrt(r1(1)^2+r1(2)^2+r1(3)^2)-r2/sqrt(r2(1)^2+r2(2)^2+r2(3)^2)];
                    UVW4 = (aa/bb)*(ro(1)*cc(1)+ro(2)*cc(2)+ro(3)*cc(3))/4/pi;
                    
                    ro = [Xr(m1,n1+1,i)-Xr(m1,n1,i),Yr(m1,n1+1,i)-Yr(m1,n1,i),Zr(m1,n1+1,i)-Zr(m1,n1,i)];
                    r1 = [Xc(m,n,i)-Xr(m1,n1,i),-Yc(m,n,i)-Yr(m1,n1,i),Zc(m,n,i)-Zr(m1,n1,i)];
                    r2 = [Xc(m,n,i)-Xr(m1,n1+1,i),-Yc(m,n,i)-Yr(m1,n1+1,i),Zc(m,n,i)-Zr(m1,n1+1,i)];
                    aa = [r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                    bb = (r1(2)*r2(3)-r1(3)*r2(2))^2+(r1(3)*r2(1)-r1(1)*r2(3))^2+(r1(1)*r2(2)-r1(2)*r2(1))^2;
                    cc = [r1/sqrt(r1(1)^2+r1(2)^2+r1(3)^2)-r2/sqrt(r2(1)^2+r2(2)^2+r2(3)^2)];
                    UVW1r = (aa/bb)*(ro(1)*cc(1)+ro(2)*cc(2)+ro(3)*cc(3))/4/pi;
                    UVW1r(2) = -UVW1r(2);
                    
                    ro = [Xr(m1+1,n1+1,i)-Xr(m1,n1+1,i),Yr(m1+1,n1+1,i)-Yr(m1,n1+1,i),Zr(m1+1,n1+1,i)-Zr(m1,n1+1,i)];
                    r1 = [Xc(m,n,i)-Xr(m1,n1+1,i),-Yc(m,n,i)-Yr(m1,n1+1,i),Zc(m,n,i)-Zr(m1,n1+1,i)];
                    r2 = [Xc(m,n,i)-Xr(m1+1,n1+1,i),-Yc(m,n,i)-Yr(m1+1,n1+1,i),Zc(m,n,i)-Zr(m1+1,n1+1,i)];
                    aa = [r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                    bb = (r1(2)*r2(3)-r1(3)*r2(2))^2+(r1(3)*r2(1)-r1(1)*r2(3))^2+(r1(1)*r2(2)-r1(2)*r2(1))^2;
                    cc = [r1/sqrt(r1(1)^2+r1(2)^2+r1(3)^2)-r2/sqrt(r2(1)^2+r2(2)^2+r2(3)^2)];
                    UVW2r = (aa/bb)*(ro(1)*cc(1)+ro(2)*cc(2)+ro(3)*cc(3))/4/pi;
                    UVW2r(2) = -UVW2r(2);
                    
                    ro = [Xr(m1+1,n1,i)-Xr(m1+1,n1+1,i),Yr(m1+1,n1,i)-Yr(m1+1,n1+1,i),Zr(m1+1,n1,i)-Zr(m1+1,n1+1,i)];
                    r1 = [Xc(m,n,i)-Xr(m1+1,n1+1,i),-Yc(m,n,i)-Yr(m1+1,n1+1,i),Zc(m,n,i)-Zr(m1+1,n1+1,i)];
                    r2 = [Xc(m,n,i)-Xr(m1+1,n1,i),-Yc(m,n,i)-Yr(m1+1,n1,i),Zc(m,n,i)-Zr(m1+1,n1,i)];
                    aa = [r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                    bb = (r1(2)*r2(3)-r1(3)*r2(2))^2+(r1(3)*r2(1)-r1(1)*r2(3))^2+(r1(1)*r2(2)-r1(2)*r2(1))^2;
                    cc = [r1/sqrt(r1(1)^2+r1(2)^2+r1(3)^2)-r2/sqrt(r2(1)^2+r2(2)^2+r2(3)^2)];
                    UVW3r = (aa/bb)*(ro(1)*cc(1)+ro(2)*cc(2)+ro(3)*cc(3))/4/pi;
                    UVW3r(2) = -UVW3r(2);
                    
                    ro = [Xr(m1,n1,i)-Xr(m1+1,n1,i),Yr(m1,n1,i)-Yr(m1+1,n1,i),Zr(m1,n1,i)-Zr(m1+1,n1,i)];
                    r1 = [Xc(m,n,i)-Xr(m1+1,n1,i),-Yc(m,n,i)-Yr(m1+1,n1,i),Zc(m,n,i)-Zr(m1+1,n1,i)];
                    r2 = [Xc(m,n,i)-Xr(m1,n1,i),-Yc(m,n,i)-Yr(m1,n1,i),Zc(m,n,i)-Zr(m1,n1,i)];
                    aa = [r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                    bb = (r1(2)*r2(3)-r1(3)*r2(2))^2+(r1(3)*r2(1)-r1(1)*r2(3))^2+(r1(1)*r2(2)-r1(2)*r2(1))^2;
                    cc = [r1/sqrt(r1(1)^2+r1(2)^2+r1(3)^2)-r2/sqrt(r2(1)^2+r2(2)^2+r2(3)^2)];
                    UVW4r = (aa/bb)*(ro(1)*cc(1)+ro(2)*cc(2)+ro(3)*cc(3))/4/pi;
                    UVW4r(2) = -UVW4r(2);
                    
                    UVW = UVW1+UVW2+UVW3+UVW4+UVW1r+UVW2r+UVW3r+UVW4r;
                    UVW_stream = UVW2+UVW4+UVW2r+UVW4r;
                    UVW_stream_and_TE = UVW2+UVW3+UVW4+UVW2r+UVW3r+UVW4r;
                    
                    %% Wing on wing influence matrix:
                    
                    C(K,L,i) = UVW(1)*outward_X(m,n,i)+UVW(2)*outward_Y(m,n,i)+UVW(3)*outward_Z(m,n,i);
                    
                    %% Wing on wing influence matrix (just induced downwash):
                    
                    if m1 < M
                        B(K,L,i) = UVW_stream(1)*n_lift_x(m,n,i)+UVW_stream(2)*n_lift_y(m,n,i)+UVW_stream(3)*n_lift_z(m,n,i);
                    else
                        B(K,L,i) = UVW_stream_and_TE(1)*n_lift_x(m,n,i)+UVW_stream_and_TE(2)*n_lift_y(m,n,i)+UVW_stream_and_TE(3)*n_lift_z(m,n,i);
                    end
                    
                end
            end
            
            %% system RHS (normal flow at each collocation point):
            
            LL(K,i) = -normal_flow(m,n,i);
            
        end
    end
    
    %% Wake on wing RHS:
    
    if i > 1
        
        K = 0;
        for m = 1:M
            for n = 1:N
                K = K+1;
                L = 0;
                for m1 = 1:i-1
                    for n1 = 1:N
                        L = L+1;
                        
                        ro = [Xwake(m1+1,n1+1)-Xwake(m1+1,n1),Ywake(m1+1,n1+1)-Ywake(m1+1,n1),Zwake(m1+1,n1+1)-Zwake(m1+1,n1)];
                        r1 = [Xc(m,n,i)-Xwake(m1+1,n1),Yc(m,n,i)-Ywake(m1+1,n1),Zc(m,n,i)-Zwake(m1+1,n1)];
                        r2 = [Xc(m,n,i)-Xwake(m1+1,n1+1),Yc(m,n,i)-Ywake(m1+1,n1+1),Zc(m,n,i)-Zwake(m1+1,n1+1)];
                        aa = [r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                        bb = (r1(2)*r2(3)-r1(3)*r2(2))^2+(r1(3)*r2(1)-r1(1)*r2(3))^2+(r1(1)*r2(2)-r1(2)*r2(1))^2;
                        if bb > cutoff
                            cc = [r1/sqrt(r1(1)^2+r1(2)^2+r1(3)^2)-r2/sqrt(r2(1)^2+r2(2)^2+r2(3)^2)];
                            UVW1 = (aa/bb)*(ro(1)*cc(1)+ro(2)*cc(2)+ro(3)*cc(3))/4/pi;
                        else
                            UVW1 = [0,0,0];
                        end
                        
                        ro = [Xwake(m1,n1+1)-Xwake(m1+1,n1+1),Ywake(m1,n1+1)-Ywake(m1+1,n1+1),Zwake(m1,n1+1)-Zwake(m1+1,n1+1)];
                        r1 = [Xc(m,n,i)-Xwake(m1+1,n1+1),Yc(m,n,i)-Ywake(m1+1,n1+1),Zc(m,n,i)-Zwake(m1+1,n1+1)];
                        r2 = [Xc(m,n,i)-Xwake(m1,n1+1),Yc(m,n,i)-Ywake(m1,n1+1),Zc(m,n,i)-Zwake(m1,n1+1)];
                        aa = [r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                        bb = (r1(2)*r2(3)-r1(3)*r2(2))^2+(r1(3)*r2(1)-r1(1)*r2(3))^2+(r1(1)*r2(2)-r1(2)*r2(1))^2;
                        if bb > cutoff
                            cc = [r1/sqrt(r1(1)^2+r1(2)^2+r1(3)^2)-r2/sqrt(r2(1)^2+r2(2)^2+r2(3)^2)];
                            UVW2 = (aa/bb)*(ro(1)*cc(1)+ro(2)*cc(2)+ro(3)*cc(3))/4/pi;
                        else
                            UVW2 = [0,0,0];
                        end
                        
                        
                        ro = [Xwake(m1,n1)-Xwake(m1,n1+1),Ywake(m1,n1)-Ywake(m1,n1+1),Zwake(m1,n1)-Zwake(m1,n1+1)];
                        r1 = [Xc(m,n,i)-Xwake(m1,n1+1),Yc(m,n,i)-Ywake(m1,n1+1),Zc(m,n,i)-Zwake(m1,n1+1)];
                        r2 = [Xc(m,n,i)-Xwake(m1,n1),Yc(m,n,i)-Ywake(m1,n1),Zc(m,n,i)-Zwake(m1,n1)];
                        aa = [r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                        bb = (r1(2)*r2(3)-r1(3)*r2(2))^2+(r1(3)*r2(1)-r1(1)*r2(3))^2+(r1(1)*r2(2)-r1(2)*r2(1))^2;
                        if bb > cutoff
                            cc = [r1/sqrt(r1(1)^2+r1(2)^2+r1(3)^2)-r2/sqrt(r2(1)^2+r2(2)^2+r2(3)^2)];
                            UVW3 = (aa/bb)*(ro(1)*cc(1)+ro(2)*cc(2)+ro(3)*cc(3))/4/pi;
                        else
                            UVW3 = [0,0,0];
                        end
                        
                        ro = [Xwake(m1+1,n1)-Xwake(m1,n1),Ywake(m1+1,n1)-Ywake(m1,n1),Zwake(m1+1,n1)-Zwake(m1,n1)];
                        r1 = [Xc(m,n,i)-Xwake(m1,n1),Yc(m,n,i)-Ywake(m1,n1),Zc(m,n,i)-Zwake(m1,n1)];
                        r2 = [Xc(m,n,i)-Xwake(m1+1,n1),Yc(m,n,i)-Ywake(m1+1,n1),Zc(m,n,i)-Zwake(m1+1,n1)];
                        aa = [r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                        bb = (r1(2)*r2(3)-r1(3)*r2(2))^2+(r1(3)*r2(1)-r1(1)*r2(3))^2+(r1(1)*r2(2)-r1(2)*r2(1))^2;
                        if bb > cutoff
                            cc = [r1/sqrt(r1(1)^2+r1(2)^2+r1(3)^2)-r2/sqrt(r2(1)^2+r2(2)^2+r2(3)^2)];
                            UVW4 = (aa/bb)*(ro(1)*cc(1)+ro(2)*cc(2)+ro(3)*cc(3))/4/pi;
                        else
                            UVW4 = [0,0,0];
                        end
                        
                        ro = [Xwake(m1+1,n1+1)-Xwake(m1+1,n1),Ywake(m1+1,n1+1)-Ywake(m1+1,n1),Zwake(m1+1,n1+1)-Zwake(m1+1,n1)];
                        r1 = [Xc(m,n,i)-Xwake(m1+1,n1),-Yc(m,n,i)-Ywake(m1+1,n1),Zc(m,n,i)-Zwake(m1+1,n1)];
                        r2 = [Xc(m,n,i)-Xwake(m1+1,n1+1),-Yc(m,n,i)-Ywake(m1+1,n1+1),Zc(m,n,i)-Zwake(m1+1,n1+1)];
                        aa = [r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                        bb = (r1(2)*r2(3)-r1(3)*r2(2))^2+(r1(3)*r2(1)-r1(1)*r2(3))^2+(r1(1)*r2(2)-r1(2)*r2(1))^2;
                        if bb > cutoff
                            cc = [r1/sqrt(r1(1)^2+r1(2)^2+r1(3)^2)-r2/sqrt(r2(1)^2+r2(2)^2+r2(3)^2)];
                            UVW1r = (aa/bb)*(ro(1)*cc(1)+ro(2)*cc(2)+ro(3)*cc(3))/4/pi;
                            UVW1r(2) = -UVW1r(2);
                        else
                            UVW1r = [0,0,0];
                        end
                        
                        ro = [Xwake(m1,n1+1)-Xwake(m1+1,n1+1),Ywake(m1,n1+1)-Ywake(m1+1,n1+1),Zwake(m1,n1+1)-Zwake(m1+1,n1+1)];
                        r1 = [Xc(m,n,i)-Xwake(m1+1,n1+1),-Yc(m,n,i)-Ywake(m1+1,n1+1),Zc(m,n,i)-Zwake(m1+1,n1+1)];
                        r2 = [Xc(m,n,i)-Xwake(m1,n1+1),-Yc(m,n,i)-Ywake(m1,n1+1),Zc(m,n,i)-Zwake(m1,n1+1)];
                        aa = [r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                        bb = (r1(2)*r2(3)-r1(3)*r2(2))^2+(r1(3)*r2(1)-r1(1)*r2(3))^2+(r1(1)*r2(2)-r1(2)*r2(1))^2;
                        if bb > cutoff
                            cc = [r1/sqrt(r1(1)^2+r1(2)^2+r1(3)^2)-r2/sqrt(r2(1)^2+r2(2)^2+r2(3)^2)];
                            UVW2r = (aa/bb)*(ro(1)*cc(1)+ro(2)*cc(2)+ro(3)*cc(3))/4/pi;
                            UVW2r(2) = -UVW2r(2);
                        else
                            UVW2r = [0,0,0];
                        end
                        
                        ro = [Xwake(m1,n1)-Xwake(m1,n1+1),Ywake(m1,n1)-Ywake(m1,n1+1),Zwake(m1,n1)-Zwake(m1,n1+1)];
                        r1 = [Xc(m,n,i)-Xwake(m1,n1+1),-Yc(m,n,i)-Ywake(m1,n1+1),Zc(m,n,i)-Zwake(m1,n1+1)];
                        r2 = [Xc(m,n,i)-Xwake(m1,n1),-Yc(m,n,i)-Ywake(m1,n1),Zc(m,n,i)-Zwake(m1,n1)];
                        aa = [r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                        bb = (r1(2)*r2(3)-r1(3)*r2(2))^2+(r1(3)*r2(1)-r1(1)*r2(3))^2+(r1(1)*r2(2)-r1(2)*r2(1))^2;
                        if bb > cutoff
                            cc = [r1/sqrt(r1(1)^2+r1(2)^2+r1(3)^2)-r2/sqrt(r2(1)^2+r2(2)^2+r2(3)^2)];
                            UVW3r = (aa/bb)*(ro(1)*cc(1)+ro(2)*cc(2)+ro(3)*cc(3))/4/pi;
                            UVW3r(2) = -UVW3r(2);
                        else
                            UVW3r = [0,0,0];
                        end
                        
                        ro = [Xwake(m1+1,n1)-Xwake(m1,n1),Ywake(m1+1,n1)-Ywake(m1,n1),Zwake(m1+1,n1)-Zwake(m1,n1)];
                        r1 = [Xc(m,n,i)-Xwake(m1,n1),-Yc(m,n,i)-Ywake(m1,n1),Zc(m,n,i)-Zwake(m1,n1)];
                        r2 = [Xc(m,n,i)-Xwake(m1+1,n1),-Yc(m,n,i)-Ywake(m1+1,n1),Zc(m,n,i)-Zwake(m1+1,n1)];
                        aa = [r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                        bb = (r1(2)*r2(3)-r1(3)*r2(2))^2+(r1(3)*r2(1)-r1(1)*r2(3))^2+(r1(1)*r2(2)-r1(2)*r2(1))^2;
                        if bb > cutoff
                            cc = [r1/sqrt(r1(1)^2+r1(2)^2+r1(3)^2)-r2/sqrt(r2(1)^2+r2(2)^2+r2(3)^2)];
                            UVW4r = (aa/bb)*(ro(1)*cc(1)+ro(2)*cc(2)+ro(3)*cc(3))/4/pi;
                            UVW4r(2) = -UVW4r(2);
                        else
                            UVW4r = [0,0,0];
                        end
                        
                        UVW = UVW1+UVW2+UVW3+UVW4+UVW1r+UVW2r+UVW3r+UVW4r;
                        UVW_stream = UVW2+UVW4+UVW2r+UVW4r;
                        
                        %% Wake on wing influence matrix:
                        
                        C_wake(K,L) = UVW(1)*outward_X(m,n,i)+UVW(2)*outward_Y(m,n,i)+UVW(3)*outward_Z(m,n,i);
                        
                        %% Wake on wing influence matrix (just induced downwash):
                        
                        B_wake(K,L) = UVW(1)*n_lift_x(m,n,i)+UVW(2)*n_lift_y(m,n,i)+UVW(3)*n_lift_z(m,n,i);
                        
                    end
                end
                
            end
        end
    end
    
    %% Solve system:
    
    if i == 1
        foo = inv(C(:,:,i))*LL(:,i);
        w_ind(:,:,i) = reshape(B(:,:,i)*foo,N,M)';  % induced downwash
    else
        foo = inv(C(:,:,i))*(LL(:,i)-C_wake*reshape(circ_wake',(i-1)*N,1));
        w_ind(:,:,i) = reshape(B(:,:,i)*foo,N,M)' + reshape(B_wake*reshape(circ_wake',(i-1)*N,1),N,M)';  % induced downwash
    end
    
    circ(:,:,i) = reshape(foo,N,M)';
    foo = circ(:,:,i);
    foo(2:M,:) = circ(2:M,:,i)-circ(1:M-1,:,i);
    circ_local(:,:,i) = foo;
    
    %% Compute forces, power, pressure:
    
    sigma1(1,:,i) = 0.5*circ(1,:,i);
    sigma1(2:M,:,i) = 0.5*(circ(1:M-1,:,i)+circ(2:M,:,i));
    if i == 1
        dfdt(:,:,i) = sigma1(:,:,i)/t_step;
    else
        dfdt(:,:,i) = (sigma1(:,:,i) - sigma1(:,:,i-1))/t_step;
    end
    
    power(1,i) = 0;
    for m = 1:M
        for n = 1:N
            
            chord(m,n) = sqrt((.5*(x(m+1,n,i)+x(m+1,n+1,i))-.5*(x(m,n,i)+x(m,n+1,i)))^2+(.5*(y(m+1,n,i)+y(m+1,n+1,i))-.5*(y(m,n,i)+y(m,n+1,i)))^2+(.5*(z(m+1,n,i)+z(m+1,n+1,i))-.5*(z(m,n,i)+z(m,n+1,i)))^2);
            span(m,n) = sqrt((.5*(x(m,n+1,i)+x(m+1,n+1,i))-.5*(x(m,n,i)+x(m+1,n,i)))^2+(.5*(y(m,n+1,i)+y(m+1,n+1,i))-.5*(y(m,n,i)+y(m+1,n,i)))^2+(.5*(z(m,n+1,i)+z(m+1,n+1,i))-.5*(z(m,n,i)+z(m+1,n,i)))^2);
            
            dL = rho*span(m,n)*(sqrt(U_kin(m,n,i)^2+V_kin(m,n,i)^2+W_kin(m,n,i)^2)*circ_local(m,n,i)+chord(m,n)*dfdt(m,n,i))*cos(alpha(m,n,i)*pi/180);
            dD = rho*span(m,n)*(-w_ind(m,n,i)*circ_local(m,n,i)+chord(m,n)*dfdt(m,n,i)*sin(alpha(m,n,i)*pi/180));
            TSec(m,n,i) = - dD;

            dL_vec = dL*[n_lift_x(m,n,i),n_lift_y(m,n,i),n_lift_z(m,n,i)];
            dD_vec = dD*[n_drag_x(m,n,i),n_drag_y(m,n,i),n_drag_z(m,n,i)];
 
            F_X(m,n,i) = dL_vec(1)+dD_vec(1);
            F_Y(m,n,i) = dL_vec(2)+dD_vec(2);
            F_Z(m,n,i) = dL_vec(3)+dD_vec(3);

            pressure = dot([outward_X(m,n,i),outward_Y(m,n,i),outward_Z(m,n,i)],[F_X(m,n,i),F_Y(m,n,i),F_Z(m,n,i)])/chord(m,n)/span(m,n);
            pres(m,n,i) = pressure;
            power(1,i) = power(1,i) + 2*pressure*dot([U_kin(m,n,i),V_kin(m,n,i),W_kin(m,n,i)],[outward_X(m,n,i),outward_Y(m,n,i),outward_Z(m,n,i)])*chord(m,n)*span(m,n);
            
        end
    end

    CD(i) = sum(sum(F_X(:,:,i)*2))/(b*c*.5*1.225*V_ref^2);
    CL(i) = sum(sum(F_Z(:,:,i)*2))/(b*c*.5*1.225*V_ref^2);

    %% Wing on wake rollup:
    
    if i > 1
        if i-1 < N_w
            wake_start = 1;
        else
            wake_start = i-N_w;
        end
        
        K = 0;
        for m = wake_start:i-1
            for n = 1:N+1
                
                K = K+1;
                L = 0;
                for m1 = 1:M
                    for n1 = 1:N
                        L = L+1;
                        
                        ro = [Xr(m1,n1+1,i)-Xr(m1,n1,i),Yr(m1,n1+1,i)-Yr(m1,n1,i),Zr(m1,n1+1,i)-Zr(m1,n1,i)];
                        r1 = [Xwake(m,n)-Xr(m1,n1,i),Ywake(m,n)-Yr(m1,n1,i),Zwake(m,n)-Zr(m1,n1,i)];
                        r2 = [Xwake(m,n)-Xr(m1,n1+1,i),Ywake(m,n)-Yr(m1,n1+1,i),Zwake(m,n)-Zr(m1,n1+1,i)];
                        aa = [r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                        bb = (r1(2)*r2(3)-r1(3)*r2(2))^2+(r1(3)*r2(1)-r1(1)*r2(3))^2+(r1(1)*r2(2)-r1(2)*r2(1))^2;
                        if bb > cutoff
                            cc = [r1/sqrt(r1(1)^2+r1(2)^2+r1(3)^2)-r2/sqrt(r2(1)^2+r2(2)^2+r2(3)^2)];
                            UVW1 = (aa/bb)*(ro(1)*cc(1)+ro(2)*cc(2)+ro(3)*cc(3))/4/pi;
                        else
                            UVW1 = [0,0,0];
                        end
                        
                        ro = [Xr(m1+1,n1+1,i)-Xr(m1,n1+1,i),Yr(m1+1,n1+1,i)-Yr(m1,n1+1,i),Zr(m1+1,n1+1,i)-Zr(m1,n1+1,i)];
                        r1 = [Xwake(m,n)-Xr(m1,n1+1,i),Ywake(m,n)-Yr(m1,n1+1,i),Zwake(m,n)-Zr(m1,n1+1,i)];
                        r2 = [Xwake(m,n)-Xr(m1+1,n1+1,i),Ywake(m,n)-Yr(m1+1,n1+1,i),Zwake(m,n)-Zr(m1+1,n1+1,i)];
                        aa = [r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                        bb = (r1(2)*r2(3)-r1(3)*r2(2))^2+(r1(3)*r2(1)-r1(1)*r2(3))^2+(r1(1)*r2(2)-r1(2)*r2(1))^2;
                        if bb > cutoff
                            cc = [r1/sqrt(r1(1)^2+r1(2)^2+r1(3)^2)-r2/sqrt(r2(1)^2+r2(2)^2+r2(3)^2)];
                            UVW2 = (aa/bb)*(ro(1)*cc(1)+ro(2)*cc(2)+ro(3)*cc(3))/4/pi;
                        else
                            UVW2 = [0,0,0];
                        end
                        
                        ro = [Xr(m1+1,n1,i)-Xr(m1+1,n1+1,i),Yr(m1+1,n1,i)-Yr(m1+1,n1+1,i),Zr(m1+1,n1,i)-Zr(m1+1,n1+1,i)];
                        r1 = [Xwake(m,n)-Xr(m1+1,n1+1,i),Ywake(m,n)-Yr(m1+1,n1+1,i),Zwake(m,n)-Zr(m1+1,n1+1,i)];
                        r2 = [Xwake(m,n)-Xr(m1+1,n1,i),Ywake(m,n)-Yr(m1+1,n1,i),Zwake(m,n)-Zr(m1+1,n1,i)];
                        aa = [r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                        bb = (r1(2)*r2(3)-r1(3)*r2(2))^2+(r1(3)*r2(1)-r1(1)*r2(3))^2+(r1(1)*r2(2)-r1(2)*r2(1))^2;
                        if bb > cutoff
                            cc = [r1/sqrt(r1(1)^2+r1(2)^2+r1(3)^2)-r2/sqrt(r2(1)^2+r2(2)^2+r2(3)^2)];
                            UVW3 = (aa/bb)*(ro(1)*cc(1)+ro(2)*cc(2)+ro(3)*cc(3))/4/pi;
                        else
                            UVW3 = [0,0,0];
                        end
                        
                        ro = [Xr(m1,n1,i)-Xr(m1+1,n1,i),Yr(m1,n1,i)-Yr(m1+1,n1,i),Zr(m1,n1,i)-Zr(m1+1,n1,i)];
                        r1 = [Xwake(m,n)-Xr(m1+1,n1,i),Ywake(m,n)-Yr(m1+1,n1,i),Zwake(m,n)-Zr(m1+1,n1,i)];
                        r2 = [Xwake(m,n)-Xr(m1,n1,i),Ywake(m,n)-Yr(m1,n1,i),Zwake(m,n)-Zr(m1,n1,i)];
                        aa = [r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                        bb = (r1(2)*r2(3)-r1(3)*r2(2))^2+(r1(3)*r2(1)-r1(1)*r2(3))^2+(r1(1)*r2(2)-r1(2)*r2(1))^2;
                        if bb > cutoff
                            cc = [r1/sqrt(r1(1)^2+r1(2)^2+r1(3)^2)-r2/sqrt(r2(1)^2+r2(2)^2+r2(3)^2)];
                            UVW4 = (aa/bb)*(ro(1)*cc(1)+ro(2)*cc(2)+ro(3)*cc(3))/4/pi;
                        else
                            UVW4 = [0,0,0];
                        end
                        
                        ro = [Xr(m1,n1+1,i)-Xr(m1,n1,i),Yr(m1,n1+1,i)-Yr(m1,n1,i),Zr(m1,n1+1,i)-Zr(m1,n1,i)];
                        r1 = [Xwake(m,n)-Xr(m1,n1,i),-Ywake(m,n)-Yr(m1,n1,i),Zwake(m,n)-Zr(m1,n1,i)];
                        r2 = [Xwake(m,n)-Xr(m1,n1+1,i),-Ywake(m,n)-Yr(m1,n1+1,i),Zwake(m,n)-Zr(m1,n1+1,i)];
                        aa = [r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                        bb = (r1(2)*r2(3)-r1(3)*r2(2))^2+(r1(3)*r2(1)-r1(1)*r2(3))^2+(r1(1)*r2(2)-r1(2)*r2(1))^2;
                        if bb > cutoff
                            cc = [r1/sqrt(r1(1)^2+r1(2)^2+r1(3)^2)-r2/sqrt(r2(1)^2+r2(2)^2+r2(3)^2)];
                            UVW1r = (aa/bb)*(ro(1)*cc(1)+ro(2)*cc(2)+ro(3)*cc(3))/4/pi;
                            UVW1r(2) = -UVW1r(2);
                        else
                            UVW1r = [0,0,0];
                        end
                        
                        
                        ro = [Xr(m1+1,n1+1,i)-Xr(m1,n1+1,i),Yr(m1+1,n1+1,i)-Yr(m1,n1+1,i),Zr(m1+1,n1+1,i)-Zr(m1,n1+1,i)];
                        r1 = [Xwake(m,n)-Xr(m1,n1+1,i),-Ywake(m,n)-Yr(m1,n1+1,i),Zwake(m,n)-Zr(m1,n1+1,i)];
                        r2 = [Xwake(m,n)-Xr(m1+1,n1+1,i),-Ywake(m,n)-Yr(m1+1,n1+1,i),Zwake(m,n)-Zr(m1+1,n1+1,i)];
                        aa = [r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                        bb = (r1(2)*r2(3)-r1(3)*r2(2))^2+(r1(3)*r2(1)-r1(1)*r2(3))^2+(r1(1)*r2(2)-r1(2)*r2(1))^2;
                        if bb > cutoff
                            cc = [r1/sqrt(r1(1)^2+r1(2)^2+r1(3)^2)-r2/sqrt(r2(1)^2+r2(2)^2+r2(3)^2)];
                            UVW2r = (aa/bb)*(ro(1)*cc(1)+ro(2)*cc(2)+ro(3)*cc(3))/4/pi;
                            UVW2r(2) = -UVW2r(2);
                        else
                            UVW2r = [0,0,0];
                        end
                        
                        ro = [Xr(m1+1,n1,i)-Xr(m1+1,n1+1,i),Yr(m1+1,n1,i)-Yr(m1+1,n1+1,i),Zr(m1+1,n1,i)-Zr(m1+1,n1+1,i)];
                        r1 = [Xwake(m,n)-Xr(m1+1,n1+1,i),-Ywake(m,n)-Yr(m1+1,n1+1,i),Zwake(m,n)-Zr(m1+1,n1+1,i)];
                        r2 = [Xwake(m,n)-Xr(m1+1,n1,i),-Ywake(m,n)-Yr(m1+1,n1,i),Zwake(m,n)-Zr(m1+1,n1,i)];
                        aa = [r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                        bb = (r1(2)*r2(3)-r1(3)*r2(2))^2+(r1(3)*r2(1)-r1(1)*r2(3))^2+(r1(1)*r2(2)-r1(2)*r2(1))^2;
                        if bb > cutoff
                            cc = [r1/sqrt(r1(1)^2+r1(2)^2+r1(3)^2)-r2/sqrt(r2(1)^2+r2(2)^2+r2(3)^2)];
                            UVW3r = (aa/bb)*(ro(1)*cc(1)+ro(2)*cc(2)+ro(3)*cc(3))/4/pi;
                            UVW3r(2) = -UVW3r(2);
                        else
                            UVW3r = [0,0,0];
                        end
                        
                        ro = [Xr(m1,n1,i)-Xr(m1+1,n1,i),Yr(m1,n1,i)-Yr(m1+1,n1,i),Zr(m1,n1,i)-Zr(m1+1,n1,i)];
                        r1 = [Xwake(m,n)-Xr(m1+1,n1,i),-Ywake(m,n)-Yr(m1+1,n1,i),Zwake(m,n)-Zr(m1+1,n1,i)];
                        r2 = [Xwake(m,n)-Xr(m1,n1,i),-Ywake(m,n)-Yr(m1,n1,i),Zwake(m,n)-Zr(m1,n1,i)];
                        aa = [r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                        bb = (r1(2)*r2(3)-r1(3)*r2(2))^2+(r1(3)*r2(1)-r1(1)*r2(3))^2+(r1(1)*r2(2)-r1(2)*r2(1))^2;
                        if bb > cutoff
                            cc = [r1/sqrt(r1(1)^2+r1(2)^2+r1(3)^2)-r2/sqrt(r2(1)^2+r2(2)^2+r2(3)^2)];
                            UVW4r = (aa/bb)*(ro(1)*cc(1)+ro(2)*cc(2)+ro(3)*cc(3))/4/pi;
                            UVW4r(2) = -UVW4r(2);
                        else
                            UVW4r = [0,0,0];
                        end
                        
                        UVW = UVW1+UVW2+UVW3+UVW4+UVW1r+UVW2r+UVW3r+UVW4r;
                        
                        %% Wing on wake influence matrix:
                        
                        C1_rollup(K,L,1) = UVW(1);
                        C1_rollup(K,L,2) = UVW(2);
                        C1_rollup(K,L,3) = UVW(3);
                        
                    end
                end
            end
        end
        
        %% Wake on wake rollup:
        
        K = 0;
        for m = wake_start:i-1
            for n = 1:N+1
                
                K = K+1;
                L = 0;
                for m1 = 1:i-1
                    for n1 = 1:N
                        L = L+1;
                        
                        ro = [Xwake(m1+1,n1+1)-Xwake(m1+1,n1),Ywake(m1+1,n1+1)-Ywake(m1+1,n1),Zwake(m1+1,n1+1)-Zwake(m1+1,n1)];
                        r1 = [Xwake(m,n)-Xwake(m1+1,n1),Ywake(m,n)-Ywake(m1+1,n1),Zwake(m,n)-Zwake(m1+1,n1)];
                        r2 = [Xwake(m,n)-Xwake(m1+1,n1+1),Ywake(m,n)-Ywake(m1+1,n1+1),Zwake(m,n)-Zwake(m1+1,n1+1)];
                        aa = [r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                        bb = (r1(2)*r2(3)-r1(3)*r2(2))^2+(r1(3)*r2(1)-r1(1)*r2(3))^2+(r1(1)*r2(2)-r1(2)*r2(1))^2;
                        if bb > cutoff
                            cc = [r1/sqrt(r1(1)^2+r1(2)^2+r1(3)^2)-r2/sqrt(r2(1)^2+r2(2)^2+r2(3)^2)];
                            UVW1 = (aa/bb)*(ro(1)*cc(1)+ro(2)*cc(2)+ro(3)*cc(3))/4/pi;
                        else
                            UVW1 = [0,0,0];
                        end
                        
                        ro = [Xwake(m1,n1+1)-Xwake(m1+1,n1+1),Ywake(m1,n1+1)-Ywake(m1+1,n1+1),Zwake(m1,n1+1)-Zwake(m1+1,n1+1)];
                        r1 = [Xwake(m,n)-Xwake(m1+1,n1+1),Ywake(m,n)-Ywake(m1+1,n1+1),Zwake(m,n)-Zwake(m1+1,n1+1)];
                        r2 = [Xwake(m,n)-Xwake(m1,n1+1),Ywake(m,n)-Ywake(m1,n1+1),Zwake(m,n)-Zwake(m1,n1+1)];
                        aa = [r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                        bb = (r1(2)*r2(3)-r1(3)*r2(2))^2+(r1(3)*r2(1)-r1(1)*r2(3))^2+(r1(1)*r2(2)-r1(2)*r2(1))^2;
                        if bb > cutoff
                            cc = [r1/sqrt(r1(1)^2+r1(2)^2+r1(3)^2)-r2/sqrt(r2(1)^2+r2(2)^2+r2(3)^2)];
                            UVW2 = (aa/bb)*(ro(1)*cc(1)+ro(2)*cc(2)+ro(3)*cc(3))/4/pi;
                        else
                            UVW2 = [0,0,0];
                        end
                        
                        ro = [Xwake(m1,n1)-Xwake(m1,n1+1),Ywake(m1,n1)-Ywake(m1,n1+1),Zwake(m1,n1)-Zwake(m1,n1+1)];
                        r1 = [Xwake(m,n)-Xwake(m1,n1+1),Ywake(m,n)-Ywake(m1,n1+1),Zwake(m,n)-Zwake(m1,n1+1)];
                        r2 = [Xwake(m,n)-Xwake(m1,n1),Ywake(m,n)-Ywake(m1,n1),Zwake(m,n)-Zwake(m1,n1)];
                        aa = [r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                        bb = (r1(2)*r2(3)-r1(3)*r2(2))^2+(r1(3)*r2(1)-r1(1)*r2(3))^2+(r1(1)*r2(2)-r1(2)*r2(1))^2;
                        if bb > cutoff
                            cc = [r1/sqrt(r1(1)^2+r1(2)^2+r1(3)^2)-r2/sqrt(r2(1)^2+r2(2)^2+r2(3)^2)];
                            UVW3 = (aa/bb)*(ro(1)*cc(1)+ro(2)*cc(2)+ro(3)*cc(3))/4/pi;
                        else
                            UVW3 = [0,0,0];
                        end
                        
                        ro = [Xwake(m1+1,n1)-Xwake(m1,n1),Ywake(m1+1,n1)-Ywake(m1,n1),Zwake(m1+1,n1)-Zwake(m1,n1)];
                        r1 = [Xwake(m,n)-Xwake(m1,n1),Ywake(m,n)-Ywake(m1,n1),Zwake(m,n)-Zwake(m1,n1)];
                        r2 = [Xwake(m,n)-Xwake(m1+1,n1),Ywake(m,n)-Ywake(m1+1,n1),Zwake(m,n)-Zwake(m1+1,n1)];
                        aa = [r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                        bb = (r1(2)*r2(3)-r1(3)*r2(2))^2+(r1(3)*r2(1)-r1(1)*r2(3))^2+(r1(1)*r2(2)-r1(2)*r2(1))^2;
                        if bb > cutoff
                            cc = [r1/sqrt(r1(1)^2+r1(2)^2+r1(3)^2)-r2/sqrt(r2(1)^2+r2(2)^2+r2(3)^2)];
                            UVW4 = (aa/bb)*(ro(1)*cc(1)+ro(2)*cc(2)+ro(3)*cc(3))/4/pi;
                        else
                            UVW4 = [0,0,0];
                        end
                        
                        ro = [Xwake(m1+1,n1+1)-Xwake(m1+1,n1),Ywake(m1+1,n1+1)-Ywake(m1+1,n1),Zwake(m1+1,n1+1)-Zwake(m1+1,n1)];
                        r1 = [Xwake(m,n)-Xwake(m1+1,n1),-Ywake(m,n)-Ywake(m1+1,n1),Zwake(m,n)-Zwake(m1+1,n1)];
                        r2 = [Xwake(m,n)-Xwake(m1+1,n1+1),-Ywake(m,n)-Ywake(m1+1,n1+1),Zwake(m,n)-Zwake(m1+1,n1+1)];
                        aa = [r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                        bb = (r1(2)*r2(3)-r1(3)*r2(2))^2+(r1(3)*r2(1)-r1(1)*r2(3))^2+(r1(1)*r2(2)-r1(2)*r2(1))^2;
                        if bb > cutoff
                            cc = [r1/sqrt(r1(1)^2+r1(2)^2+r1(3)^2)-r2/sqrt(r2(1)^2+r2(2)^2+r2(3)^2)];
                            UVW1r = (aa/bb)*(ro(1)*cc(1)+ro(2)*cc(2)+ro(3)*cc(3))/4/pi;
                            UVW1r(2) = -UVW1r(2);
                        else
                            UVW1r = [0,0,0];
                        end
                        
                        ro = [Xwake(m1,n1+1)-Xwake(m1+1,n1+1),Ywake(m1,n1+1)-Ywake(m1+1,n1+1),Zwake(m1,n1+1)-Zwake(m1+1,n1+1)];
                        r1 = [Xwake(m,n)-Xwake(m1+1,n1+1),-Ywake(m,n)-Ywake(m1+1,n1+1),Zwake(m,n)-Zwake(m1+1,n1+1)];
                        r2 = [Xwake(m,n)-Xwake(m1,n1+1),-Ywake(m,n)-Ywake(m1,n1+1),Zwake(m,n)-Zwake(m1,n1+1)];
                        aa = [r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                        bb = (r1(2)*r2(3)-r1(3)*r2(2))^2+(r1(3)*r2(1)-r1(1)*r2(3))^2+(r1(1)*r2(2)-r1(2)*r2(1))^2;
                        if bb > cutoff
                            cc = [r1/sqrt(r1(1)^2+r1(2)^2+r1(3)^2)-r2/sqrt(r2(1)^2+r2(2)^2+r2(3)^2)];
                            UVW2r = (aa/bb)*(ro(1)*cc(1)+ro(2)*cc(2)+ro(3)*cc(3))/4/pi;
                            UVW2r(2) = -UVW2r(2);
                        else
                            UVW2r = [0,0,0];
                        end
                        
                        ro = [Xwake(m1,n1)-Xwake(m1,n1+1),Ywake(m1,n1)-Ywake(m1,n1+1),Zwake(m1,n1)-Zwake(m1,n1+1)];
                        r1 = [Xwake(m,n)-Xwake(m1,n1+1),-Ywake(m,n)-Ywake(m1,n1+1),Zwake(m,n)-Zwake(m1,n1+1)];
                        r2 = [Xwake(m,n)-Xwake(m1,n1),-Ywake(m,n)-Ywake(m1,n1),Zwake(m,n)-Zwake(m1,n1)];
                        aa = [r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                        bb = (r1(2)*r2(3)-r1(3)*r2(2))^2+(r1(3)*r2(1)-r1(1)*r2(3))^2+(r1(1)*r2(2)-r1(2)*r2(1))^2;
                        if bb > cutoff
                            cc = [r1/sqrt(r1(1)^2+r1(2)^2+r1(3)^2)-r2/sqrt(r2(1)^2+r2(2)^2+r2(3)^2)];
                            UVW3r = (aa/bb)*(ro(1)*cc(1)+ro(2)*cc(2)+ro(3)*cc(3))/4/pi;
                            UVW3r(2) = -UVW3r(2);
                        else
                            UVW3r = [0,0,0];
                        end
                        
                        ro = [Xwake(m1+1,n1)-Xwake(m1,n1),Ywake(m1+1,n1)-Ywake(m1,n1),Zwake(m1+1,n1)-Zwake(m1,n1)];
                        r1 = [Xwake(m,n)-Xwake(m1,n1),-Ywake(m,n)-Ywake(m1,n1),Zwake(m,n)-Zwake(m1,n1)];
                        r2 = [Xwake(m,n)-Xwake(m1+1,n1),-Ywake(m,n)-Ywake(m1+1,n1),Zwake(m,n)-Zwake(m1+1,n1)];
                        aa = [r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                        bb = (r1(2)*r2(3)-r1(3)*r2(2))^2+(r1(3)*r2(1)-r1(1)*r2(3))^2+(r1(1)*r2(2)-r1(2)*r2(1))^2;
                        if bb > cutoff
                            cc = [r1/sqrt(r1(1)^2+r1(2)^2+r1(3)^2)-r2/sqrt(r2(1)^2+r2(2)^2+r2(3)^2)];
                            UVW4r = (aa/bb)*(ro(1)*cc(1)+ro(2)*cc(2)+ro(3)*cc(3))/4/pi;
                            UVW4r(2) = -UVW4r(2);
                        else
                            UVW4r = [0,0,0];
                        end
                        
                        UVW = UVW1+UVW2+UVW3+UVW4+UVW1r+UVW2r+UVW3r+UVW4r;
                        
                        %% Wake on wake influence matrix:
                        
                        C2_rollup(K,L,1) = UVW(1);
                        C2_rollup(K,L,2) = UVW(2);
                        C2_rollup(K,L,3) = UVW(3);
                        
                    end
                end
            end
        end
        
        %% perturb the wake (force-free wake):
        
        del_X = C1_rollup(:,:,1)*reshape(circ(:,:,i)',M*N,1)*t_step+C2_rollup(:,:,1)*reshape(circ_wake',(i-1)*N,1)*t_step;
        del_Y = C1_rollup(:,:,2)*reshape(circ(:,:,i)',M*N,1)*t_step+C2_rollup(:,:,2)*reshape(circ_wake',(i-1)*N,1)*t_step;
        del_Z = C1_rollup(:,:,3)*reshape(circ(:,:,i)',M*N,1)*t_step+C2_rollup(:,:,3)*reshape(circ_wake',(i-1)*N,1)*t_step;
        
        del_X = reshape(del_X,N+1,i-wake_start)';
        del_Y = reshape(del_Y,N+1,i-wake_start)';
        del_Z = reshape(del_Z,N+1,i-wake_start)';
        
        Xwake(wake_start:i-1,:) = Xwake(wake_start:i-1,:) + del_X;
        Ywake(wake_start:i-1,:) = Ywake(wake_start:i-1,:) + del_Y;
        Zwake(wake_start:i-1,:) = Zwake(wake_start:i-1,:) + del_Z;
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Storing instanteneous pressure and circulations on the
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% wing
%     count = 1;
%     for m = 1:M
%         for n = 1:N
%             WingX(:,count) = [count;Xr(m,n,i)-path_X(i);Xr(m,n+1,i)-path_X(i);Xr(m+1,n+1,i)-path_X(i);Xr(m+1,n,i)-path_X(i)];
%             WingY(:,count) = [count;Yr(m,n,i);Yr(m,n+1,i);Yr(m+1,n+1,i);Yr(m+1,n,i)];
%             WingZ(:,count) = [count;Zr(m,n,i);Zr(m,n+1,i);Zr(m+1,n+1,i);Zr(m+1,n,i)];
%             WingCirc(:,count) = [count;circ(m,n,i)];
%             WingPress(:,count) = [count;pres(m,n,i)/(.5*1.225*V_ref^2)];
%             WingAOA(:,count) = [count;alpha(m,n,i)];
%             count = count + 1;
%         end
%     end
%     
%     
%     dlmwrite(sprintf('WingX%d.dat',i),WingX','delimiter','\t','precision', '%.6f')
%     dlmwrite(sprintf('WingY%d.dat',i),WingY','delimiter','\t','precision', '%.6f')
%     dlmwrite(sprintf('WingZ%d.dat',i),WingZ','delimiter','\t','precision', '%.6f')
%     dlmwrite(sprintf('WingCirc%d.dat',i),WingCirc','delimiter','\t','precision', '%.6f')
%     dlmwrite(sprintf('WingPress%d.dat',i),WingPress','delimiter','\t','precision', '%.6f')
%     dlmwrite(sprintf('WingAOA%d.dat',i),WingAOA','delimiter','\t','precision', '%.6f')
%     
%     if (i>1)
%         count = 1;
%         for q = 1:i-1
%             for n = 1:N
%                 WakeX(:,count) = [count;Xwake(q,n)-path_X(i);Xwake(q,n+1)-path_X(i);Xwake(q+1,n+1)-path_X(i);Xwake(q+1,n)-path_X(i)];
%                 WakeY(:,count) = [count;Ywake(q,n);Ywake(q,n+1);Ywake(q+1,n+1);Ywake(q+1,n)];
%                 WakeZ(:,count) = [count;Zwake(q,n);Zwake(q,n+1);Zwake(q+1,n+1);Zwake(q+1,n)];
%                 WakeCirc(:,count) = [count;circ_wake(q,n)];
%                 count = count + 1;
%             end
%         end
%         
%     dlmwrite(sprintf('WakeX%d.dat',i),WakeX','delimiter','\t','precision', '%.6f')
%     dlmwrite(sprintf('WakeY%d.dat',i),WakeY','delimiter','\t','precision', '%.6f')
%     dlmwrite(sprintf('WakeZ%d.dat',i),WakeZ','delimiter','\t','precision', '%.6f')
%     dlmwrite(sprintf('WakeCirc%d.dat',i),WakeCirc','delimiter','\t','precision', '%.6f')
%     end
%     
%     
    
    %%%%%%%%%%%%%%%%%%%%%% end storing
     count = 1;
    for m = 1:M
        for n = 1:N
            fooX(:,count) = [Xr(m,n,i);Xr(m,n+1,i);Xr(m+1,n+1,i);Xr(m+1,n,i)];
            fooY(:,count) = [Yr(m,n,i);Yr(m,n+1,i);Yr(m+1,n+1,i);Yr(m+1,n,i)];
            fooZ(:,count) = [Zr(m,n,i);Zr(m,n+1,i);Zr(m+1,n+1,i);Zr(m+1,n,i)];
            count = count + 1;
        end
    end
    if (mod(i,2)==0)
    patch(fooX,fooY,fooZ,reshape(pres(:,:,i)'/(.5*1.225*V_ref^2),M*N,1)')
    hold on
    patch(fooX,-fooY,fooZ,reshape(pres(:,:,i)'/(.5*1.225*V_ref^2),M*N,1)')
    axis equal
    axis([-35 3 -4 4 -2 2])
    grid on
    view([-30 30])
    colorbar
    drawnow
    end
end



 %% Plotting:
    
    count = 1;
    for m = 1:M
        for n = 1:N
            fooX(:,count) = [Xr(m,n,i);Xr(m,n+1,i);Xr(m+1,n+1,i);Xr(m+1,n,i)];
            fooY(:,count) = [Yr(m,n,i);Yr(m,n+1,i);Yr(m+1,n+1,i);Yr(m+1,n,i)];
            fooZ(:,count) = [Zr(m,n,i);Zr(m,n+1,i);Zr(m+1,n+1,i);Zr(m+1,n,i)];
            count = count + 1;
        end
    end
    
    
        count = 1;
        for q = 1:i-1
            for n = 1:N
                fooX_wake(:,count) = [Xwake(q,n);Xwake(q,n+1);Xwake(q+1,n+1);Xwake(q+1,n)];
                fooY_wake(:,count) = [Ywake(q,n);Ywake(q,n+1);Ywake(q+1,n+1);Ywake(q+1,n)];
                fooZ_wake(:,count) = [Zwake(q,n);Zwake(q,n+1);Zwake(q+1,n+1);Zwake(q+1,n)];
                count = count + 1;
            end
        end
        
   
    
    
%     figure(1)
%     patch(fooX-path_X(i),fooY,fooZ,reshape(circ(:,:,i)',M*N,1)')
%     hold on
%     patch(fooX_wake-path_X(i),fooY_wake,fooZ_wake,reshape(circ_wake',(i-1)*N,1)')
%     hold on
%     patch(fooX-path_X(i),-fooY,fooZ,reshape(circ(:,:,i)',M*N,1)')
%     hold on
%     patch(fooX_wake-path_X(i),-fooY_wake,fooZ_wake,reshape(circ_wake',(i-1)*N,1)')
%         
%     axis equal
%     axis([0 30 -4 4 -2 2])
%     grid on
%     view([-30 30])
%     colorbar
%     
%     
%     count = 1;
%     for m = 1:M
%         for n = 1:N
%             WingX(:,count) = [count;Xr(m,n,i)-path_X(i);Xr(m,n+1,i)-path_X(i);Xr(m+1,n+1,i)-path_X(i);Xr(m+1,n,i)-path_X(i)];
%             WingY(:,count) = [count;Yr(m,n,i);Yr(m,n+1,i);Yr(m+1,n+1,i);Yr(m+1,n,i)];
%             WingZ(:,count) = [count;Zr(m,n,i);Zr(m,n+1,i);Zr(m+1,n+1,i);Zr(m+1,n,i)];
%             WingCirc(:,count) = [count;circ(m,n,i)];
%             WingPress(:,count) = [count;pres(m,n)/(.5*1.225*V_ref^2)];
%             WingAOA(:,count) = [count;alpha(m,n,i)];
%             count = count + 1;
%         end
%     end
%     
%     
%     dlmwrite(sprintf('WingX%d.dat',i),WingX','delimiter','\t','precision', '%.6f')
%     dlmwrite(sprintf('WingY%d.dat',i),WingY','delimiter','\t','precision', '%.6f')
%     dlmwrite(sprintf('WingZ%d.dat',i),WingZ','delimiter','\t','precision', '%.6f')
%     dlmwrite(sprintf('WingCirc%d.dat',i),WingCirc','delimiter','\t','precision', '%.6f')
%     dlmwrite(sprintf('WingPress%d.dat',i),WingPress','delimiter','\t','precision', '%.6f')
%     dlmwrite(sprintf('WingAOA%d.dat',i),WingAOA','delimiter','\t','precision', '%.6f')
%     
%     count = 1;
%         for q = 1:i-1
%             for n = 1:N
%                 WakeX(:,count) = [count;Xwake(q,n)-path_X(i);Xwake(q,n+1)-path_X(i);Xwake(q+1,n+1)-path_X(i);Xwake(q+1,n)-path_X(i)];
%                 WakeY(:,count) = [count;Ywake(q,n);Ywake(q,n+1);Ywake(q+1,n+1);Ywake(q+1,n)];
%                 WakeZ(:,count) = [count;Zwake(q,n);Zwake(q,n+1);Zwake(q+1,n+1);Zwake(q+1,n)];
%                 WakeCirc(:,count) = [count;circ_wake(q,n)];
%                 count = count + 1;
%             end
%         end
%         
%     dlmwrite('WakeX.dat',WakeX','delimiter','\t','precision', '%.6f')
%     dlmwrite('WakeY.dat',WakeY','delimiter','\t','precision', '%.6f')
%     dlmwrite('WakeZ.dat',WakeZ','delimiter','\t','precision', '%.6f')
%     dlmwrite('WakeCirc.dat',WakeCirc','delimiter','\t','precision', '%.6f')
%     
%     figure(2)
%     patch(fooX-path_X(i),fooY,fooZ,reshape(pres(:,:,i)'/(.5*1.225*V_ref^2),M*N,1)')
%     hold on
%     patch(fooX-path_X(i),-fooY,fooZ,reshape(pres(:,:,i)'/(.5*1.225*V_ref^2),M*N,1)')
%     
%     axis equal
%     axis([-1.5 1.5 -3.5 3.5 -1 1])
%     grid on
%     view([-50 50])
%     colorbar
%     
%     figure(3)
%     patch(fooX-path_X(i),fooY,fooZ,reshape(alpha(:,:,i)',M*N,1)')
%     hold on
%     patch(fooX-path_X(i),-fooY,fooZ,reshape(alpha(:,:,i)',M*N,1)')
%     
%     axis equal
%     axis([-1.5 1.5 -3.5 3.5 -1 1])
%     grid on
%     view([-50 50])
%     colorbar
%     figure(4)

% Thr = zeros(N,1);
% for j=1:i
%     Thr(:) = Thr(:) + sum(TSec(:,:,j))'
% end
% 
% Thr(:) = Thr(:)/i;
% hold on
% plot(Thr,'b')
    


% Compute time-averaged forces, power, effiency:

% Q = [trapz(t(N_Q_begin:N_steps+1),CL(N_Q_begin:N_steps+1)),trapz(t(N_Q_begin:N_steps+1),-CD(N_Q_begin:N_steps+1)),trapz(t(N_Q_begin:N_steps+1),power(N_Q_begin:N_steps+1))]/(t(N_steps+1)-t(N_Q_begin));
% 
% CL_ave = Q(1)
% CT_ave = Q(2)
% CP_ave = Q(3)/(.5*1.225*V_ref^3*b*c)
% eff = .5*Q(2)*1.225*V_ref^2*b*c*V_ref/Q(3)
% Q(4) = eff;

% figure(3)
% subplot(3,1,1)
% plot(t/(2*pi/omega)-.25,-CD)
% title('Thrust coefficient')
% subplot(3,1,2)
% plot(t/(2*pi/omega)-.25,CL)
% title('Lift coefficient')
% subplot(3,1,3)
% plot(t/(2*pi/omega)-.25,power/(.5*1.225*V_ref^3*b*c))
% title('Power coefficient')
