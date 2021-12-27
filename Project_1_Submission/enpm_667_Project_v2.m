cf = 117000;
cr = 108000;
lf = 1.4;
lr = 1.65;
Iz = 3234;
m = 1650; %mass of the vehicle
rw = 0.35179;
vxd = 25; 
vx = vxd; %vehicle velocity

delta_t = 0.4; % preview time
W = 1.861;
dT = 0.01; % Time step 
Np = 2.0/dT; % Prediction Horizon
Nc = 0.2*Np;
nu = 0.7; % Friction Coeff
m_w = 24.5; % mass of wheel
I_w = m_w*rw*rw/2; % Inertia of wheel
g = 9.81;

C = [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0 ; 0 0 0 1 0];
Cd = C;
X = zeros(5*Np, 1);
Y = zeros(4*Np, 1);
Y_des = zeros(4*Np, 1);
U = zeros(5*Nc, 1);



% disp(size(Dp))
Q = 0.01*eye([4*Np, 4*Np]);
R = 0.0001*eye([5*Nc, 5*Nc]);

Ydes = [0; 0; 0; 0];
rho = 0.0001;
epsilon = 0.01;
epsilon_max = 0.1;
k_curv = 0;
t = 0;
x = zeros(5, 10/dT);
Gk = zeros(5*Nc+1, 5*Nc+1);
Hk = zeros(1, 5*Nc+1);
Wex0 = 0.01;
Wf0 = 0.0001;
Wey = 0.01;
Wdpsi = 0.01;
Wvy = 0.01;
as = 0.001;
bs = -0.0008;
cs = 20;
ks = 20;
kw = 20;
Kr = zeros(4,1);
Wf = zeros(4,1);
ey_the = 0.5;
dpsi_the = 5;
beta_the = 15*pi/180;
umin = [-0.7; 0; 0; 0; 0];
umax = [0.7; 0.25*(nu*m*g + m*1); 0.25*(nu*m*g + m*1); 0.25*(nu*m*g + m*1); 0.25*(nu*m*g + m*1)];
dumax = [0.03; 0.43*0.25*(nu*m*g + m*1); 0.43*0.25*(nu*m*g + m*1); 0.43*0.25*(nu*m*g + m*1); 0.43*0.25*(nu*m*g + m*1)];
ymin = [-5; -5; -0.1; -5];
ymax = [5; 5; 0.1; 5];

up = [0 ; nu*m*g/4; nu*m*g/4; nu*m*g/4; nu*m*g/4];
Umax = umax;
Umin = umin;
Ymin = ymin;
Ymax = ymax;
y = zeros(4, 10/dT);
y(:,1) = Cd*x(:,1);
for i = 2:Nc
    Umax = cat(1,Umax, umax);
    Umin = cat(1,Umin, umin);
end
Umax = cat(1,Umax, epsilon_max);
Umin = cat(1,Umin, 0);

for i=2:Np
    Ydes = cat(1, Ydes, [0; 0; 0; 0]);
    Ymax = cat(1,Ymax,ymax);
    Ymin = cat(1, Ymin, ymin);
end
for k = 1:1:((10/dT) - 1)
   
    
    if t<2
       k_curv = 0;
    elseif t<3
        k_curv = 0.01*(t-2);
    elseif t<5
        k_curv = 0.01;
    elseif t<6
        k_curv = -0.01*t + 0.06;
    else
        k_curv = 0;
    end
    t = t + dT;
    ex = y(1,k);
    eyp = y(2,k);
    dpsi = y(3,k);
    vy = y(4,k);
    vx = ex + vxd;
    beta = atan2(vy, vx);
    DL = vx * delta_t; % preview distance
    w = [(4*vx)/(rw); (4*vx)/(rw); (4*vx)/(rw); (4*vx)/(rw)]; % angular wheel speed
    A = [0 0 0 0 0 ; 0 0 vx 1 DL ; 0 0 0 0 1; 0 0 0 -((cf+cr)/(m*vx)) -vx - (cf*lf - cr*lr)/(m*vx); 0 0 0 (cr*lr -cf*lf)/(Iz*vx) -(cf*lf*lf + cr*lr*lr)/(Iz*vx)];

    B = [0 1/m 1/m 1/m 1/m ; 0 0 vx 1 DL; 0 0 0 0 1; cf/m 0 0 0 0; (cf*lf)/Iz -W/(2*Iz) W/(2*Iz) -W/(2*Iz) W/(2*Iz)];
    Ad = eye([5,5]) + A*dT;
    Bd = B*dT;
    Cp = Cd*Ad;
    for i=2:1:Np
        Cp = cat(1, Cp, Cd*(Ad^i));
    end
    %disp(size(Cp));

    Ep = cat(2, Cd, zeros(4,5*(Np - 1)));

    for i = 2:Np
        Ep_row = Cd*Ad^(i-1);
        for j = 2:Np
            if j >i
                Ep_row = cat(2, Ep_row, zeros(4,5)); 
        
            elseif j == i
                Ep_row = cat(2, Ep_row, Cd);
        
            else
                Ep_row = cat(2, Ep_row, Cd*Ad^(i-j));
            end
        end
        Ep = cat(1, Ep, Ep_row);
    end

    %disp(Ep(9:12, 1:15));

    Dp = cat(2, Cd*Bd, zeros(4, 5*(Nc - 1)));   
    for i = 2:Np
        Dp_row = Cd*(Ad^(i-1));
        for j = 2:Nc
            if i <= Nc
                if j >i
                    Dp_row = cat(2, Dp_row, zeros(4,5));
                elseif j == i
                    Dp_row = cat(2, Dp_row, Cd*Bd);
                else
                    Dp_row = cat(2, Dp_row, Cd*(Ad^(i-j)));
                end
            else
                Dp_row = cat(2, Dp_row, Cd*(Ad^(i-j)));
            end
        end
        Dp = cat(1, Dp, Dp_row);
    end
    Acons = cat(1,Dp, -Dp);
    Acons = cat(2, Acons, -1*ones(8*Np, 1));
    
    
    Qs = max(abs(eyp)/ey_the, max(abs(dpsi)/dpsi_the, abs(beta)/beta_the));
    if k ~= 1
        w(1) = w(1) + up(2)*rw*dT/(I_w*cos(up(1)*pi/180)); % Angular Wheel Speeds
        w(2) = w(2) + up(3)*rw*dT/(I_w*cos(up(1)*pi/180));
        w(3) = w(3) + up(4)*rw*dT/I_w;
        w(4) = w(4) + up(5)*rw*dT/I_w;
    end
    
    for z = 1:4
        if vx > rw*w(z)
            Kr(z) = abs(rw*w(z)/vx - 1);
        else
            Kr(z) = abs(1 - vx/(rw*w(z)));
        end
        if(Kr(z) <= 0.1)
            Wf(z) = Wf0;
        else
            Wf(z) = Wf0*exp(kw*(Kr(z) -0.1));
        end
    end
    
    for j = 1:Np
        if Qs <=1
            Q(4*j-3,4*j-3) = Wex0;
        else
            Q(4*j-3,4*j-3) = as*bs*tanh((ks/Qs)^cs);
        end
    end
    
    for o = 1:Nc
        R(5*o - 3, 5*o - 3) = Wf(1);
        R(5*o - 2, 5*o - 2) = Wf(2);
        R(5*o - 1, 5*o - 1) = Wf(3);
        R(5*o, 5*o) = Wf(4);
    end
    Gk(1:(5*Nc), 1:(5*Nc)) = 2*(transpose(Dp)*Q*Dp + R);
    Gk(5*Nc+1, 5*Nc+1) = 2*rho;
    Hk(1:5*Nc) = 2*transpose(Cp*x(:,k)-Ydes)*Q*Dp;
    Umin(1:5) = max(umin, up-dumax);
    Umax(1:5) = min(umax, up+dumax);
    Bcons = Ymax - Cp*x(:,k);
    Bcons = cat(1,Bcons, -Ymin+Cp*x(:,k)); 
    disp(size(Gk));
    U = quadprog(Gk, transpose(Hk), Acons, Bcons, [], [], [], Umax);
    disp(size(U));
    up = U(1:5);
    
    w_dist = [(vxd - vx)/(Np*dT); vx*k_curv; 0; 0; 0];
    
    x(:,k+1) = Ad*x(:,k) + Bd*U(1:5) + dT*w_dist;
    y(:, k+1) = Cd*x(:, k+1);
end