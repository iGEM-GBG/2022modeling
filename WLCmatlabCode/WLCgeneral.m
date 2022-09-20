clc, clear, close all

Compound = "DNA";          % "prot" or "DNA" to model protein or DNA respectively
disp("The chosen compound is " + Compound)

N = 50;                 % Number of segments (less effect on smoothness than aT) %Doesn't seem to affect the result very much           
aT = 10^5;              % Number of chains    N=20  & 10^5 : 17 sek  /  10^6:169 sek

if Compound == "prot"
    L = 10;                 % Protein chain length (nm)
    lp=0.6;                 % Protein persistence length (nm)
elseif Compound == "DNA"
    bp = 30;              % Number of base pairs
    L = 0.34*bp;            % DNA chain length (1 base pair: 0.34 nm, ~30 base pairs)
    lp=39;                  % DNA persistence length (nm). Source: https://www.nature.com/articles/nphys2002
else
    L = 10;                % Chain length (nm)
    lp=0.45;               % Persistence length (0.45 nm for some protein linkers, 40 nm for DNA)
end

sigma = sqrt(L/(N*lp)); % standard deviation

R = zeros(1,aT);        % Radius vector
Rx=zeros(aT,N);          % x-direction
Ry=zeros(aT,N);          % y-direction
Rz=zeros(aT,N);          % z-direction

for j=1:aT    

    % Create initial segment with given angles
    if j <= aT/2    % Starting direction for first linker
        theta_0 = pi/2-0.001;
        phi_0 = 0+0.001;
    else            % Starting direction for second linker
        theta_0 = pi/2-0.001;
        phi_0 = pi-0.001; % pi
    end

    % Giving position of first vector in linker
    Rx(j, 1) = L/N *cos(phi_0).*sin(theta_0);
    Ry(j, 1) = L/N *sin(phi_0).*sin(theta_0);
    Rz(j, 1) = L/N *cos(theta_0);
    
    rho = [normrnd(0,sigma,[1,N])];     %Angle between consecutive segments in linker. Diffrent angles have diffrent values
    theta=2*pi*(rand(1,N)-1/2);         %Rotation between consecutive segments in lnker. All angles have equal probability

    % Making a linker from generated angles
    u = [Rx(j,1); Ry(j,1); Rz(j,1)]; % Position vector of one (initialy the first) segment in linker
    for i=1:N-1
        % Rotation of vectors is done with Rodrigues' formula, using the matrix notation u_rot = R*u.
        % See https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula

        uNorm=u/norm(u); %Normalized u
        
        % Rotation by rho
        n = [1,-u(1)/u(2),0]/norm([1,-u(1)/u(2),0]); % Make an arbitrary normalized orthogonal vector to u           
        W=[ 0,    -n(3),   n(2)
            n(3),  0,     -n(1)
           -n(2),  n(1),   0   ]; % Rotation matrice for angle rho between new section and previous section
        RotMat_rho=eye(3) + sin(rho(i))*W + (1-cos(rho(i)))*W^2;    % Rodrigue's rotation matrix
        u=RotMat_rho*u;    % u after rotation by rho


        % Rotation by theta
        W=[ 0,    -uNorm(3),   uNorm(2)
            uNorm(3),  0,     -uNorm(1)
           -uNorm(2),  uNorm(1),   0  ]; % Rotation matrice for angle phi between new section and previous section
        RotMat_theta=eye(3)+sin(theta(i))*W+(1-cos(theta(i))).*W^2;  % Rodrigue's rotation matrix
        u=RotMat_theta*u; % u after rotation by theta

        %Store new vector
        Rx(j,i+1)= u(1);
        Ry(j,i+1)= u(2);
        Rz(j,i+1)= u(3);
    end    
    
    R(j) = norm([sum(Rx(j,:)), sum(Ry(j,:)), sum(Rz(j,:))]);
end

% Plot the chain generated last
figure(2)
% cumsum is used for adding elements cumulatively for the coordinates of
% all segment chain joints
Rx_cumsum = cumsum(Rx(end,:));
Ry_cumsum = cumsum(Ry(end,:));
Rz_cumsum = cumsum(Rz(end,:));

plot3(Rx_cumsum,Ry_cumsum,Rz_cumsum,'k') %Plot the chain!
axis image
xlabel("X")
ylabel("Y")
zlabel("Z")


% Create histogram
ti = "Estimated probability density";
figure(3)
lb = 0; % Lower boundary for the histogram. >0 to avoid noise
ub = L; % Upper boundary for the histogram
step = (ub-lb)/(2*sqrt(aT)); % Step size for the points in the histogram
[y,edges] = histcounts(R,lb:step:ub); % Creating points for histogram
x=lb:step:(ub-step); % Calculate the x-value for each data point in the histogram
k=1/trapz(x,y);  % trapz: Trapezoidal numerical integration
plot(x,y*k,'b') % Plot the data points in the histogram
xlabel("Chain reach [$nm$]", 'interpreter','latex','FontSize',14)
ylabel("Density function", 'interpreter','latex','FontSize',14)
title(ti)
%eval(sprintf("title('%s %s','interpreter','latex','FontSize',14)",ti,num2str(N)))
