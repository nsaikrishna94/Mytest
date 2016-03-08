c = 1;                                                   % Chord length 
t = 0.15*c;                                              % Maximum Thickness(in terms of chord length)
n = 100;                                                 % No. of points
first = 0.2*c;                                           % left limit of porous region
last = 0.8*c;                                            % right limit of porous region
phi = 0.8;                                               % Porosity 
depth = 0.05*t;                                          % Depth of the porous region
noc = 2;                                                 % No.of circles in the porous region
                                    
div = 5000;                                              % div is used for no. of interpolants
mult = 10;                                               % mult*(n2-n1) is for the no. of circles in a row

%%%%%%%%%%%%%%%%%%%%%%%%%%% Airfoil Building %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% x co-ordinates of points on airfoil
xu = 0:c/n:c;                                            % Array of x co-ordinates of points of top surface
xd = 0:c/n:c;                                            % Array of x co-ordinates of points of bottom surface

% y co-ordinates of points on airfoil
yu =  (t/.2)*(0.2969*sqrt(xu) - 0.126*xu + (-0.3516 + (0.2843 - 0.1015*xu).*xu).*xu.*xu);  % Array of y co-ordinates of points of top surface
yd =  -(t/.2)*(0.2969*sqrt(xd)- 0.126*xd + (-0.3516 + (0.2843 - 0.1015*xd).*xd).*xd.*xd);  % Array of y co-ordinates of points of bottom surface

% Cubic Spine Interpolation
xx = 0:c/div:c;
yy = spline(xu,yu,xx);

[w1 w2] = size(xx);                                      % Size of xx

% slopes and normals of top surface of airfoil
slope = zeros(2,div);
norm  = zeros(2,div);

for i = 1:div
    
    slope(1,i) = (xx(1,i+1)-xx(1,i))/sqrt((yy(1,i+1)-yy(1,i))^2 + (xx(1,i+1)-xx(1,i))^2);
    slope(2,i) = (yy(1,i+1)-yy(1,i))/sqrt((yy(1,i+1)-yy(1,i))^2 + (xx(1,i+1)-xx(1,i))^2);
    
    norm(1,i)  = -(yy(1,i+1)-yy(1,i))/sqrt((yy(1,i+1)-yy(1,i))^2 + (xx(1,i+1)-xx(1,i))^2); 
    norm(2,i)  = (xx(1,i+1)-xx(1,i))/sqrt((yy(1,i+1)-yy(1,i))^2 + (xx(1,i+1)-xx(1,i))^2);  

end;

mx = zeros(1,div);
my = zeros(2,div);

for i = 1:div
    
    mx(1,i) = 0.5*(xx(1,i) + xx(1,i+1));   % Co-ordinates of mid points of Interpolants of surface of airfoil
    my(1,i) = 0.5*(yy(1,i) + yy(1,i+1));

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Circles Formation %%%%%%%%%%%%%%%%%%%%%%%%%%%%

i = 1;
while first > xx(1,i)
    i = i+1;                               % This block of code is used for finding n1(left limit)
end;
n1 = i;

while last > xx(1,i)
    i = i+1;                               % This block of code is used for finding n1(right limit)
end;
n2 = i;

% The parameters defined in this block are useful for finding the co-ordinates of the centers and radii of first column circles
alpha = slope(2,n1)/slope(1,n1);
disc = -alpha*mx(1,n1) + my(1,n1) + first*sqrt(1+alpha^2);
a = disc/(sqrt(1+alpha^2) - alpha);
b = disc;
            
ratio = 2*sqrt((1-phi)/pi);
S = depth/(noc-1+ratio);
D = ratio*S;

cen = zeros(2,mult*(n2-n1),noc);           % Declaring matrix to store co-ordinates of the centers of circles

cen(1,1,1) = first + 0.5*D;
cen(2,1,1) = -(b/a)*cen(1,1,1) + disc;

prevx = cen(1,1,1) + 0.5*D*norm(1,n1);
prevy = cen(2,1,1) + 0.5*D*norm(2,n1);

if(noc ~= 1)
    for i = 2:noc
        cen(1,1,i) = cen(1,1,1) - (i-1)*S*norm(1,n1);
        cen(2,1,i) = cen(2,1,1) - (i-1)*S*norm(2,n1);
    end;
end; 

k = 1;
j = n1;

while xx(1,j) < last
    
    k = k+1;
    dS1 = sqrt((prevx-xx(1,j+1))^2 + (prevy-yy(1,j+1))^2);
                
    if(S-dS1 > 0)                          % In this block j and next are updated based on this condition (S-dS1 > 0)
        j = j+1;
        next = S - dS1;
        dS2 = sqrt((xx(1,j+1)-xx(1,j))^2 + (yy(1,j+1)-yy(1,j))^2);
        
        while next-dS2 > 0
            j = j+1;
            next = next - dS2;
            dS2 = sqrt((xx(1,j+1)-xx(1,j))^2 + (yy(1,j+1)-yy(1,j))^2);
        end;
        
        prevx = xx(1,j) + next*slope(1,j); % Updating the prevx,prevy 
        prevy = yy(1,j) + next*slope(2,j);
    else    
        next = S;
        prevx = prevx + next*slope(1,j);   % Updating the prevx,prevy 
        prevy = prevy + next*slope(2,j);
    end;
                
    cen(1,k,1) = prevx - 0.5*norm(1,j)*D;
    cen(2,k,1) = prevy - 0.5*norm(2,j)*D;
                
    if(noc ~= 1) 
        for i = 2:noc
            cen(1,k,i) = cen(1,k,1) - (i-1)*S*norm(1,j);
            cen(2,k,i) = cen(2,k,1) - (i-1)*S*norm(2,j);
        end;
    end; 

end;


display(k);

np = 51;                                   % No. of points on a circle
cir = zeros(2,np,mult*(n2-n1),noc);

for i = 1:mult*(n2-n1)
    
    for j = 1:noc
                        
        for rot = 1:np
            cir(1,rot,i,j) = cen(1,i,j) + 0.5*D*cos(2*pi*(rot-1)/(np-1));
            cir(2,rot,i,j) = cen(2,i,j) + 0.5*D*sin(2*pi*(rot-1)/(np-1));
        end;
                        
    end;
     
end;
            
       
%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%

%plot(xx(1,:),yy(1,:),'-');           % Plot of Airfoil
plot(xx(1,1:n1),yy(1,1:n1),'-');
hold on
plot(xx(1,n2:div),yy(1,n2:div),'-');
hold on 
plot(xd(1,:),yd(1,:),'-');
hold on 

for j = 1:noc           
    
    for col = 1:k
        
        plot(cir(1,:,col,j),cir(2,:,col,j));   % Plot of a circle (which lies completely inside the airfoil)
        hold on
                
    end;
            
end;
        
daspect([1 1 1]);
title('Airfoil with a porous region');      % Title labelling
xlabel('Along the chord (x)');              % x label
ylabel('Thickness (y)');                    % y label
        