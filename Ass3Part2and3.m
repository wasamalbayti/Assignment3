%setting up variables and more
qCharge = 1.602e-19;
massE = 0.26*9.109e-31;
l = 200e-9; 
w = 100e-9;
Temperature = 300;
kb = 1.380649e-23;
avmeanTime = 0.2e-12;
numSim = 10000; 
Length = 2;
Width = 1;
Vo = 0.8;
bottleW = 0.4;
bottleL = 0.4;
nx = 100*Length;
ny = 100*Width;
sigma = zeros(ny,nx);

G = sparse(nx*ny,nx*ny);
F = zeros(nx*ny,1);


for x = 1:ny
    for y= 1:nx
        if(y >= nx*bottleW && y <= nx-nx*bottleW && (x >= ny-ny*bottleL || x <= ny*bottleL))
            sigma(x,y) = 1e-2;
        else
            sigma(x,y) = 1;
        end
    end
end
for y = 1:nx
    for x = 1:ny
        n = x + (y - 1)*ny;
        nxx = x + (y - 1-1)*ny;
        nxp = x + y*ny;
        nyy = x - 1+ (y - 1)*ny;
        nyp = x + 1+ (y - 1)*ny;
        
        if y == 1  
            G(n,n) = 1;
        elseif y == nx
            G(n,n) = 1;
        elseif (x == 1) 
            upperY = (sigma(x,y)+sigma(x+1,y))/2;
            upperX = (sigma(x,y)+sigma(x,y+1))/2;
            lowerX = (sigma(x,y)+sigma(x,y-1))/2;
            G(n,nyp) = upperY;
            G(n,nxp) = upperX;
            G(n,nxx) = lowerX;  
            G(n,n) = -(upperY + upperX + lowerX);
   
        elseif (x == ny) 
            lowerY = (sigma(x,y)+sigma(x-1,y))/2;
            upperX = (sigma(x,y)+sigma(x,y+1))/2;
            lowerX = (sigma(x,y)+sigma(x,y-1))/2;
            G(n,nyy) = lowerY;
            G(n,nxp) = upperX;
            G(n,nxx) = lowerX;
            G(n,n) = -(lowerY + upperX + lowerX);

        else 
            upperY = (sigma(x,y)+sigma(x+1,y))/2;
            lowerY = (sigma(x,y)+sigma(x-1,y))/2;
            upperX = (sigma(x,y)+sigma(x,y+1))/2;
            lowerX = (sigma(x,y)+sigma(x,y-1))/2;
            
            G(n,nyp) = upperY;
            G(n,nyy) = lowerY;
            G(n,nxp) = upperX;
            G(n,nxx) = lowerX;
            G(n,n) = -(upperY + lowerY + upperX + lowerX);

        end
    end
end

for y = 1:nx
    for x = 1:ny
        n = x + (y - 1)*ny;        
        if y == 1 
            F(n) = Vo;
        end
    end
end

V = G\F;

for y = 1:nx
    for x = 1:ny
        n = x + (y - 1)*ny;
        Emap(x,y) = V(n);
    end
end

figure(1)
surf(Emap)
xlabel('Length')
ylabel('Width')
zlabel('Voltage')
title('Surface plot of V(x,y)')
pause (1)
[Ex,Ey] = gradient(-Emap,1e-9);
figure(2)
quiver(Ex,Ey)
xlabel('Length')
ylabel('Width')
title('Electric field vector plot')
ylim([0,100])
xlim([0,200])
pause (1)
VelocityTherm = sqrt((2*kb*Temperature)/massE);
%randomize x and y position within region
xPosition = l.*rand(numSim,2);
yPosition = w.*rand(numSim,2);

for y = 1:numSim
    while ((xPosition(y,2) > 0.8e-7 && xPosition(y,2) < 1.2e-7) && ...
            (yPosition(y,2) > 0.6e-7 ||  yPosition(y,2) < 0.4e-7))
        xPosition(y,2) = l.*rand(1,1);
        yPosition(y,2) = w.*rand(1,1);
    end

end
tStep = 0.22e-9/VelocityTherm/2;

% assign random velocity to electrons
xVelocity = randn(numSim,2)*sqrt((kb*Temperature)/massE);
yVelocity = randn(numSim,2)*sqrt((kb*Temperature)/massE);
xPosition(:,1) = xPosition(:,2);
yPosition(:,1) = yPosition(:,2);
sum = 0;
% the new position with velocity and time step
newXPos = tStep*xVelocity(:,1);
newYPos = tStep*yVelocity(:,1);
% scatter
PScat = 1 - exp(-(tStep/avmeanTime));
%position before colliding and after
XoldPosition = xPosition(:,1);
YoldPosition = yPosition(:,1);

for t = 1:1000
    for n = 1:numSim 
        if(round(1e9*xPosition(n,1)) <= 0 || round(yPosition(n,1)*1e9) <= 0)
            Forcex = Ex(ceil(yPosition(n,1)*1e9),ceil(1e9*xPosition(n,1))).*qCharge; 
            Forcey = Ey(ceil(yPosition(n,1)*1e9),ceil(1e9*xPosition(n,1))).*qCharge; 
        elseif(round(yPosition(n,1)*1e9) > 100 || round(1e9*xPosition(n,1)) > 200)
            Forcex = Ex(floor(yPosition(n,1)*1e9),floor(1e9*xPosition(n,1))).*qCharge; 
            Forcey = Ey(floor(yPosition(n,1)*1e9),floor(1e9*xPosition(n,1))).*qCharge; 
        else
            Forcex = Ex(round(yPosition(n,1)*1e9),round(1e9*xPosition(n,1))).*qCharge; 
            Forcey = Ey(round(yPosition(n,1)*1e9),round(1e9*xPosition(n,1))).*qCharge; 
        end 
        % calcualte the acclerations now with the forces
        accelerationx = Forcex/massE;
        accelerationy = Forcey/massE;
        
        % get the velocities 
        xVelocity(n,1) = xVelocity(n,1) + accelerationx*tStep;
        yVelocity(n,1) = yVelocity(n,1) + accelerationy*tStep;
        
        if (rand < PScat)
            xVelocity(n,:) = randn(1)*sqrt((kb*Temperature)/massE);
            yVelocity(n,:) = randn(1)*sqrt((kb*Temperature)/massE);
            newXPos(n,:) = tStep*xVelocity(n,1);
            newYPos(n,:) = tStep*yVelocity(n,1);                       
            XoldPosition = xPosition(:,1);
            YoldPosition = yPosition(:,1);
        else
            newXPos(n,:) = tStep*xVelocity(n,1);
            newYPos(n,:) = tStep*yVelocity(n,1);
        end
        if (xPosition(n,1) + newXPos(n) > 2e-7)
            xPosition(n,2) = newXPos(n) + xPosition(n,1) - 2e-7;
        elseif (xPosition(n,1) + newXPos(n) < 0)
            xPosition(n,2) = newXPos(n) + xPosition(n,1) + 2e-7;
        else
            xPosition(n,2) = newXPos(n) + xPosition(n,1);
        end
       
        if ((yPosition(n,1) + newYPos(n) > 1e-7)|| (yPosition(n,1) + newYPos(n) < 0))
            newYPos(n) = -newYPos(n);
            yPosition(n,2) = newYPos(n) + yPosition(n,1);
        elseif(xPosition(n,1) > 1.2e-7 && xPosition(n,1) < 2e-7 &&  yPosition(n,1) > 0.6e-7 && yPosition(n,1) < 1e-7 &&(xPosition(n,1) + newXPos(n) < 1.2e-7))
            newXPos(n) = -newXPos(n);
            xPosition(n,2) = newXPos(n) + xPosition(n,1);
        elseif (xPosition(n,1) > 0 && xPosition(n,1) < 0.8e-7 &&  yPosition(n,1) > 0.6e-7 && yPosition(n,1) < 1e-7 && (xPosition(n,1) + newXPos(n) > 0.8e-7))
            newXPos(n) = -newXPos(n);
            xPosition(n,2) = newXPos(n) + xPosition(n,1);
        elseif (xPosition(n,1) > 0 && xPosition(n,1) < 0.8e-7 && (xPosition(n,1) + newXPos(n) > 0.8e-7) && yPosition(n,1)> 0 && yPosition(n,1) < 0.4e-7)
            newXPos(n) = -newXPos(n);
            xPosition(n,2) = newXPos(n) + xPosition(n,1);
        elseif(xPosition(n,1) > 1.2e-7 && xPosition(n,1) < 2e-7 &&  yPosition(n,1) > 0 && yPosition(n,1) < 0.4e-7 &&(xPosition(n,1) + newXPos(n) < 1.2e-7))
            newXPos(n) = -newXPos(n);
            xPosition(n,2) = newXPos(n) + xPosition(n,1);
        elseif(xPosition(n,1) > 0.8e-7 && xPosition(n,1) < 1.2e-7 && yPosition(n,1) > 0.4e-7 && yPosition(n,1) < 0.6e-7 && (yPosition(n,1) + newYPos(n) < 0.4e-7))
            newYPos(n) = -newYPos(n);
            yPosition(n,2) = newYPos(n) + yPosition(n,1); 
        elseif(xPosition(n,1) > 0.8e-7 && xPosition(n,1) < 1.2e-7 && yPosition(n,1) > 0.4e-7 && yPosition(n,1) < 0.6e-7 &&(yPosition(n,1) + newYPos(n) > 0.6e-7))
            newYPos(n) = -newYPos(n);
            yPosition(n,2) = newYPos(n) + yPosition(n,1);            
        else
            yPosition(n,2) = newYPos(n) + yPosition(n,1);
        end
    end    
    timeVector(t) = tStep*t;
    
    if (sum == 0)
        figure(3)
        scatter(xPosition([1:10],2),yPosition([1:10],2),1)
        hold on
        plot([0.8e-7,0.8e-7],[0,0.4e-7],'k',[0.8e-7,1.2e-7],[0.4e-7,0.4e-7],'k',[1.2e-7,1.2e-7],[0.4e-7,0],'k',[0.8e-7,1.2e-7],[0,0],'k')
        plot([0.8e-7,0.8e-7],[1e-7,0.6e-7],'k',[0.8e-7,1.2e-7],[0.6e-7,0.6e-7],'k',[1.2e-7,1.2e-7],[0.6e-7,1e-7],'k',[0.8e-7,1.2e-7],[1e-7,1e-7],'k')
       
        title(['particle trajector'])
        xlabel('X')
        ylabel('Y')
        xlim([0 200e-9])
        ylim([0 100e-9])
    elseif (sum > 0 && sum < 1000)
        title(['particle trajectory'])
        scatter(xPosition([1:10],2),yPosition([1:10],2),1)
        hold on
    else
        title(['particle trajector'])
        scatter(xPosition([1:10],2),yPosition([1:10],2),1)
        hold off
    end 
    pause(0.0001)
    xPosition(:,1) = xPosition(:,2);
    yPosition(:,1) = yPosition(:,2);
    sum = sum + 1;
end

%%%%%Question 3
figure(5)
hist3([xPosition(:,1),yPosition(:,1)],[10,20])
title('electron density map')
xlabel('X')
ylabel('Y')
zlabel('Num of particles')

for numberSolution = 1:5
    bottleW = 0.4;
    bottleL = 0.4*(1+(numberSolution/20));
    nx = 100*Length;
    ny = 100*Width;
    sigma = zeros(ny,nx);
    for x = 1:ny
        for y = 1:nx
            if(y >= nx*bottleW && y <= nx-nx*bottleW &&  (x >= ny-ny*bottleL || x <= ny*bottleL))
                sigma(x,y) = 1e-2;
            else
                sigma(x,y) = 1;
            end
        end
    end
    G = sparse(nx*ny,nx*ny);
    F = zeros(nx*ny,1);
    for y = 1:nx
        for x = 1:ny
            n = x + (y - 1)*ny;
            nxx = x + (y - 1-1)*ny;
            nxp = x + y*ny;
            nyy = x - 1 + (y - 1)*ny;
            nyp = x + 1 + (y - 1)*ny;

            if y == 1 
                G(n,n) = 1;
            elseif y == nx
                G(n,n) = 1;
            elseif (x == 1) 
                upperY = (sigma(x,y)+sigma(x+1,y))/2;
                upperX = (sigma(x,y)+sigma(x,y+1))/2;
                lowerX = (sigma(x,y)+sigma(x,y-1))/2;

                G(n,n) = -(upperY + upperX + lowerX);
                G(n,nyp) = upperY;
                G(n,nxp) = upperX;
                G(n,nxx) = lowerX;

            elseif (x == ny) 
                lowerY = (sigma(x,y)+sigma(x-1,y))/2;
                upperX = (sigma(x,y)+sigma(x,y+1))/2;
                lowerX = (sigma(x,y)+sigma(x,y-1))/2;

                G(n,n) = -(lowerY + upperX + lowerX);
                G(n,nyy) = lowerY;
                G(n,nxp) = upperX;
                G(n,nxx) = lowerX;
            else 
                upperY = (sigma(x,y)+sigma(x+1,y))/2;
                lowerY = (sigma(x,y)+sigma(x-1,y))/2;
                upperX = (sigma(x,y)+sigma(x,y+1))/2;
                lowerX = (sigma(x,y)+sigma(x,y-1))/2;

                G(n,n) = -(upperY + lowerY + upperX + lowerX);
                G(n,nyp) = upperY;
                G(n,nyy) = lowerY;
                G(n,nxp) = upperX;
                G(n,nxx) = lowerX;
            end
        end
    end
    for y = 1:nx
        for x = 1:ny
            n = x + (y - 1)*ny;
            if (y == 1) 
                F(n) = Vo;
            end
        end
    end
    V = G\F;
    for y = 1:nx
        for x = 1:ny
            n = x + (y - 1)*ny;
            Emap(x,y) = V(n);
        end
    end
    [Ex,Ey] = gradient(-Emap,1e-9);
    currentY(numberSolution) = mean(sigma(:,1).*Ex(:,1));
    currentX(numberSolution) = mean(sigma(:,1).*Ex(:,1));
    currentDensity= sqrt(currentX.^2 + currentY.^2);
end
figure(6)
plot(linspace(1,5,5),current)
xlabel('bottleneck width')
ylabel('current Density')
title('current Density vs bottleneck')
grid on
