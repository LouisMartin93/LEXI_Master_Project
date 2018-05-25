%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   Fit a 2D Gaussian Function to Data
%
%  Created by G. Nootz, May 2012. and modif by: Manuel A. Diaz, Jan 2016.
%  and modif by: Louis Martin, May 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE:  
%   Fit a 2D gaussian centroid to data.
% 
% INPUT:
%   Data: two-dimensional array of size nxn.
%   x0 = [Amp,x0,wx,y0,wy,theta]: Inital guess parameters.
%   x = [Amp,x0,wx,y0,wy,theta]: simulated gaussian parameters.
%
%
% OUTPUT: 
%   Gaussian function parameters.
%
% NOTE:
%   1.This routine uses Matlab's 'lsqcurvefit' function to fit.
%   2.The initial values in x0 must be close to x in order for the fit
%   to converge to the values of x (especially if noise is added).
% 
% MODIFS
%   -Code is re-formulated into a single function for simpler implementation.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function centre = Centroid2DGaussian(image,varargin)
    %% ---Fitting Functions---
    %
    % Coeficients A convention:
    %	A = [Amplitude, x0, x-Width, y0, y-Width, Angle(in Radians)]
    %
    % X-data convention:
    %	X is of size(n,n,2) where 
    %	X(:,:,1) : x-coordinates,  
    %	X(:,:,2) : y-coordinates.
    %
    % In this numerical test we use two-dimensional fitting functions:
    % 1. 2D Rotated Gaussian function ( A requires 6 coefs ).
    % essentially an ellipse fit.
    f = @(A,X) A(1)*exp( -(...
        ( X(:,:,1)*cos(A(6))-X(:,:,2)*sin(A(6)) - A(2)*cos(A(6))+A(4)*sin(A(6)) ).^2/(2*A(3)^2) + ... 
        ( X(:,:,1)*sin(A(6))+X(:,:,2)*cos(A(6)) - A(2)*sin(A(6))-A(4)*cos(A(6)) ).^2/(2*A(5)^2) ) );
    %% Checking and seting Varargin
    numvarargs = length(varargin);
    if numvarargs > 1
        error('Centroid2DGaussian.m: Too many outputs requires only 1 aditional output')
    end
    
    optargs = {false};
    optargs(1:numvarargs) = varargin;
    PlotPlots = optargs{:};

    %%  Crop ROI
    % crops the main region around the central PSF so a gaussian can be fit
    % easier to that. 

    image = image - min(min(image));
    %figure(1);
    %imagesc(log10(image));
    %colorbar();
    [val,xindex]  = max(max(image));
    [val,yindex] = max(image(:,xindex));

    dims = 100;
    %figure(2);
    %imagesc(log10(image(yindex-dims/2:yindex+dims/2,xindex-dims/2:xindex+dims/2)));
    %colorbar();
    CropImage = image(yindex-dims/2:yindex+dims/2,xindex-dims/2:xindex+dims/2);


    %% ---Build numerical Grids---

    [n,m]=size(CropImage);
    % Numerical Grid
    [x,y]=meshgrid(-n/2:n/2-1,-m/2:m/2-1);
    X=zeros(m,n,2);
    X(:,:,1)=x; X(:,:,2)=y;
    % High Resolution Grid
    h=3; [xh,yh]=meshgrid(-n/2:1/h:n/2,-m/2:1/h:m/2); Xh=zeros(h*m+1,h*n+1,2); Xh(:,:,1)=xh; Xh(:,:,2)=yh;

    S = CropImage;
    %% ---Fit---
    % inital guess [Amp,xo,wx,yo,wy,fi]
    A0 = [1,0,50,0,50,0];   
    % Define lower and upper bounds [Amp,xo,wx,yo,wy,fi]
    lb = [0,-n/2,0,-n/2,0,-pi];
    ub = [realmax('double'),n/2,(n/2)^2,n/2,(n/2)^2,pi];

    % Fit sample data
   
    [A,resnorm,res,flag,output] = lsqcurvefit(f,A0,X,S,lb,ub);
  
    disp(output); % display summary of LSQ algorithm

    centre = [yindex+A(4),xindex+A(2)];
    
    %% ---Plot Data---
    if PlotPlots == true
        %{
        % Plot 3D Data and Fitted curve
        hf1=figure(1); set(hf1,'Position',[1000 600 800 500]); 
        C=del2(f(A,Xh)); 
        mesh(xh,yh,f(A,Xh),C);
        hold on


        surface(x,y,S,'EdgeColor','none'); alpha(0.5); 
        colormap('pink'); view(-60,20); grid on; hold off
        %}
        % Plot Sample Pixels data
        hf2=figure(2); set(hf2,'Position',[20 20 800 800]); 
        subplot(4,4,[5,6,7,9,10,11,13,14,15]); imagesc(x(1,:),y(:,1),S); 
        colormap('hot');

        % Output and compare data and fitted function coefs
        text(-n/2-7,m/2+5.5,sprintf('\t Amplitude \t X-Coord \t X-Width \t Y-Coord \t Y-Width \t Angle'),'Color','black');
        text(-n/2-7,m/2+7.5,sprintf('Set \t %1.3f \t %1.3f \t %1.3f \t %1.3f \t %1.3f \t %1.3f',A0),'Color','blue');
        text(-n/2-7,m/2+9.5,sprintf('Fit \t %1.3f \t %1.3f \t %1.3f \t %1.3f \t %1.3f \t %1.3f',A),'Color','red');

        % Plot vertical and horizontal axis
        vx_h=x(1,:); vy_v=y(:,1);
        M=-tan(A(6));
        % generate points along _horizontal & _vertical axis
        vy_h=M*(vx_h-A(2))+A(4); hPoints = interp2(x,y,S,vx_h,vy_h,'nearest');
        vx_v=M*(A(4)-vy_v)+A(2); vPoints = interp2(x,y,S,vx_v,vy_v,'nearest');
      
        

        % plot lines 
        hold on; plot(A(2),A(4),'+b',vx_h,vy_h,'.r',vx_v,vy_v,'.g'); hold off;

        % Plot cross sections 
        dmin=1.1*min(S(:)); xfit=xh(1,:); hfit=A(1)*exp(-(xfit-A(2)).^2/(2*A(3)^2));
        dmax=1.1*max(S(:)); yfit=yh(:,1); vfit=A(1)*exp(-(yfit-A(4)).^2/(2*A(5)^2));
        subplot(4,4,[1,2,3]); xposh = (vx_h-A(2))/cos(A(6))+A(2); 
        plot(xposh,hPoints,'r.',xfit,hfit,'black'); grid on; axis([-n/2,n/2,dmin,dmax]);
        subplot(4,4,[8,12,16]); xposv = (vy_v-A(4))/cos(A(6))+A(4); 
        plot(vPoints,xposv,'g.',vfit,yfit,'black'); grid on; axis([dmin,dmax,-m/2,m/2]); 
        set(gca,'YDir','reverse');
    end
end
