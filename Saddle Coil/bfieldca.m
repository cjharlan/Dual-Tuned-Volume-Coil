function [B]=bfieldca(thc,thx,ar,al,Ro,tol);

echo on
%This function computes the Biot-Savart integral for array elements
% on the surface of a cylinder.  
%
%Usage: [B]=BField(thc,thx,ar,al,Ro);
%
%
% B      = [Bx, By, Bz] at the observation point(s). (A/m^2)
%
% thc    = theta_center, the center of the array element conductor
% thx    = theta_extent, the angle subtended by the coil.
% ar     = array radius, the radius of the cylinder on which the array
%          exists
% al     = length of the array element (along z)
% Ro     = [x, y, z] coordinates of the observation point.
%          [n, NPE, NOP] n=1..3 for x,y,z positions of imaging plane
%
% ** ALL DIMENSIONS ARE METERS AND RADIANS! **
% 
echo off

mu=4*pi*1e-7;
[sxo syo szo] = size(Ro);
intfac=max([64 round(thx*180/pi)]);

if ~exist('tol','var')
    tol=1e-3;
end

for kk=1:4        % Four segments of the conductor
    switch kk
        case 1    % the +z branch 
            rxp=ar*cos(linspace(thc-thx/2,thc+thx/2,intfac));
            ryp=ar*sin(linspace(thc-thx/2,thc+thx/2,intfac));
            rzp=ones(1,intfac)*al/2;
            dlx=-ryp/ar;          % vector containing -sin(theta)
            dly=rxp/ar;           % vector containing cos(theta)
            dlz=zeros(1,intfac);  % dl vector is already normalized        
            clength = ar*thx;     % length of the conductor
        case 2    % the -z branch
            rxp=ar*cos(linspace(thc+thx/2,thc-thx/2,intfac));
            ryp=ar*sin(linspace(thc+thx/2,thc-thx/2,intfac));
            rzp=-ones(1,intfac)*al/2;
            dlx=ryp/ar;             % vector containing sin(theta)
            dly=-rxp/ar;            % vector containing -cos(theta)
            dlz=zeros(1,intfac);
            clength=ar*thx;
        case 3    % the thc-thx/2 branch
            rxp=ones(1,intfac)*ar*cos(thc-thx/2);
            ryp=ones(1,intfac)*ar*sin(thc-thx/2);
            rzp=linspace(-al/2,al/2,intfac);
            dlx=zeros(1,intfac);    % zero x
            dly=dlx;                % zero y
            dlz=ones(1,intfac);     % dz = 1;
            clength=al;
        case 4    % the thc+thx/2 branch
            rxp=ones(1,intfac)*ar*cos(thc+thx/2);
            ryp=ones(1,intfac)*ar*sin(thc+thx/2);
            rzp=linspace(al/2,-al/2,intfac);
            dlx=zeros(1,intfac);    % zero x
            dly=dlx;                % zero y
            dlz=-ones(1,intfac);    % dz = -1
            clength=al;
    end
    if sxo==szo
        if kk==1
            B=zeros(1,3);
        end
        rx  = Ro(1)-rxp;
        ry  = Ro(2)-ryp;
        rz  = Ro(3)-rzp;
        R   = sqrt(rx.^2+ry.^2+rz.^2);
        Bx  = sum(((dly.*rz-dlz.*ry)./(R.^3)).*(abs(R)>tol));
        By  = sum(((dlz.*rx-dlx.*rz)./(R.^3)).*(abs(R)>tol));
        Bz  = sum(((dlx.*ry-dly.*rx)./(R.^3)).*(abs(R)>tol));
        B   = B + [Bx By Bz]*mu/4/pi*clength/intfac;
    else
        if kk==1
            B=zeros(3,syo,szo);
        end
        Bx=zeros(syo,szo);
        By=Bx;
        Bz=Bx;
        for ii=1:intfac
            rx = squeeze(Ro(1,:,:))-rxp(ii);
            ry = squeeze(Ro(2,:,:))-ryp(ii);
            rz = squeeze(Ro(3,:,:))-rzp(ii);
            R   = sqrt((rx.^2)+(ry.^2)+(rz.^2));
            Bx  = Bx + ((dly(ii)*rz-dlz(ii)*ry)./(R.^3)).*(abs(R)>tol);
            By  = By + ((dlz(ii)*rx-dlx(ii)*rz)./(R.^3)).*(abs(R)>tol);
            Bz  = Bz + ((dlx(ii)*ry-dly(ii)*rx)./(R.^3)).*(abs(R)>tol);
        end
        B(1,:,:)=Bx*mu/4/pi*clength/intfac+squeeze(B(1,:,:));
        B(2,:,:)=By*mu/4/pi*clength/intfac+squeeze(B(2,:,:));
        B(3,:,:)=Bz*mu/4/pi*clength/intfac+squeeze(B(3,:,:));
    end
end