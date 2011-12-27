function [p seg]= samba_mcplot3d(d, n, p,varargin)
% VERSION: Alpha 0.1
% DESCRIPTION:
%          Plot mocap joint/segment file in 3D space
%          
% USAGE:
%          p = samba_mcplotframe3D(mocap, frame, animpar);
% PARAMETERS:
%          mocap:mocap structure*
%          frame: frame number
%          anaimpar: animation structure*
% OUTPUT:
%          m2jpar parameter structure
% DETAILS:
%           Based on mcplotframe
% EXAMPLE:
%          [m2jpar]=samba_mcinitm2jpar('optitrack34})
% SEE ALSO:
%           
% DEPENDENCIES:
%       Mocap Toolbox [complete]
% CREDITS AND CITATIONS    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% AUTHORS:
%
%   Luiz NAVEDA and Marc LEMAN
%  
%       luiznaveda[at]gmail.com
%       http://www.ipem.ugent.be/
%
%      This project is supported by a grant from University of Ghent
%      and partially supported by CAPES/BRAZIL
%      Copyright (C) 2008 Ghent University
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% This function is part of the Samba Lib., a library for dance and music 
% analysis in Matlab platform.
% This library is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, see <http://www.gnu.org/licenses/>

% Plots frames of motion capture data.
%
% syntax
% function p = mcplotframe(d, n, par);
% function p = mcplotframe(d, n);
%
% input parameters
% d: MoCap data structure
% n: vector containing the numbers of the frames to be plotted
% par: animparams structure
%
% output
% p: animparams structure used in plotting the frames
% data.conn: 3d connection points [x1 y1 z1 x2 y2 z2] 
%
% examples
% p = mcplotframe(d, 1);
% p = mcplotframe(d, 500:10:600, par);
%
% comments
% If the animparams structure is not given, the function calls
% mcinitanimpar and sets the .limits field of the animparams structure
% automatically so that all the markers fit into all frames.
%
% ? Part of the Motion Capture Toolbox, Copyright ?2008,
% University of Jyvaskyla, Finland

cf=struct('FaceAlpha',0.3,'SegRadius',3,'SegFaces',6,'Gmode','square','Ground',0,'Galpha',.5,'Gcolor','red',...
    'EdgeAlpha',0,'FaceColor',[0 0 0],'Visible','on','az',44,'el',33,'FontSize',10,...
    'PlotFlag',1,'colormap','lines','Gx',[],'Gy',[],'Gz',0,'ShowMarkerNumbers',0,'ShowMarkers',0,...
    'SaveFrame',[],'hoff',1,'ShowMarkerNames',0,'TextLag',20);
if ~isempty(varargin)
    for k=1:2:length(varargin)
        cf.(varargin{k})=varargin{k+1};
    end
end

h=newplot;
set(h,'visible',cf.Visible)
h=newplot;
%figure(h)

if nargin==2
    p = mcinitanimpar;
end
data=[];
% if p.animate
%     mkdir(p.folder)
%     currdir = cd; % store current directory
%     cd(p.folder)
% end

bgcol=p.colors(1); mcol=p.colors(2); ccol=p.colors(3); tcol=p.colors(4); ncol=p.colors(5);

[X,Y,Z]=sphere;
az=p.az;
el=p.el;

% d1 = mcrotate(d, az, [0 0 1]);
% d2 = mcrotate(d1, el, [1 0 0]) %%%%%%%
d2=d;
x=d2.data(n,1:3:end);
y=d2.data(n,2:3:end);
z=d2.data(n,3:3:end);
xa=d2.data(:,1:3:end);
ya=d2.data(:,2:3:end);
za=d2.data(:,3:3:end);

if isempty(p.limits)
    % find ranges of coordinates
    if p.animate
        tmp=xa(:); maxx=prctile(tmp(find(~isnan(tmp))),99); minx=prctile(tmp(find(~isnan(tmp))),1);
        tmp=ya(:); maxy=prctile(tmp(find(~isnan(tmp))),99); miny=prctile(tmp(find(~isnan(tmp))),1);
        tmp=za(:); maxz=prctile(tmp(find(~isnan(tmp))),99); minz=prctile(tmp(find(~isnan(tmp))),1);
    else
        tmp=x(:); maxx=max(tmp(find(~isnan(tmp)))); minx=min(tmp(find(~isnan(tmp))));
        tmp=y(:); maxy=max(tmp(find(~isnan(tmp)))); miny=min(tmp(find(~isnan(tmp))));
        tmp=z(:); maxz=max(tmp(find(~isnan(tmp)))); minz=min(tmp(find(~isnan(tmp))));
    end
    midx = (maxx+minx)/2;
    midy = (maxy+miny)/2;
    midz = (maxz+minz)/2;

    scrratio = p.scrsize(1)/p.scrsize(2);
    range = max((maxx-minx)/scrratio, maxz-minz)/2;
    zrange = (maxy-miny)/2;
    % axis limits for plot
    p.limits = [midx-scrratio*1.2*range midx+scrratio*1.2*range midz-1.2*range midz+1.2*range];
end
minxx = p.limits(1);
maxxx = p.limits(2);
minzz = p.limits(3);
maxzz = p.limits(4);

for k=1:size(x,1) % main loop
    if p.animate
        clf
    else
        newplot;
        set(gcf,'Position',[50 50 p.scrsize(1) p.scrsize(2)]) ; % DVD: w=720 h=420

    end

%     axes('position', [0 0 1 1]); 
    hold on
    set(gcf, 'color', bgcol);
    view(45,45);
    colormap([ones(64,1) zeros(64,1) zeros(64,1)]);

    % plot marker-to-marker connections
    if ~isempty(p.conn)
        for m=1:size(p.conn,1)
            if x(k,p.conn(m,1))*x(k,p.conn(m,2))~=0
%                 seg{m}=plotmesh([x(k,p.conn(m,1)) x(k,p.conn(m,2)) y(k,p.conn(m,1)) y(k,p.conn(m,2)) z(k,p.conn(m,1)) z(k,p.conn(m,2))],cf)
%                 plot3([x(k,p.conn(m,1)) x(k,p.conn(m,2))],[y(k,p.conn(m,1)) y(k,p.conn(m,2))], [z(k,p.conn(m,1)) z(k,p.conn(m,2))], [ccol '-'], 'LineWidth', p.cwidth(min(p.conn(m,2),length(p.cwidth))));
                data.conn1(m,:)=[x(k,p.conn(m,1)) y(k,p.conn(m,1))  z(k,p.conn(m,1)) x(k,p.conn(m,2)) y(k,p.conn(m,2)) z(k,p.conn(m,2))];
                data.ProximalJoint{m}=d.markerName{p.conn(m,1)};
                data.DistalJoint{m}=d.markerName{p.conn(m,2)};
            end
        end
    end

    if ~isempty(p.conn)
        for m=1:size(p.conn,1)
            if x(k,p.conn(m,1))*x(k,p.conn(m,2))~=0
                seg{m}=plotmesh(data.conn1(m,:),cf);
%                 plot3([x(k,p.conn(m,1)) x(k,p.conn(m,2))],[y(k,p.conn(m,1)) y(k,p.conn(m,2))], [z(k,p.conn(m,1)) z(k,p.conn(m,2))], [ccol '-'], 'LineWidth', p.cwidth(min(p.conn(m,2),length(p.cwidth))));

            end
        end
    end
    
    % plot midpoint-to-midpoint connections
    if ~isempty(p.conn2)
        for m=1:size(p.conn2,1)
            if x(k,p.conn2(m,1))*x(k,p.conn2(m,2))*x(k,p.conn2(m,3))*x(k,p.conn2(m,4))~=0
                tmpx1 = (x(k,p.conn2(m,1))+x(k,p.conn2(m,2)))/2;
                tmpx2 = (x(k,p.conn2(m,3))+x(k,p.conn2(m,4)))/2;
                tmpz1 = (y(k,p.conn2(m,1))+y(k,p.conn2(m,2)))/2;
                tmpz2 = (y(k,p.conn2(m,3))+y(k,p.conn2(m,4)))/2;
                tmpy1 = (z(k,p.conn2(m,1))+z(k,p.conn2(m,2)))/2;
                tmpy2 = (z(k,p.conn2(m,3))+z(k,p.conn2(m,4)))/2;
                plot3([tmpx1 tmpz1 tmpx2], [tmpy1 tmpz1 tmpy2], [ccol '-'], 'LineWidth', p.cwidth);
                
            end
        end
    end
    
    % plot traces if animation
    if p.animate & ~isempty(p.trm)
        trlen = round(p.fps * p.trl);
        start=max(1,k-trlen);
        ind = start-1+find(~isnan(x(start:k)));
        for m=1:length(p.trm)
            plot(x(ind,p.trm(m)),z(ind,p.trm(m)),[tcol '-'],'Linewidth',p.cwidth);
        end
    end


    % plot markers
    for m=1:size(x,2)
        if x(k,m)~=0 && ~isnan(x(k,m)) % if marker visible
            if cf.ShowMarkers
                scatter3(x(k,m),y(k,m),z(k,m), 'filled');hold on
            end
            if p.showmnum || cf.ShowMarkerNumbers
                h2=text(x(k,m),y(k,m),z(k,m)+50,num2str(m),'FontSize',cf.FontSize);
                set(h2,'FontSize',16);
                set(h2,'Color',ncol)
            end
            % Plot marker names
            if cf.ShowMarkerNames
                
                nam=mcgetmarkername(d2);
                text(x(k,m)+cf.TextLag,y(k,m)+cf.TextLag,z(k,m),nam{m},'FontSize',cf.FontSize) ;
            end
%             axis off
        end
%         axis([minxx maxxx minzz maxzz]); axis off
    end
    if p.showfnum 
        text(minxx+0.95*(maxxx-minxx), minzz+0.97*(maxzz-minzz), num2str(k),...
            'Color',ncol,'HorizontalAlignment','Right','FontSize',p.msize,'FontWeight','bold');
    end
    if cf.Ground
        if isempty(cf.Gx)
            gr=makeground([minxx maxxx],[minzz maxzz],cf.Gz,cf.Gmode,0,cf.Gcolor);
        else
            gr=makeground(cf.Gx,cf.Gy,cf.Gz,cf.Gmode,cf.Galpha,cf.Gcolor);
        end
        axis off
    else
        axis on
        grid on
        set(gca,'DataAspectRatio',[1 1 1],...
        'PlotBoxAspectRatio',[1 1 1])
    end
    drawnow
    if cf.hoff
        hold off
    end
%     if p.animate
%         fn=['frame', sprintf('%0.4d',k),'.jpg'];
%         imwrite(frame2im(getframe),fn,'jpg','Quality',100);
%     end
end

% if p.animate
%     close
%     cd(currdir);
% end


return;


function seg=plotmesh(p,cf)
% plot mesh
% p=dataj.conn1(k,:); % rknee to rankle
% anglePoints3d(p1(1:3),p1(4:6))
% seg=segment.segm(k).r
% PLANE = createPlane(P1(1:3), N);
h=newplot;
tg=[p(1:3);p(4:6)];
euc=samba_tgaeucl(tg);
e=euc.eucl(1);

% cf.SegRadius=5;
% cf.SegFaces=10;
p1=[0 0 0]; 
p2=[0 0 e];
model=[p1;p2];
% figure(2)
[hz c1]=drawCircle([p1(1) p1(2) p1(3) cf.SegRadius 0 0]);

[hz c2]=drawCircle([p2(1) p2(2) p2(3) cf.SegRadius 0 0]); 
% close all
polyP1=c1(1:round(length(c2)/cf.SegFaces):(end-1),:);
polyP2=c2(1:round(length(c2)/cf.SegFaces):(end-1),:);
seg=[polyP1;polyP2];
   
[d,Z,tr] = procrustes(tg,model);

mx=(seg*tr.T*tr.b)+repmat(tr.c,[cf.SegFaces 1]);
% figure(2)
% newmx=[newmx; mx];
% scatter3(newmx(:,1),newmx(:,2),newmx(:,3)); hold on
% plot3(tg(:,1),tg(:,3),tg(:,3));hold on
dt = DelaunayTri(mx);
[ch v] = convexHull(dt);
% figure(1)
trisurf(ch, mx(:,1),mx(:,2),mx(:,3), 'FaceColor', cf.FaceColor,'EdgeAlpha',cf.EdgeAlpha,...
    'FaceAlpha',cf.FaceAlpha);hold on

% q = [1 0 1 0; 1 0.5 0.3 0.1];
% [pitch, roll, yaw] = quat2angle(q, 'YXZ')

function varargout = drawCircle(varargin)
%DRAWCIRCLE3D draw a 3D circle
%
%   Possible calls for the function :
%   drawCircle3d([XC YC ZC R THETA PHI])
%   drawCircle3d([XC YC ZC R THETA PHI PSI])
%   drawCircle3d([XC YC ZC R], [THETA PHI])
%   drawCircle3d([XC YC ZC R], [THETA PHI PSI])
%   drawCircle3d([XC YC ZC R], THETA, PHI)
%   drawCircle3d([XC YC ZC], R, THETA, PHI)
%   drawCircle3d([XC YC ZC R], THETA, PHI, PSI)
%   drawCircle3d([XC YC ZC], R, THETA, PHI, PSI)
%   drawCircle3d(XC, YC, ZC, R, THETA, PHI)
%   drawCircle3d(XC, YC, ZC, R, THETA, PHI, PSI)
%
%   where XC, YC, ZY are coordinate of circle center, R is the radius of he
%   circle, PHI and THETA are 3D angle of the normal to the plane
%   containing the circle (PHI between 0 and 2xPI corresponding to
%   longitude, and THETA from 0 to PI, corresponding to angle with
%   vertical).
%   
%   H = drawCircle3d(...)
%   return handle on the created LINE object
%   
%   See also:
%   circles3d
%
%   ------
%   Author: David Legland
%   e-mail: david.legland@jouy.inra.fr
%   Created: 2005-02-17
%   Copyright 2005 INRA - CEPIA Nantes - MIAJ (Jouy-en-Josas).

%   HISTORY
%   14/12/2006: allows unspecified PHI and THETA
%   04/01/2007: update doc, add todo for angle convention

%   Possible calls for the function, with number of arguments :
%   drawCircle3d([XC YC ZC R THETA PHI])            1
%   drawCircle3d([XC YC ZC R THETA PHI PSI])        1
%   drawCircle3d([XC YC ZC R], [THETA PHI])         2
%   drawCircle3d([XC YC ZC R], [THETA PHI PSI])     2
%   drawCircle3d([XC YC ZC R], THETA, PHI)          3
%   drawCircle3d([XC YC ZC], R, THETA, PHI)         4
%   drawCircle3d([XC YC ZC R], THETA, PHI, PSI)     4
%   drawCircle3d([XC YC ZC], R, THETA, PHI, PSI)    5
%   drawCircle3d(XC, YC, ZC, R, THETA, PHI)         6
%   drawCircle3d(XC, YC, ZC, R, THETA, PHI, PSI)    7


if length(varargin)==1
    % get center and radius
    circle = varargin{1};
    xc = circle(:,1);
    yc = circle(:,2);
    zc = circle(:,3);
    r  = circle(:,4);
    
    % get colatitude of normal
    if size(circle, 2)>=5
        theta = circle(:,5);
    else
        theta = zeros(size(circle, 1), 1);
    end

    % get azimut of normal
    if size(circle, 2)>=6
        phi     = circle(:,6);
    else
        phi = zeros(size(circle, 1), 1);
    end
    
    % get roll
    if size(circle, 2)==7
        psi = circle(:,7);
    else
        psi = zeros(size(circle, 1), 1);
    end
    
elseif length(varargin)==2
    % get center and radius
    circle = varargin{1};
    xc = circle(:,1);
    yc = circle(:,2);
    zc = circle(:,3);
    r  = circle(:,4);
    
    % get angle of normal
    angle = varargin{2};
    theta   = angle(:,1);
    phi     = angle(:,2);
    
    % get roll
    if size(angle, 2)==3
        psi = angle(:,3);
    else
        psi = zeros(size(angle, 1), 1);
    end

elseif length(varargin)==3    
    % get center and radius
    circle = varargin{1};
    xc = circle(:,1);
    yc = circle(:,2);
    zc = circle(:,3);
    r  = circle(:,4);
    
    % get angle of normal and roll
    theta   = varargin{2};
    phi     = varargin{3};
    psi     = zeros(size(phi, 1), 1);
    
elseif length(varargin)==4
    % get center and radius
    circle = varargin{1};
    xc = circle(:,1);
    yc = circle(:,2);
    zc = circle(:,3);
    
    if size(circle, 2)==4
        r   = circle(:,4);
        theta   = varargin{2};
        phi     = varargin{3};
        psi     = varargin{4};
    else
        r   = varargin{2};
        theta   = varargin{3};
        phi     = varargin{4};
        psi     = zeros(size(phi, 1), 1);
    end
    
elseif length(varargin)==5
    % get center and radius
    circle = varargin{1};
    xc = circle(:,1);
    yc = circle(:,2);
    zc = circle(:,3);
    r  = varargin{2};
    theta   = varargin{3};
    phi     = varargin{4};
    psi     = varargin{5};

elseif length(varargin)==6
    xc      = varargin{1};
    yc      = varargin{2};
    zc      = varargin{3};
    r       = varargin{4};
    theta   = varargin{5};
    phi     = varargin{6};
    psi     = zeros(size(phi, 1), 1);
  
elseif length(varargin)==7   
    xc      = varargin{1};
    yc      = varargin{2};
    zc      = varargin{3};
    r       = varargin{4};
    theta   = varargin{5};
    phi     = varargin{6};
    psi     = varargin{7};


    elseif length(varargin)==8   
    xc      = varargin{1};
    yc      = varargin{2};
    zc      = varargin{3};
    r       = varargin{4};
    theta   = varargin{5};
    phi     = varargin{6};
    psi     = varargin{7};
    plotf= varargin{8};
else
    error('DRAWCIRCLE3D: please specify center and radius');
end

N = 64;
t = [0:2*pi/N:2*pi*(1-1/N) 2*pi];


x = r*cos(t)';
y = r*sin(t)';
z = zeros(length(t), 1);

circle0 = [x y z];

tr      = translation3d(xc, yc, zc);
rot1    = rotationOz(-psi);
rot2    = rotationOy(-theta);
rot3    = rotationOz(-phi);
trans   = tr*rot3*rot2*rot1;

circle = transformPoint3d(circle0, trans);

% h = drawCurve3d(circle,plotf);
h=[];
if nargout>0
    varargout{1}=h;
end

if nargout>1
    varargout{2}=circle;
end

function gr=makeground(xl,yl,z,mode,alphaf,fcolor)

 switch mode
     case 'square'
         xx=[xl(1) xl(1) xl(2) xl(2)];
         yy=[yl(1) yl(2) yl(2) yl(1)];
         zz=[z z z z];
         patch(xx,yy,zz,fcolor,'FaceAlpha',alphaf,'FaceColor',fcolor,'EdgeAlpha',0)
         gr=[xx,yy,zz];
 end
 
