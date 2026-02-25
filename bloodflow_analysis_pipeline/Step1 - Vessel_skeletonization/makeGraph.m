%skel = GU_extractVesselness_Skel(...
%    'F:\Vessel\99hz blood flow data\kf29_14z_roi1\kf29_14z_roi1_gauss_1umz_1p3238x.tif');
skel = readtiff('F:\Vessel\99hz blood flow data\kf29_14z_roi1\kf29_14z_roi1_gauss_1umz_1p3238x_skel.tif');
%{
f = 'F:\Vessel\99hz blood flow data\kf19_20z_roi1\kf19_20z_roi1_gauss_2umz_1p35umx_FIJIskeleton.tif';
im = readtiff(f);
skel2 = logical(readtiff(f));
%interpolate
spacing = [1.324; 1.324; 1];
[ny,nx,nz] = size(im);
[y,x,z] = ndgrid(1:ny,1:nx,1:nz);
[Y,X,Z] = ndgrid(1:1/spacing(1):ny,1:1/spacing(2):nx,1:1/spacing(3):nz); % adjust this line if you need to adjust anisotropy
skel = single(interp3(x,y,z,double(skel2),X,Y,Z,'nearest'));
clear x y z X Y Z
% cut this vessel
skel(70,66,13) = 0;
skel(14,83,17) = 0;

writetiff(uint8(skel), [f(1:end-4) '_Interp.tif']);
%}
%skel(200,369,10) = 0;

[x,y,z] = ind2sub(size(skel),find(skel == 1));
Position = [x y z];

d = pdist(Position);
d(d>1.8)=0;

dd = squareform(d);
G = graph(dd);
G.Nodes = array2table(Position,'VariableNames',{'X','Y','Z'});

figure
h = plot(G,'XData',G.Nodes.X,'YData',G.Nodes.Y,'ZData',G.Nodes.Z);
h.NodeColor = 'none';
h.EdgeColor = 'k';
h.LineWidth = 0.5;
daspect([1 1 1])
%%
% x y z from fiji
a = [273 238 82];b = [145 190 67];
% find 2 location indices to calculate shortest path
%loc1 = Position(:, 1) == 54 & Position(:, 2) == 18 & Position(:, 3) == 42;
%loc2 = Position(:, 1) == 45 & Position(:, 2) == 28 & Position(:, 3) == 42;
loc1 = Position(:, 1) == (a(2)+1) & Position(:, 2) == (a(1)+1) & Position(:, 3) == (a(3)+1);
loc2 = Position(:, 1) == (b(2)+1) & Position(:, 2) == (b(1)+1) & Position(:, 3) == (b(3)+1);
ind1 = find(loc1==1);
ind2 = find(loc2==1);
[P,L] = shortestpath(G,ind1,ind2);

% second ROI
a = [223 83 72];b = [227 16 61];
loc1 = Position(:, 1) == (a(2)+1) & Position(:, 2) == (a(1)+1) & Position(:, 3) == (a(3)+1);
loc2 = Position(:, 1) == (b(2)+1) & Position(:, 2) == (b(1)+1) & Position(:, 3) == (b(3)+1);
ind1 = find(loc1==1);
ind2 = find(loc2==1);
[P2,L2] = shortestpath(G,ind1,ind2);


% Visualize the result
highlight(h,P,'EdgeColor','r','LineWidth',1)
alpha(1)
highlight(h,P2,'EdgeColor','b','LineWidth',1)
alpha(1)
% Display the path length
disp(L)
disp(L2)

xlim([1 308])
ylim([1 308])
zlim([8 100])
set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[]);
set(gca,'XTick',[],'YTick',[],'ZTick',[]);

%set(gca,'XColor','none','YColor','none','ZColor','none')
% get splined length
%i = Position(P,:);
%[arclen,seglen] = arclength(i(:,1),i(:,2),i(:,3),'s');
%disp(arclen)

export_fig('fig5c.eps')

