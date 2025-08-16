clear; close all;
addpath ./..
% fileID = fopen('M:\Nubes escuela\Lab18\lab18_1_params5.txt','r');
fileID = fopen('M:\bimFile_params5.txt','r');
fgetl(fileID);
NdP1 = textscan(fileID,'%f %f %f %f %f %f %*f %*f %*f %*f %*f %*f %*f %f %*f %*f %*f %f %*f %*f','delimiter',' ','CollectOutput',1);
NdP1 = single(cell2mat(NdP1));
normals1 = NdP1(:,4:6);
planarity1 = NdP1(:,7);
entropy1 = NdP1(:,8);
NdP1 = NdP1(:,1:3);

fileID = fopen('M:\Nubes escuela\Lab18\lab18_2_params5.txt','r');
fgetl(fileID);
NdP2 = textscan(fileID,'%f %f %f %f %f %f %*f %*f %*f %*f %*f %*f %*f %f %*f %*f %*f %f %*f %*f','delimiter',' ','CollectOutput',1);
NdP2 = single(cell2mat(NdP2));
normals2 = NdP2(:,4:6);
planarity2 = NdP2(:,7);
entropy2 = NdP2(:,8);
NdP2 = NdP2(:,1:3);

min_z = min(NdP1(:,3));

NdP_aux = NdP1;
tresh = 0.05;
slab_norm = 0.02;
n_iter = length(NdP_aux)/50;
i = 1;
tic;
while ((i < n_iter) && ~isempty(NdP_aux))
    disp(['Iter ',num2str(i),' of ',num2str(n_iter),' // ', num2str(size(NdP_aux,1)),' points left']);
    [~, idx] = max(planarity1);
    planarity1(idx) = -10000;
    p1 = NdP_aux(idx,:);
    n1 = normals1(idx,:);
    [points,~] = find((normals1(:,1)-n1(1)>-tresh & normals1(:,2)-n1(2)>-tresh & normals1(:,3)-n1(3)>-tresh) & ...
                      (normals1(:,1)-n1(1)< tresh & normals1(:,2)-n1(2)< tresh & normals1(:,3)-n1(3)< tresh));
    PoI = NdP_aux(points,1:3);
    
    A = n1(1);
    B = n1(2);
    C = n1(3);
    D = -p1*n1';
    dist_plane = abs(PoI*n1' + D)/(sqrt(A^2 + B^2 + C^2));
    aux = dist_plane<slab_norm;
    plano1 = NdP_aux(points(aux),:);
    if length(plano1)>1000
        [n,V,p] = affine_fit(double(plano1));
        plane_array1{i} = plano1;
        points_array1{i} = p;
        normals_array1{i} = n;
        NdP_aux(points(aux),:) = [];
        normals1(points(aux),:) = [];
        planarity1(points(aux)) = [];
    end
    i = i+1;
end
time_plane_segm = toc;
disp(['Plane segmentation elapsed time: ',num2str(time_plane_segm),' seconds']);
%# find empty cells
emptyCells = cellfun(@isempty,plane_array1);
%# remove empty cells
plane_array1(emptyCells) = [];
emptyCells = cellfun(@isempty,points_array1);
points_array1(emptyCells) = [];
emptyCells = cellfun(@isempty,normals_array1);
normals_array1(emptyCells) = [];

% sort the cells arrays by planes with more number of points inside
s = cellfun(@size,plane_array1,'uniform',false);
[~,is] = sortrows(cat(1,s{:}),-[1 2]);
plane_array1 = plane_array1(is);
points_array1 = points_array1(is);
normals_array1 = normals_array1(is);

%% flip all the normals towars the viewpoint vp = (0,0,1.4)
vp = [0 ; 0 ; min_z + 1.4];
for i = 1:length(normals_array1)
    aux = normals_array1{i}' * (vp - points_array1{i}');
    if aux < 0
        normals_array1{i} = -1 * normals_array1{i};
    end
end

figure(1);
hold off;
% plot3(NdP_aux(:,1),NdP_aux(:,2),NdP_aux(:,3),'.','MarkerSize',1);
grid on
hold on
colores = rand(length(plane_array1),3);
for i=1:length(plane_array1)
    plot3(plane_array1{i}(:,1),plane_array1{i}(:,2),plane_array1{i}(:,3),'.','Color',colores(i,:),'MarkerSize',1);
    quiver3(points_array1{i}(1),points_array1{i}(2),points_array1{i}(3),normals_array1{i}(1),normals_array1{i}(2),normals_array1{i}(3),'k');
end
hold off;
xlabel('X coordinate (m.)');
ylabel('Y coordinate (m.)');
zlabel('Z coordinate (m.)');

%%
min_z = min(NdP2(:,3));
NdP_aux = NdP2;
n_iter = length(NdP_aux)/50;
i = 1;
tic;
while ((i < n_iter) && ~isempty(NdP_aux))
    [~, idx] = max(planarity2);
    planarity2(idx) = -10000;
    p1 = NdP_aux(idx,:);
    n1 = normals2(idx,:);
    [points,~] = find((normals2(:,1)-n1(1)>-tresh & normals2(:,2)-n1(2)>-tresh & normals2(:,3)-n1(3)>-tresh) & ...
                      (normals2(:,1)-n1(1)< tresh & normals2(:,2)-n1(2)< tresh & normals2(:,3)-n1(3)< tresh));
    PoI = NdP_aux(points,1:3);
    
    A = n1(1);
    B = n1(2);
    C = n1(3);
    D = -p1*n1';
    dist_plane = abs(PoI*n1' + D)/(sqrt(A^2 + B^2 + C^2));
    aux = dist_plane<slab_norm;
    plano1 = NdP_aux(points(aux),:);
    if length(plano1)>1000
        [n,V,p] = affine_fit(double(plano1));
        plane_array2{i} = plano1;
        points_array2{i} = p;
        normals_array2{i} = n;
        NdP_aux(points(aux),:) = [];
        normals2(points(aux),:) = [];
        planarity2(points(aux)) = [];
    end
    i = i+1;
end
time_plane_segm = toc;
disp(['Plane segmentation elapsed time: ',num2str(time_plane_segm),' seconds']);
%# find empty cells
emptyCells = cellfun(@isempty,plane_array2);
%# remove empty cells
plane_array2(emptyCells) = [];
emptyCells = cellfun(@isempty,points_array2);
points_array2(emptyCells) = [];
emptyCells = cellfun(@isempty,normals_array2);
normals_array2(emptyCells) = [];

% sort the cells arrays by planes with more number of points inside
s = cellfun(@size,plane_array2,'uniform',false);
[~,is] = sortrows(cat(1,s{:}),-[1 2]);
plane_array2 = plane_array2(is);
points_array2 = points_array2(is);
normals_array2 = normals_array2(is);

%% flip all the normals towars the viewpoint vp = (0,0,1.4)
vp = [0 ; 0 ; min_z + 1.4];
for i = 1:length(normals_array2)
    aux = normals_array2{i}' * (vp - points_array2{i}');
    if aux < 0
        normals_array2{i} = -1 * normals_array2{i};
    end
end

figure(2);
hold off;
% plot3(NdP_aux(:,1),NdP_aux(:,2),NdP_aux(:,3),'.','MarkerSize',1);
grid on
hold on
colores = rand(length(plane_array2),3);
for i=1:length(plane_array2)
    plot3(plane_array2{i}(:,1),plane_array2{i}(:,2),plane_array2{i}(:,3),'.','Color',colores(i,:),'MarkerSize',1);
    quiver3(points_array2{i}(1),points_array2{i}(2),points_array2{i}(3),normals_array2{i}(1),normals_array2{i}(2),normals_array2{i}(3),'k');
end
hold off;
xlabel('X coordinate (m.)');
ylabel('Y coordinate (m.)');
zlabel('Z coordinate (m.)');

tic; trans = compute_planar_reg(normals_array1, normals_array2, plane_array1, plane_array2, points_array1, points_array2, 3, 0); toc

[B, I] = sort([trans.rmse],'ascend');
best10_idx = I(1:10);
best_trans = [trans(best10_idx).Tmax];
best_rmse = [trans(best10_idx).rmse];


