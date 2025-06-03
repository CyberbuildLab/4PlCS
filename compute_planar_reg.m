function trans = compute_planar_reg(normals_array1, normals_array2, plane_array1, plane_array2, points_array1, points_array2, area_array1, area_array2, num_planes, sup_thresh, vis_flag)
%% compute the final distances and angles between the planes of point cloud 1
plane_area1 = cell2mat(area_array1);
planes_weight = plane_area1/norm(plane_area1);
proy_dist1 = zeros(size(points_array1,2),size(points_array1,2),1);
angle_dist1 = zeros(size(points_array1,2));
parallel_planes1 = zeros(size(points_array1,2));
NaN_aux=[];
NaN_aux(:,:) = eye(size(points_array1,2),size(points_array1,2));
NaN_aux(NaN_aux==1) = NaN;
for i = 1:size(points_array1,2)
    for j = 1:size(points_array1,2)
        %         mat_dist1(i,j) = sqrt((points_array1{i}(1)-points_array1{j}(1))^2 + (points_array1{i}(2)-points_array1{j}(2))^2 + (points_array1{i}(3)-points_array1{j}(3))^2);
        angle_dist1(i,j) = real(acosd(dot(normals_array1{i},normals_array1{j})));
        %         A1 = normals_array1{i}(1);
        %         B1 = normals_array1{i}(2);
        %         C1 = normals_array1{i}(3);
        %         D1 = -points_array1{i}*normals_array1{i};
        A2 = normals_array1{j}(1);
        B2 = normals_array1{j}(2);
        C2 = normals_array1{j}(3);
        D2 = -points_array1{j}*normals_array1{j};
        proy_dist1(i,j) = abs(points_array1{i}*normals_array1{j} + D2)/(sqrt(A2^2 + B2^2 + C2^2));
        if (angle_dist1(i,j) <= 5 || angle_dist1(i,j) >= 175) && i ~= j
            parallel_planes1(i,j) = 1;
        else
            parallel_planes1(i,j) = 0;
        end
        
    end
    
end
% mat_dist1(:,:) = mat_dist1(:,:) + NaN_aux;
angle_dist1 = angle_dist1 + NaN_aux;
% parallel_planes1 = parallel_planes1 + NaN_aux;
proy_dist1 = proy_dist1 + NaN_aux;

%% compute the final distances and angles between the planes of point cloud 2
% mat_dist2 = zeros(size(points_array2,2),size(points_array2,2),1);
plane_area2 = cell2mat(area_array2);
proy_dist2 = zeros(size(points_array2,2),size(points_array2,2),1);
angle_dist2 = zeros(size(points_array2,2));
parallel_planes2 = zeros(size(points_array2,2));
NaN_aux=[];
NaN_aux(:,:) = eye(size(points_array2,2),size(points_array2,2));
NaN_aux(NaN_aux==1) = NaN;
for i = 1:size(points_array2,2)
    for j = 1:size(points_array2,2)
        %         mat_dist2(i,j) = sqrt((points_array2{i}(1)-points_array2{j}(1))^2 + (points_array2{i}(2)-points_array2{j}(2))^2 + (points_array2{i}(3)-points_array2{j}(3))^2);
        angle_dist2(i,j) = real(acosd(dot(normals_array2{i},normals_array2{j})));
        %         A1 = normals_array2{i}(1);
        %         B1 = normals_array2{i}(2);
        %         C1 = normals_array2{i}(3);
        %         D1 = -points_array2{i}*normals_array2{i};
        A2 = normals_array2{j}(1);
        B2 = normals_array2{j}(2);
        C2 = normals_array2{j}(3);
        D2 = -points_array2{j}*normals_array2{j};
        proy_dist2(i,j) = abs(points_array2{i}*normals_array2{j} + D2)/(sqrt(A2^2 + B2^2 + C2^2));
        if (angle_dist2(i,j) <= 5 || angle_dist2(i,j) >= 175) && i ~= j
            parallel_planes2(i,j) = 1;
        else
            parallel_planes2(i,j) = 0;
        end
        
    end
    
end
% mat_dist2(:,:) = mat_dist2(:,:) + NaN_aux;
angle_dist2 = angle_dist2 + NaN_aux;
% parallel_planes2 = parallel_planes2 + NaN_aux;
proy_dist2 = proy_dist2 + NaN_aux;

%%
iter = 0;
cont = 1;
normals1 = cell2mat(normals_array1)';
normals2 = cell2mat(normals_array2)';
trans = struct('Tmax', [], 'Supporting_planes', [], 'Supporting_points', [], 'RMSE', [], 'Corresponding_planes', []);
% normals1 = cell2mat(normals_array1)';
% normals2 = cell2mat(normals_array2)';
% centroids1 = reshape(cell2mat(points_array1),3,[])';
centroids2 = reshape(cell2mat(points_array2),3,[])';
% NS = KDTreeSearcher(normals1');
% NS = KDTreeSearcher(centroids1);
tested_bases = 0;

if num_planes == 3
    
    model_congruent_matrix = [NaN NaN NaN];
    cant_cand = 0;
    plane_sup_bases = 0;
    centroid_sup_bases = 0;
    time_plane = 0;
    time_angle = 0;
    time_point = 0;
    time_rmse = 0;
    while iter < 100
        
        plane1_idx = 1;
        plane2_idx = 1;
        plane3_idx = 1;
        coincident_flag = 0;
        normalZ = [0 ; 0 ; 1];
        
        while (length(unique([plane1_idx plane2_idx plane3_idx])) < 3 ) && (coincident_flag == 0) && ~ismember(sort([plane1_idx plane2_idx plane3_idx]), model_congruent_matrix,'rows')
            plane1_idx = randsample(size(angle_dist1,1), 1, true, planes_weight);
            pos = find(angle_dist1(plane1_idx,:) >= 10 & angle_dist1(plane1_idx,:) <= 170);
            plane2_idx = pos(randsample(length(pos), 1, true, planes_weight(pos)));
            
            pos = find(parallel_planes1(plane1_idx,:)==0 & parallel_planes1(plane2_idx,:)==0);
            plane3_idx = pos(randsample(length(pos), 1, true, planes_weight(pos)));
            
            coincident_flag = (proy_dist1(plane1_idx, plane3_idx) > 0.2) && (proy_dist1(plane2_idx, plane3_idx) > 0.2) && (proy_dist1(plane1_idx, plane2_idx) > 0.2);
            iter = iter + 1;
            model_congruent_matrix = [model_congruent_matrix ; sort([plane1_idx plane2_idx plane3_idx])];
            
        end
        
        if vis_flag
            figure(1);
            hold off;
            plot3(plane_array1{plane1_idx}(:,1),plane_array1{plane1_idx}(:,2),plane_array1{plane1_idx}(:,3),'b');
            hold on;
            plot3(plane_array1{plane2_idx}(:,1),plane_array1{plane2_idx}(:,2),plane_array1{plane2_idx}(:,3),'r');
            plot3(plane_array1{plane3_idx}(:,1),plane_array1{plane3_idx}(:,2),plane_array1{plane3_idx}(:,3),'g');
            grid on;
        end
        
        ang12 = angle_dist1(plane1_idx, plane2_idx);
        ang13 = angle_dist1(plane1_idx, plane3_idx);
        ang23 = angle_dist1(plane2_idx, plane3_idx);
        ang1Z = real(acosd(dot(normals_array1{plane1_idx},normalZ)));
        ang2Z = real(acosd(dot(normals_array1{plane2_idx},normalZ)));
        ang3Z = real(acosd(dot(normals_array1{plane3_idx},normalZ)));
        
        [row12,col12] = find(angle_dist2 >= (ang12-5) & angle_dist2 <= (ang12+5));
%         area1_idx = find(plane_area2/plane_area1(plane1_idx)<= 1.3 & plane_area2/plane_area1(plane1_idx)>= 0.001);
%         area1_idx = find(plane_area2/plane_area1(plane1_idx)<= 1.5 & plane_area2/plane_area1(plane1_idx)>= 0.5); %casita & casota
        area1_idx = find(plane_area2/plane_area1(plane1_idx)<= 7 & plane_area2/plane_area1(plane1_idx)>= 0.13); %mercury
%         area1_idx = find(plane_area2/plane_area1(plane1_idx)<= 1.5 & plane_area2/plane_area1(plane1_idx)>= 0.5); %steel struct
%         area1_idx = find(plane_area2/plane_area1(plane1_idx)<= 36 & plane_area2/plane_area1(plane1_idx)>= 0.02); %waterloo
        area2_idx = find(plane_area2/plane_area1(plane2_idx)<= 7 & plane_area2/plane_area1(plane2_idx)>= 0.13);
        cols12 = col12(ismember(row12, area2_idx) & ismember(col12, area1_idx));
        rows12 = row12(ismember(row12, area2_idx) & ismember(col12, area1_idx));
        
        aux_mat = sort([cols12 rows12],2);
        [~, index] = unique(aux_mat,'rows','first');  %# Finds indices of unique rows
        repeatedIndex = setdiff(1:size(aux_mat,1),index);
        cols12(repeatedIndex) = [];
        rows12(repeatedIndex) = [];
        cant_cand = cant_cand + length(rows12);
        for i = 1:length(rows12)
            disp(['Iteration ', num2str(i),' of ', num2str(length(rows12))]);
            cand1 = cols12(i);
            cand2 = rows12(i);
            
            [~,col3] = find(angle_dist2(cand2,:) >= (ang23-5) & angle_dist2(cand2,:) <= (ang23+5) & angle_dist2(cand1,:) >= (ang13-5) & angle_dist2(cand1,:) <= (ang13+5));
            area3_idx = find(plane_area2/plane_area1(plane3_idx)<= 36 & plane_area2/plane_area1(plane3_idx)>= 0.02);
            col3 = col3(ismember(col3, area3_idx));
            cant_cand = cant_cand + length(col3);
            for j = 1:length(col3)
                
                cand3 = col3(j);
                
                time1 = tic;
                if (length(unique([cand1 cand2 cand3])) == 3) && ...
                        (angle_dist2(cand1,cand2)>=(ang12-5)  && angle_dist2(cand1,cand2)<=(ang12+5)) && ...
                        (angle_dist2(cand1,cand3)>=(ang13-5)  && angle_dist2(cand1,cand3)<=(ang13+5)) && ...
                        (angle_dist2(cand2,cand3)>=(ang23-5)  && angle_dist2(cand2,cand3)<=(ang23+5)) && ...
                        (real(acosd(dot(normals_array2{cand1},normalZ)))>=(ang1Z-5)  && real(acosd(dot(normals_array2{cand1},normalZ)))<=(ang1Z+5)) && ...
                        (real(acosd(dot(normals_array2{cand2},normalZ)))>=(ang2Z-5)  && real(acosd(dot(normals_array2{cand2},normalZ)))<=(ang2Z+5)) && ...
                        (real(acosd(dot(normals_array2{cand3},normalZ)))>=(ang3Z-5)  && real(acosd(dot(normals_array2{cand3},normalZ)))<=(ang3Z+5))
                    time_angle = time_angle + toc(time1);
                    if vis_flag
                        figure(2);
                        hold off;
                        plot3(plane_array2{cand1}(:,1),plane_array2{cand1}(:,2),plane_array2{cand1}(:,3),'.b','MarkerSize',1);
                        hold on;
                        plot3(plane_array2{cand2}(:,1),plane_array2{cand2}(:,2),plane_array2{cand2}(:,3),'.r','MarkerSize',1);
                        plot3(plane_array2{cand3}(:,1),plane_array2{cand3}(:,2),plane_array2{cand3}(:,3),'.g','MarkerSize',1);
                        grid on;
                    end
                    tested_bases = tested_bases + 1;
                    n11 = normals_array1{plane1_idx};
                    n21 = normals_array2{cand1};
                    n12 = normals_array1{plane2_idx};
                    n22 = normals_array2{cand2};
                    n13 = normals_array1{plane3_idx};
                    n23 = normals_array2{cand3};
                    
                    [~,~,transform] = procrustes([n11' ; n12' ; n13' ; normalZ'], [n21' ; n22' ; n23' ; normalZ'],'scaling',0,'reflection',0);
                    
                    R = transform.T;
                    
                    if R(3,3) >=0.9
                        
                        n1 = normals_array1{plane1_idx};
                        n2 = normals_array1{plane2_idx};
                        n3 = normals_array1{plane3_idx};
                        P1_aux = points_array1{plane1_idx};
                        P2_aux = points_array1{plane2_idx};
                        P3_aux = points_array1{plane3_idx};
                        A = [n1(1) n1(2) n1(3) ; n2(1) n2(2) n2(3) ; n3(1) n3(2) n3(3)];
                        B = [P1_aux*n1 ; P2_aux*n2 ; P3_aux*n3];
                        P1 = A\B;
                        P1 = P1';
                        
                        n1 = normals_array2{cand1};
                        n2 = normals_array2{cand2};
                        n3 = normals_array2{cand3};
                        P1_aux = points_array2{cand1};
                        P2_aux = points_array2{cand2};
                        P3_aux = points_array2{cand3};
                        A = [n1(1) n1(2) n1(3) ; n2(1) n2(2) n2(3) ; n3(1) n3(2) n3(3)];
                        B = [P1_aux*n1 ; P2_aux*n2 ; P3_aux*n3];
                        P2 = A\B;
                        
                        P2a = P2' * R;
                        T = P1-P2'- (P2a-P2');
                        
                        time2 = tic;
                        centroids2_reg = centroids2 * R + repmat(T,size(centroids2,1),1);
                        normals2_reg = normals2 * R;
                        supp_plane_dist = zeros(size(plane_array2,2),5);
                        supp_plane_idx = zeros(size(plane_array2,2),5);
                        supp_plane_orient = zeros(size(plane_array2,2),5);
                        for yy = 1:size(plane_array2,2)
                            dist = zeros(size(plane_array1,2),1);
                            for xx = 1:size(plane_array1,2)
                                
                                A1 = normals_array1{xx}(1);
                                B1 = normals_array1{xx}(2);
                                C1 = normals_array1{xx}(3);
                                D1 = -points_array1{xx}*normals_array1{xx};
                                dist(xx) = abs(A1*centroids2_reg(yy,1) + B1*centroids2_reg(yy,2) + C1*centroids2_reg(yy,3) + D1)/sqrt(A1*A1 + B1*B1 + C1*C1);
                                
                            end
%                             [aux_dist, aux_idx] = sort(dist,'ascend');
%                             supp_plane_dist(yy,:) = aux_dist(1:5);
%                             supp_plane_idx(yy,:) = aux_idx(1:5);
                            orient = (normals1(:,1)-normals2_reg(yy,1)>-0.1 & normals1(:,2)-normals2_reg(yy,2)>-0.1 & normals1(:,3)-normals2_reg(yy,3)>-0.1) & (normals1(:,1)-normals2_reg(yy,1)< 0.1 & normals1(:,2)-normals2_reg(yy,2)< 0.1 & normals1(:,3)-normals2_reg(yy,3)< 0.1);
                            [aux_dist, aux_idx] = sort(dist,'ascend');
                            orient = orient(aux_idx);
                            supp_plane_dist(yy,:) = aux_dist(1:5);
                            supp_plane_idx(yy,:) = aux_idx(1:5);
                            supp_plane_orient(yy,:) = orient(1:5);
                        end
%                         n_supporting_planes = sum(supp_plane_dist(:,1)<=sup_thresh);
                        n_supporting_planes = sum(supp_plane_dist(:,1)<=sup_thresh & supp_plane_orient(:,1)==1);
                        time_plane = time_plane + toc(time2);
                        
                        if  (n_supporting_planes*100/size(points_array2,2))>=20
                            
                            time3 = tic;
                            plane_sup_bases = plane_sup_bases +1;
                            centroids = [points_array2{cand1} ; points_array2{cand2} ; points_array2{cand3}];
                            centroids_reg = centroids * R + repmat(T,size(centroids,1),1);
                            alfa = acosd(dot(normalZ,normals_array1{plane1_idx})/sqrt(normals_array1{plane1_idx}(1)*normals_array1{plane1_idx}(1) + normals_array1{plane1_idx}(2)*normals_array1{plane1_idx}(2) + normals_array1{plane1_idx}(3)*normals_array1{plane1_idx}(3)));
                            RR = [cosd(alfa) -sind(alfa) 0 ; sind(alfa) cosd(alfa) 0 ; 0 0 cosd(alfa)+(1-cosd(alfa))];
                            plane_aux = plane_array1{plane1_idx} * RR;
                            coord_idx = var(plane_aux) ~= min(var(plane_aux));
                            plane_aux1 = sortrows(plane_aux(:,coord_idx ));
                            borde = boundary(double(plane_aux1));
                            plane_aux = centroids_reg(1,:) * RR;
                            plane_aux2 = plane_aux(:,coord_idx );
                            in1 = inpolygon(plane_aux2(:,1),plane_aux2(:,2),plane_aux1(borde,1),plane_aux1(borde,2));
                            
                            alfa = acosd(dot(normalZ,normals_array1{plane2_idx})/sqrt(normals_array1{plane2_idx}(1)*normals_array1{plane2_idx}(1) + normals_array1{plane2_idx}(2)*normals_array1{plane2_idx}(2) + normals_array1{plane2_idx}(3)*normals_array1{plane2_idx}(3)));
                            RR = [cosd(alfa) -sind(alfa) 0 ; sind(alfa) cosd(alfa) 0 ; 0 0 cosd(alfa)+(1-cosd(alfa))];
                            plane_aux = plane_array1{plane2_idx} * RR;
                            coord_idx = var(plane_aux) ~= min(var(plane_aux));
                            plane_aux1 = sortrows(plane_aux(:,coord_idx));
                            borde = boundary(double(plane_aux1));
                            plane_aux = centroids_reg(2,:) * RR;
                            plane_aux2 = plane_aux(:,coord_idx);
                            in2 = inpolygon(plane_aux2(:,1),plane_aux2(:,2),plane_aux1(borde,1),plane_aux1(borde,2));
                            
                            alfa = acosd(dot(normalZ,normals_array1{plane3_idx})/sqrt(normals_array1{plane3_idx}(1)*normals_array1{plane3_idx}(1) + normals_array1{plane3_idx}(2)*normals_array1{plane3_idx}(2) + normals_array1{plane3_idx}(3)*normals_array1{plane3_idx}(3)));
                            RR = [cosd(alfa) -sind(alfa) 0 ; sind(alfa) cosd(alfa) 0 ; 0 0 cosd(alfa)+(1-cosd(alfa))];
                            plane_aux = plane_array1{plane3_idx} * RR;
                            coord_idx = var(plane_aux) ~= min(var(plane_aux));
                            plane_aux1 = sortrows(plane_aux(:,coord_idx));
                            borde = boundary(double(plane_aux1));
                            plane_aux = centroids_reg(3,:) * RR;
                            plane_aux2 = plane_aux(:,coord_idx);
                            in3 = inpolygon(plane_aux2(:,1),plane_aux2(:,2),plane_aux1(borde,1),plane_aux1(borde,2));
                            time_point = time_point + toc(time3);
                            
                            if sum([in1 in2 in3]) == 3
                                
                                centroid_sup_bases = centroid_sup_bases + 1;
                                Tmax = [R' T' ; 0 0 0 1];
                                
                                time4 = tic;
                                plane_array_aux = plane_array2;
                                supp_plane_idx = supp_plane_idx(supp_plane_dist(:,1)<=sup_thresh,:);
                                points2 = plane_array_aux(supp_plane_dist(:,1)<=sup_thresh);
                                dist = cell(1,size(supp_plane_idx,1));
                                acum_points_in = 0;
                                acum_points_tot = 0;
                                
                                for zz = 1 : size(supp_plane_idx,1)
                                    
                                    points_reg = points2{zz} * R + repmat(T,length(points2{zz}),1);
                                    in = zeros(5,size(points_reg,1));
                                    
                                    for tt = 1:5
                                        
                                        A1 = normals_array1{supp_plane_idx(zz,tt)}(1);
                                        B1 = normals_array1{supp_plane_idx(zz,tt)}(2);
                                        C1 = normals_array1{supp_plane_idx(zz,tt)}(3);
                                        alfa = acosd(dot(normalZ,normals_array1{supp_plane_idx(zz,tt)})/sqrt(A1*A1 + B1*B1 + C1*C1));
                                        RR = [cosd(alfa) -sind(alfa) 0 ; sind(alfa) cosd(alfa) 0 ; 0 0 cosd(alfa)+(1-cosd(alfa))];
                                        plane_aux = plane_array1{supp_plane_idx(zz,tt)} * RR;
                                        coord_idx = var(plane_aux) ~= min(var(plane_aux));
                                        plane_aux1 = sortrows(plane_aux(:,coord_idx));
                                        borde = boundary(double(plane_aux1));
                                        plane_aux = points_reg * RR;
                                        plane_aux2 = plane_aux(:,coord_idx);
                                        in(tt,:) = inpolygon(plane_aux2(:,1),plane_aux2(:,2),plane_aux1(borde,1),plane_aux1(borde,2));
                                        %                                         cont_in(tt) = sum(in);
                                        
                                    end
                                    
                                    [~, in_idx] = max(sum(in,2));
                                    A1 = normals_array1{supp_plane_idx(zz,1)}(1);
                                    B1 = normals_array1{supp_plane_idx(zz,1)}(2);
                                    C1 = normals_array1{supp_plane_idx(zz,1)}(3);
                                    D1 = -points_array1{supp_plane_idx(zz,1)}*normals_array1{supp_plane_idx(zz,1)};
                                    dist{zz} = abs(A1*points_reg(:,1) + B1*points_reg(:,2) + C1*points_reg(:,3) + D1)/sqrt(A1*A1 + B1*B1 + C1*C1);
                                    dist{zz}(in(in_idx,:)==false) = dist{zz}(in(in_idx,:)==false)*NaN;
                                    acum_points_in = acum_points_in + sum(in(in_idx,:));
                                    acum_points_tot = acum_points_tot + size(points_reg,1);
                                    
                                end
                                dist = cell2mat(dist');
                                dist(isnan(dist)) = [];
                                rmse = mean(dist);
                                time_rmse = time_rmse + toc(time4);
                                points_in = acum_points_in*100/acum_points_tot;
                                trans(cont) = struct('Tmax', Tmax, 'Supporting_planes', n_supporting_planes, 'Supporting_points', points_in, 'RMSE', rmse, 'Corresponding_planes', [plane1_idx plane2_idx plane3_idx cand1 cand2 cand3]);
                                cont = cont + 1;
                                disp(['Cont : ',num2str(cont)]);
                                
                            end
                        end
                    end
                end
            end
        end
    end
    
elseif num_planes == 4
    
    model_congruent_matrix = [NaN NaN NaN NaN];
    cant_cand = 0;
    plane_sup_bases = 0;
    centroid_sup_bases = 0;
    time_plane = 0;
    time_angle = 0;
    time_point = 0;
    time_rmse = 0;
    while iter < 100
        
        plane1_idx = 1;
        plane2_idx = 1;
        plane3_idx = 1;
        plane4_idx = 1;
        coincident_flag = 0;
        normalZ = [0 ; 0 ; 1];
        
        while (length(unique([plane1_idx plane2_idx plane3_idx plane4_idx])) < 4 ) && (coincident_flag == 0) && ~ismember(sort([plane1_idx plane2_idx plane3_idx plane4_idx]), model_congruent_matrix,'rows')
            plane1_idx = randsample(size(angle_dist1,1), 1, true, planes_weight);
            pos = find(angle_dist1(plane1_idx,:) >= 10 & angle_dist1(plane1_idx,:) <= 170);
            plane2_idx = pos(randsample(length(pos), 1, true, planes_weight(pos)));
            
            pos = find(parallel_planes1(plane1_idx,:)==0 & parallel_planes1(plane2_idx,:)==0);
            plane3_idx = pos(randsample(length(pos), 1, true, planes_weight(pos)));
            
            %             pos = find(parallel_planes1(plane1_idx,:) & proy_dist1(plane1_idx,:) >= 1);
            pos = find(angle_dist1(plane3_idx,:) >= 10 & angle_dist1(plane3_idx,:) <= 170);
            if ~isempty(pos)
                plane4_idx = pos(randsample(length(pos), 1, true, planes_weight(pos)));
            else
                plane4_idx = plane1_idx;
            end
            
            coincident_flag = (proy_dist1(plane1_idx, plane3_idx) > 0.2) && (proy_dist1(plane1_idx, plane4_idx) > 0.2) && (proy_dist1(plane2_idx, plane3_idx) > 0.2) && (proy_dist1(plane2_idx, plane4_idx) > 0.2 && (length(unique([plane1_idx plane2_idx plane3_idx plane4_idx])) < 4));
            iter = iter + 1;
            model_congruent_matrix = [model_congruent_matrix ; sort([plane1_idx plane2_idx plane3_idx plane4_idx])];
        end
        
        if vis_flag
            figure(1);
            hold off;
            plot3(plane_array1{plane1_idx}(:,1),plane_array1{plane1_idx}(:,2),plane_array1{plane1_idx}(:,3),'b');
            hold on;
            plot3(plane_array1{plane2_idx}(:,1),plane_array1{plane2_idx}(:,2),plane_array1{plane2_idx}(:,3),'r');
            plot3(plane_array1{plane3_idx}(:,1),plane_array1{plane3_idx}(:,2),plane_array1{plane3_idx}(:,3),'g');
            plot3(plane_array1{plane4_idx}(:,1),plane_array1{plane4_idx}(:,2),plane_array1{plane4_idx}(:,3),'y');
            grid on;
        end
        
        ang12 = angle_dist1(plane1_idx, plane2_idx);
        ang13 = angle_dist1(plane1_idx, plane3_idx);
        ang14 = angle_dist1(plane1_idx, plane4_idx);
        ang23 = angle_dist1(plane2_idx, plane3_idx);
        ang24 = angle_dist1(plane2_idx, plane4_idx);
        ang34 = angle_dist1(plane3_idx, plane4_idx);
        ang1Z = real(acosd(dot(normals_array1{plane1_idx},normalZ)));
        ang2Z = real(acosd(dot(normals_array1{plane2_idx},normalZ)));
        ang3Z = real(acosd(dot(normals_array1{plane3_idx},normalZ)));
        ang4Z = real(acosd(dot(normals_array1{plane4_idx},normalZ)));
        
        [row12,col12] = find(angle_dist2 >= (ang12-5) & angle_dist2 <= (ang12+5));
%         area1_idx = find(plane_area2/plane_area1(plane1_idx)<= 1.3 & plane_area2/plane_area1(plane1_idx)>= 0.001);
%         area1_idx = find(plane_area2/plane_area1(plane1_idx)<= 1.5 & plane_area2/plane_area1(plane1_idx)>= 0.5); %casita & casota
        area1_idx = find(plane_area2/plane_area1(plane1_idx)<= 7 & plane_area2/plane_area1(plane1_idx)>= 0.13); %mercury
%         area1_idx = find(plane_area2/plane_area1(plane1_idx)<= 1.5 & plane_area2/plane_area1(plane1_idx)>= 0.5); %steel struct
%         area1_idx = find(plane_area2/plane_area1(plane1_idx)<= 36 & plane_area2/plane_area1(plane1_idx)>= 0.02); %waterloo
        area2_idx = find(plane_area2/plane_area1(plane2_idx)<= 7 & plane_area2/plane_area1(plane2_idx)>= 0.13);
        cols12 = col12(ismember(row12, area2_idx) & ismember(col12, area1_idx));
        rows12 = row12(ismember(row12, area2_idx) & ismember(col12, area1_idx));
        
        aux_mat = sort([cols12 rows12],2);
        [~, index] = unique(aux_mat,'rows','first');  %# Finds indices of unique rows
        repeatedIndex = setdiff(1:size(aux_mat,1),index);
        cols12(repeatedIndex) = [];
        rows12(repeatedIndex) = [];
        
        %         toalign_congruent_sets = [NaN NaN NaN NaN];
        cant_cand = cant_cand + length(rows12);
        for i = 1:length(rows12)
            
            disp(['Iteration ', num2str(i),' of ', num2str(length(rows12))]);
            cand1 = cols12(i);
            cand2 = rows12(i);
            
            [row34,col34] = find(angle_dist2 >= (ang34-5) & angle_dist2 <= (ang34+5));
            area3_idx = find(plane_area2/plane_area1(plane3_idx)<= 7 & plane_area2/plane_area1(plane3_idx)>= 0.13);
            area4_idx = find(plane_area2/plane_area1(plane4_idx)<= 7 & plane_area2/plane_area1(plane4_idx)>= 0.13);
            cols34 = col34(ismember(row34, area4_idx) & ismember(col34, area3_idx) & ~ismember(row34, cand1) & ~ismember(row34, cand2) & ~ismember(col34, cand1) & ~ismember(col34, cand2));
            rows34 = row34(ismember(row34, area4_idx) & ismember(col34, area3_idx) & ~ismember(row34, cand1) & ~ismember(row34, cand2) & ~ismember(col34, cand1) & ~ismember(col34, cand2));
            cant_cand = cant_cand + length(rows34);
            for j = 1:length(rows34)
                
                cand3 = cols34(j);
                cand4 = rows34(j);
                
                time1 = tic;
                if (length(unique([cand1 cand2 cand3 cand4])) == 4) && ...
                        (angle_dist2(cand1,cand2)>=(ang12-5)  && angle_dist2(cand1,cand2)<=(ang12+5)) && ...
                        (angle_dist2(cand1,cand3)>=(ang13-5)  && angle_dist2(cand1,cand3)<=(ang13+5)) && ...
                        (angle_dist2(cand1,cand4)>=(ang14-5)  && angle_dist2(cand1,cand4)<=(ang14+5)) && ...
                        (angle_dist2(cand2,cand3)>=(ang23-5)  && angle_dist2(cand2,cand3)<=(ang23+5)) && ...
                        (angle_dist2(cand2,cand4)>=(ang24-5)  && angle_dist2(cand2,cand4)<=(ang24+5)) && ...
                        (angle_dist2(cand3,cand4)>=(ang34-5)  && angle_dist2(cand3,cand4)<=(ang34+5)) && ...
                        (real(acosd(dot(normals_array2{cand1},normalZ)))>=(ang1Z-5)  && real(acosd(dot(normals_array2{cand1},normalZ)))<=(ang1Z+5)) && ...
                        (real(acosd(dot(normals_array2{cand2},normalZ)))>=(ang2Z-5)  && real(acosd(dot(normals_array2{cand2},normalZ)))<=(ang2Z+5)) && ...
                        (real(acosd(dot(normals_array2{cand3},normalZ)))>=(ang3Z-5)  && real(acosd(dot(normals_array2{cand3},normalZ)))<=(ang3Z+5)) && ...
                        (real(acosd(dot(normals_array2{cand4},normalZ)))>=(ang4Z-5)  && real(acosd(dot(normals_array2{cand4},normalZ)))<=(ang4Z+5))
                    time_angle = time_angle + toc(time1);
                    if vis_flag
                        figure(2);
                        hold off;
                        plot3(plane_array2{cand1}(:,1),plane_array2{cand1}(:,2),plane_array2{cand1}(:,3),'.b','MarkerSize',1);
                        hold on;
                        plot3(plane_array2{cand2}(:,1),plane_array2{cand2}(:,2),plane_array2{cand2}(:,3),'.r','MarkerSize',1);
                        plot3(plane_array2{cand3}(:,1),plane_array2{cand3}(:,2),plane_array2{cand3}(:,3),'.g','MarkerSize',1);
                        plot3(plane_array2{cand4}(:,1),plane_array2{cand4}(:,2),plane_array2{cand4}(:,3),'.y','MarkerSize',1);
                        grid on;
                    end
                    tested_bases = tested_bases + 1;
                    n11 = normals_array1{plane1_idx};
                    n21 = normals_array2{cand1};
                    n12 = normals_array1{plane2_idx};
                    n22 = normals_array2{cand2};
                    n13 = normals_array1{plane3_idx};
                    n23 = normals_array2{cand3};
                    n14 = normals_array1{plane4_idx};
                    n24 = normals_array2{cand4};
                    
                    [~,~,transform] = procrustes([n11' ; n12' ; n13' ; n14' ; normalZ'], [n21' ; n22' ; n23' ; n24' ; normalZ'],'scaling',0,'reflection',0);
                    
                    R = transform.T;
                    
                    if R(3,3) >=0.9
                        
                        n1 = normals_array1{plane1_idx};
                        n2 = normals_array1{plane2_idx};
                        n3 = normals_array1{plane3_idx};
                        P1_aux = points_array1{plane1_idx};
                        P2_aux = points_array1{plane2_idx};
                        P3_aux = points_array1{plane3_idx};
                        A = [n1(1) n1(2) n1(3) ; n2(1) n2(2) n2(3) ; n3(1) n3(2) n3(3)];
                        B = [P1_aux*n1 ; P2_aux*n2 ; P3_aux*n3];
                        P1 = A\B;
                        P1 = P1';
                        
                        n1 = normals_array2{cand1};
                        n2 = normals_array2{cand2};
                        n3 = normals_array2{cand3};
                        P1_aux = points_array2{cand1};
                        P2_aux = points_array2{cand2};
                        P3_aux = points_array2{cand3};
                        A = [n1(1) n1(2) n1(3) ; n2(1) n2(2) n2(3) ; n3(1) n3(2) n3(3)];
                        B = [P1_aux*n1 ; P2_aux*n2 ; P3_aux*n3];
                        P2 = A\B;
                        
                        P2a = P2' * R;
                        T = P1-P2'- (P2a-P2');
                        
                        time2 = tic;
                        centroids2_reg = centroids2 * R + repmat(T,size(centroids2,1),1);
                        normals2_reg = normals2 * R;
                        supp_plane_dist = zeros(size(plane_array2,2),5);
                        supp_plane_idx = zeros(size(plane_array2,2),5);
                        supp_plane_orient = zeros(size(plane_array2,2),5);
                        for yy = 1:size(plane_array2,2)
                            dist = zeros(size(plane_array1,2),1);
                            for xx = 1:size(plane_array1,2)
                                
                                A1 = normals_array1{xx}(1);
                                B1 = normals_array1{xx}(2);
                                C1 = normals_array1{xx}(3);
                                D1 = -points_array1{xx}*normals_array1{xx};
                                dist(xx) = abs(A1*centroids2_reg(yy,1) + B1*centroids2_reg(yy,2) + C1*centroids2_reg(yy,3) + D1)/sqrt(A1*A1 + B1*B1 + C1*C1);
                                
                            end
%                             [aux_dist, aux_idx] = sort(dist,'ascend');
%                             supp_plane_dist(yy,:) = aux_dist(1:5);
%                             supp_plane_idx(yy,:) = aux_idx(1:5);
                            orient = (normals1(:,1)-normals2_reg(yy,1)>-0.1 & normals1(:,2)-normals2_reg(yy,2)>-0.1 & normals1(:,3)-normals2_reg(yy,3)>-0.1) & (normals1(:,1)-normals2_reg(yy,1)< 0.1 & normals1(:,2)-normals2_reg(yy,2)< 0.1 & normals1(:,3)-normals2_reg(yy,3)< 0.1);
                            [aux_dist, aux_idx] = sort(dist,'ascend');
                            orient = orient(aux_idx);
                            supp_plane_dist(yy,:) = aux_dist(1:5);
                            supp_plane_idx(yy,:) = aux_idx(1:5);
                            supp_plane_orient(yy,:) = orient(1:5);
                        end
%                         n_supporting_planes = sum(supp_plane_dist(:,1)<=sup_thresh);
                        n_supporting_planes = sum(supp_plane_dist(:,1)<=sup_thresh & supp_plane_orient(:,1)==1);
                        time_plane = time_plane + toc(time2);
                        if  (n_supporting_planes*100/size(points_array2,2))>=20
                            
                            time3 = tic;
                            plane_sup_bases = plane_sup_bases +1;
                            centroids = [points_array2{cand1} ; points_array2{cand2} ; points_array2{cand3} ; points_array2{cand4}];
                            centroids_reg = centroids * R + repmat(T,size(centroids,1),1);
                            alfa = acosd(dot(normalZ,normals_array1{plane1_idx})/sqrt(normals_array1{plane1_idx}(1)*normals_array1{plane1_idx}(1) + normals_array1{plane1_idx}(2)*normals_array1{plane1_idx}(2) + normals_array1{plane1_idx}(3)*normals_array1{plane1_idx}(3)));
                            RR = [cosd(alfa) -sind(alfa) 0 ; sind(alfa) cosd(alfa) 0 ; 0 0 cosd(alfa)+(1-cosd(alfa))];
                            plane_aux = plane_array1{plane1_idx} * RR;
                            coord_idx = var(plane_aux) ~= min(var(plane_aux));
                            plane_aux1 = sortrows(plane_aux(:,coord_idx ));
                            borde = boundary(double(plane_aux1));
                            plane_aux = centroids_reg(1,:) * RR;
                            plane_aux2 = plane_aux(:,coord_idx );
                            in1 = inpolygon(plane_aux2(:,1),plane_aux2(:,2),plane_aux1(borde,1),plane_aux1(borde,2));
                            
                            alfa = acosd(dot(normalZ,normals_array1{plane2_idx})/sqrt(normals_array1{plane2_idx}(1)*normals_array1{plane2_idx}(1) + normals_array1{plane2_idx}(2)*normals_array1{plane2_idx}(2) + normals_array1{plane2_idx}(3)*normals_array1{plane2_idx}(3)));
                            RR = [cosd(alfa) -sind(alfa) 0 ; sind(alfa) cosd(alfa) 0 ; 0 0 cosd(alfa)+(1-cosd(alfa))];
                            plane_aux = plane_array1{plane2_idx} * RR;
                            coord_idx = var(plane_aux) ~= min(var(plane_aux));
                            plane_aux1 = sortrows(plane_aux(:,coord_idx));
                            borde = boundary(double(plane_aux1));
                            plane_aux = centroids_reg(2,:) * RR;
                            plane_aux2 = plane_aux(:,coord_idx);
                            in2 = inpolygon(plane_aux2(:,1),plane_aux2(:,2),plane_aux1(borde,1),plane_aux1(borde,2));
                            
                            alfa = acosd(dot(normalZ,normals_array1{plane3_idx})/sqrt(normals_array1{plane3_idx}(1)*normals_array1{plane3_idx}(1) + normals_array1{plane3_idx}(2)*normals_array1{plane3_idx}(2) + normals_array1{plane3_idx}(3)*normals_array1{plane3_idx}(3)));
                            RR = [cosd(alfa) -sind(alfa) 0 ; sind(alfa) cosd(alfa) 0 ; 0 0 cosd(alfa)+(1-cosd(alfa))];
                            plane_aux = plane_array1{plane3_idx} * RR;
                            coord_idx = var(plane_aux) ~= min(var(plane_aux));
                            plane_aux1 = sortrows(plane_aux(:,coord_idx));
                            borde = boundary(double(plane_aux1));
                            plane_aux = centroids_reg(3,:) * RR;
                            plane_aux2 = plane_aux(:,coord_idx);
                            in3 = inpolygon(plane_aux2(:,1),plane_aux2(:,2),plane_aux1(borde,1),plane_aux1(borde,2));
                            
                            alfa = acosd(dot(normalZ,normals_array1{plane4_idx})/sqrt(normals_array1{plane4_idx}(1)*normals_array1{plane4_idx}(1) + normals_array1{plane4_idx}(2)*normals_array1{plane4_idx}(2) + normals_array1{plane4_idx}(3)*normals_array1{plane4_idx}(3)));
                            RR = [cosd(alfa) -sind(alfa) 0 ; sind(alfa) cosd(alfa) 0 ; 0 0 cosd(alfa)+(1-cosd(alfa))];
                            plane_aux = plane_array1{plane4_idx} * RR;
                            coord_idx = var(plane_aux) ~= min(var(plane_aux));
                            plane_aux1 = sortrows(plane_aux(:,coord_idx));
                            borde = boundary(double(plane_aux1));
                            plane_aux = centroids_reg(4,:) * RR;
                            plane_aux2 = plane_aux(:,coord_idx);
                            in4 = inpolygon(plane_aux2(:,1),plane_aux2(:,2),plane_aux1(borde,1),plane_aux1(borde,2));
                            time_point = time_point + toc(time3);
                            
                            if sum([in1 in2 in3 in4]) == 4
                                
                                centroid_sup_bases = centroid_sup_bases + 1;
                                Tmax = [R' T' ; 0 0 0 1];
                                
                                time4 = tic;
                                plane_array_aux = plane_array2;
                                %                                 plane_array_aux([cand1 cand2 cand3 cand4]) = [];
                                %                                 supp_plane_idx([cand1 cand2 cand3 cand4]) = [];
                                %                                 supp_plane_dist([cand1 cand2 cand3 cand4]) = [];
                                supp_plane_idx = supp_plane_idx(supp_plane_dist(:,1)<=sup_thresh,:);
                                points2 = plane_array_aux(supp_plane_dist(:,1)<=sup_thresh);
                                dist = cell(1,size(supp_plane_idx,1));
                                acum_points_in = 0;
                                acum_points_tot = 0;
                                
                                for zz = 1 : size(supp_plane_idx,1)
                                   
                                    points_reg = points2{zz} * R + repmat(T,length(points2{zz}),1);
                                    in = zeros(5,size(points_reg,1));
                                    
                                    for tt = 1:5
                                        
                                        A1 = normals_array1{supp_plane_idx(zz,tt)}(1);
                                        B1 = normals_array1{supp_plane_idx(zz,tt)}(2);
                                        C1 = normals_array1{supp_plane_idx(zz,tt)}(3);
                                        alfa = acosd(dot(normalZ,normals_array1{supp_plane_idx(zz,tt)})/sqrt(A1*A1 + B1*B1 + C1*C1));
                                        RR = [cosd(alfa) -sind(alfa) 0 ; sind(alfa) cosd(alfa) 0 ; 0 0 cosd(alfa)+(1-cosd(alfa))];
                                        plane_aux = plane_array1{supp_plane_idx(zz,tt)} * RR;
                                        coord_idx = var(plane_aux) ~= min(var(plane_aux));
                                        plane_aux1 = sortrows(plane_aux(:,coord_idx));
                                        borde = boundary(double(plane_aux1));
                                        plane_aux = points_reg * RR;
                                        plane_aux2 = plane_aux(:,coord_idx);
                                        in(tt,:) = inpolygon(plane_aux2(:,1),plane_aux2(:,2),plane_aux1(borde,1),plane_aux1(borde,2));
                                        %                                         cont_in(tt) = sum(in);
                                        
                                    end
                                    
                                    [~, in_idx] = max(sum(in,2));
                                    A1 = normals_array1{supp_plane_idx(zz,1)}(1);
                                    B1 = normals_array1{supp_plane_idx(zz,1)}(2);
                                    C1 = normals_array1{supp_plane_idx(zz,1)}(3);
                                    D1 = -points_array1{supp_plane_idx(zz,1)}*normals_array1{supp_plane_idx(zz,1)};
                                    dist{zz} = abs(A1*points_reg(:,1) + B1*points_reg(:,2) + C1*points_reg(:,3) + D1)/sqrt(A1*A1 + B1*B1 + C1*C1);
                                    dist{zz}(in(in_idx,:)==false) = dist{zz}(in(in_idx,:)==false)*NaN;
                                    acum_points_in = acum_points_in + sum(in(in_idx,:));
                                    acum_points_tot = acum_points_tot + size(points_reg,1);
                                    
                                    %                                     A1 = normals_array1{supp_plane_idx(zz)}(1);
                                    %                                     B1 = normals_array1{supp_plane_idx(zz)}(2);
                                    %                                     C1 = normals_array1{supp_plane_idx(zz)}(3);
                                    %                                     D1 = -points_array1{supp_plane_idx(zz)}*normals_array1{supp_plane_idx(zz)};
                                    %                                     dist{zz} = abs(A1*points_reg(:,1) + B1*points_reg(:,2) + C1*points_reg(:,3) + D1)/sqrt(A1*A1 + B1*B1 + C1*C1);
                                    %
                                    %                                     alfa = acosd(dot(normalZ,normals_array1{supp_plane_idx(zz)})/sqrt(A1*A1 + B1*B1 + C1*C1));
                                    %                                     RR = [cosd(alfa) -sind(alfa) 0 ; sind(alfa) cosd(alfa) 0 ; 0 0 cosd(alfa)+(1-cosd(alfa))];
                                    %                                     plane_aux = plane_array1{supp_plane_idx(zz)} * RR;
                                    %                                     coord_idx = var(plane_aux) ~= min(var(plane_aux));
                                    %                                     plane_aux1 = sortrows(plane_aux(:,coord_idx));
                                    %                                     borde = boundary(double(plane_aux1));
                                    %                                     plane_aux = points_reg * RR;
                                    %                                     plane_aux2 = plane_aux(:,coord_idx);
                                    %                                     in = inpolygon(plane_aux2(:,1),plane_aux2(:,2),plane_aux1(borde,1),plane_aux1(borde,2));
                                    %                                     dist{zz}(in==false) = dist{zz}(in==false)*2;
                                    %                                     acum_points_in = acum_points_in + sum(in);
                                    %                                     acum_points_tot = acum_points_tot + size(points_reg,1);
                                    
                                end
                                dist = cell2mat(dist');
                                dist(isnan(dist)) = [];
                                rmse = mean(dist);
                                time_rmse = time_rmse + toc(time4);
                                points_in = acum_points_in*100/acum_points_tot;
                                
                                
                                %                                 rand_index = randi(length(plane_array2{cand1}),100,1);
                                %                                 plane1_reg = plane_array2{cand1}(rand_index,:) * R + repmat(T,length(plane_array2{cand1}(rand_index,:)),1);
                                %                                 A1 = normals_array1{plane1_idx}(1);
                                %                                 B1 = normals_array1{plane1_idx}(2);
                                %                                 C1 = normals_array1{plane1_idx}(3);
                                %                                 D1 = -points_array1{plane1_idx}*normals_array1{plane1_idx};
                                %                                 dist1 = abs(A1*plane1_reg(:,1) + B1*plane1_reg(:,2) + C1*plane1_reg(:,3) + D1)/sqrt(A1*A1 + B1*B1 + C1*C1);
                                %
                                %                                 rand_index = randi(length(plane_array2{cand2}),100,1);
                                %                                 plane2_reg = plane_array2{cand2}(rand_index,:) * R + repmat(T,length(plane_array2{cand2}(rand_index,:)),1);
                                %                                 A1 = normals_array1{plane2_idx}(1);
                                %                                 B1 = normals_array1{plane2_idx}(2);
                                %                                 C1 = normals_array1{plane2_idx}(3);
                                %                                 D1 = -points_array1{plane2_idx}*normals_array1{plane2_idx};
                                %                                 dist2 = abs(A1*plane2_reg(:,1) + B1*plane2_reg(:,2) + C1*plane2_reg(:,3) + D1)/sqrt(A1*A1 + B1*B1 + C1*C1);
                                %
                                %                                 rand_index = randi(length(plane_array2{cand3}),100,1);
                                %                                 plane3_reg = plane_array2{cand3}(rand_index,:) * R + repmat(T,length(plane_array2{cand3}(rand_index,:)),1);
                                %                                 A1 = normals_array1{plane3_idx}(1);
                                %                                 B1 = normals_array1{plane3_idx}(2);
                                %                                 C1 = normals_array1{plane3_idx}(3);
                                %                                 D1 = -points_array1{plane3_idx}*normals_array1{plane3_idx};
                                %                                 dist3 = abs(A1*plane3_reg(:,1) + B1*plane3_reg(:,2) + C1*plane3_reg(:,3) + D1)/sqrt(A1*A1 + B1*B1 + C1*C1);
                                %
                                %                                 rand_index = randi(length(plane_array2{cand4}),100,1);
                                %                                 plane4_reg = plane_array2{cand4}(rand_index,:) * R + repmat(T,length(plane_array2{cand4}(rand_index,:)),1);
                                %                                 A1 = normals_array1{plane4_idx}(1);
                                %                                 B1 = normals_array1{plane4_idx}(2);
                                %                                 C1 = normals_array1{plane4_idx}(3);
                                %                                 D1 = -points_array1{plane4_idx}*normals_array1{plane4_idx};
                                %                                 dist4 = abs(A1*plane4_reg(:,1) + B1*plane4_reg(:,2) + C1*plane4_reg(:,3) + D1)/sqrt(A1*A1 + B1*B1 + C1*C1);
                                %                                 rmse = mean([dist1 ; dist2 ; dist3 ; dist4]);
                                
                                trans(cont) = struct('Tmax', Tmax, 'Supporting_planes', n_supporting_planes, 'Supporting_points', points_in, 'RMSE', rmse, 'Corresponding_planes', [plane1_idx plane2_idx plane3_idx plane4_idx cand1 cand2 cand3 cand4]);
                                cont = cont + 1;
                                disp(['Cont : ',num2str(cont)]);
                            end
                            
                        end
                    end
                    
                end
                %                 toalign_congruent_sets = [toalign_congruent_sets ; cand_congruent_sets];
                
            end
           
            
        end
        
    end
end


disp(['Number of candidates pairs of planes: ', num2str(cant_cand)]);
disp(['Congruent sets: ', num2str(tested_bases)]);
disp(['Congruent sets prime: ', num2str(plane_sup_bases)]);
disp(['Number of transformations: ', num2str(centroid_sup_bases)]);
disp(['Time spent in angle verification: ', num2str(time_angle)]);
disp(['Time spent in plane support: ', num2str(time_plane)]);
disp(['Time spent in point support: ', num2str(time_point)]);
disp(['Time spent in RMSE: ', num2str(time_rmse)]);
end