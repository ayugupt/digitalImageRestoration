image = imread('banana_holes2', 'png');

[m, n, q] = size(image);

hole = [0, 255, 0];

holes = {};
noHoles = 0;


for i = 2:m-1
    for j = 2:n-1
       temp1 = image(i, j, :);
       temp2 = image(i, j-1, :);
       temp3 = image(i-1, j, :);
       temp4 = image(i, j+1, :);
       temp5 = image(i+1, j, :);
        if isequal(hole, temp1(1, :))
            if ~isequal(hole, temp2(1, :)) && ~isequal(hole, temp3(1, :))
                noHoles = noHoles + 1;
                holes{noHoles}(1, 1) = i;
                holes{noHoles}(1, 2) = j;
            elseif ~isequal(hole, temp4(1, :)) && ~isequal(hole, temp3(1, :))
                holes{noHoles}(2, 1) = i;
                holes{noHoles}(2, 2) = j;
            elseif ~isequal(hole, temp2(1, :)) && ~isequal(hole, temp5(1, :))
                holes{noHoles}(3, 1) = i;
                holes{noHoles}(3, 2) = j;
            elseif ~isequal(hole, temp4(1, :)) && ~isequal(hole, temp5(1, :))
                holes{noHoles}(4, 1) = i;
                holes{noHoles}(4, 2) = j;
            end
        end
    end
end

% pos = {};
% for o = 1:noHoles
%     index = 1;
%     for rt = holes{o}(1, 1):holes{o}(4, 1)
%         for ct = holes{o}(1, 2):holes{o}(4, 2)
%             pos{o}(index, 1) = rt;
%             pos{o}(index, 2) = ct;
%             index = index + 1;
%         end
%     end
%     
% end


for hol = 1:noHoles
h1 = [3, 4, holes{hol}(1, 2), holes{hol}(2, 2), holes{hol}(4, 2), holes{hol}(3, 2) , holes{hol}(1, 1), holes{hol}(2, 1), holes{hol}(4, 1),holes{hol}(3, 1)]';

g = decsg(h1);
model1 = createpde;
model2 = createpde;
model3 = createpde;

geometryFromEdges(model1, g);
geometryFromEdges(model2, g);
geometryFromEdges(model3, g);
% pdegplot(model1, 'EdgeLabels', 'on')
% axis equal
% xlim([0 1024])
% ylim([0 512])

applyBoundaryCondition(model1, 'dirichlet', 'Edge', 1, 'u',image(holes{hol}(1, 1)-1, holes{hol}(1, 2), 1));
applyBoundaryCondition(model1, 'dirichlet', 'Edge', 2, 'u',image(holes{hol}(2, 1), holes{hol}(2, 2)+1, 1));
applyBoundaryCondition(model1, 'dirichlet', 'Edge', 3, 'u',image(holes{hol}(4, 1)+1, holes{hol}(4, 2), 1));
applyBoundaryCondition(model1, 'dirichlet', 'Edge', 4, 'u',image(holes{hol}(3, 1), holes{hol}(3, 2)-1, 1));

applyBoundaryCondition(model2, 'dirichlet', 'Edge', 1, 'u',image(holes{hol}(1, 1)-1, holes{hol}(1, 2), 2));
applyBoundaryCondition(model2, 'dirichlet', 'Edge', 2, 'u',image(holes{hol}(2, 1), holes{hol}(2, 2)+1, 2));
applyBoundaryCondition(model2, 'dirichlet', 'Edge', 3, 'u',image(holes{hol}(4, 1)+1, holes{hol}(4, 2), 2));
applyBoundaryCondition(model2, 'dirichlet', 'Edge', 4, 'u',image(holes{hol}(3, 1), holes{hol}(3, 2)-1, 2));

applyBoundaryCondition(model3, 'dirichlet', 'Edge', 1, 'u',image(holes{hol}(1, 1)-1, holes{hol}(1, 2), 3));
applyBoundaryCondition(model3, 'dirichlet', 'Edge', 2, 'u',image(holes{hol}(2, 1), holes{hol}(2, 2)+1, 3));
applyBoundaryCondition(model3, 'dirichlet', 'Edge', 3, 'u',image(holes{hol}(4, 1)+1, holes{hol}(4, 2), 3));
applyBoundaryCondition(model3, 'dirichlet', 'Edge', 4, 'u',image(holes{hol}(3, 1), holes{hol}(3, 2)-1, 3));



specifyCoefficients(model1,'m',0,'d',0,'c',1,'a',0,'f',0);
generateMesh(model1, 'Hmax', 2);
results1 = solvepde(model1); 
%pdeplot(model1,'XYData',results1.NodalSolution);

specifyCoefficients(model2,'m',0,'d',0,'c',1,'a',0,'f',0);
generateMesh(model2, 'Hmax', 2);
results2 = solvepde(model2); 
%pdeplot(model2,'XYData',results2.NodalSolution);

specifyCoefficients(model3,'m',0,'d',0,'c',1,'a',0,'f',0);
generateMesh(model3, 'Hmax', 2);
results3 = solvepde(model3); 
%pdeplot(model3,'XYData',results3.NodalSolution);

[X, Y] = meshgrid(holes{hol}(1, 2):holes{hol}(4, 2), holes{hol}(1, 1):holes{hol}(4, 1));

vals1 = interpolateSolution(results1, X, Y);
vals2 = interpolateSolution(results2, X, Y);
vals3 = interpolateSolution(results3, X, Y);

v1 = reshape(vals1, holes{hol}(4, 1) - holes{hol}(1, 1) + 1, holes{hol}(4, 2) - holes{hol}(1, 2) + 1, []);
v2 = reshape(vals2, holes{hol}(4, 1) - holes{hol}(1, 1) + 1, holes{hol}(4, 2) - holes{hol}(1, 2) + 1, []);
v3 = reshape(vals3, holes{hol}(4, 1) - holes{hol}(1, 1) + 1, holes{hol}(4, 2) - holes{hol}(1, 2) + 1, []);

for ri = holes{hol}(1, 1):holes{hol}(4, 1)
    for cj = holes{hol}(1, 2):holes{hol}(4, 2)
        image(ri, cj, 1) = v1(ri - holes{hol}(1, 1) + 1, cj - holes{hol}(1, 2) + 1);
        image(ri, cj, 2) = v2(ri - holes{hol}(1, 1) + 1, cj - holes{hol}(1, 2) + 1);
        image(ri, cj, 3) = v3(ri - holes{hol}(1, 1) + 1, cj - holes{hol}(1, 2) + 1);
    end
end
end

imshow(image)