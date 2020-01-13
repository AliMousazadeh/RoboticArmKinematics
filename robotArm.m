clear; close all; clc;

numArms = 5;
length = 10;
defAngle = pi;

lengths = length * ones(1, numArms);
angles = defAngle * ones(1, numArms);

goal = [30; -10];

dt = 0.01;
timeLimit = 20;
speed = 5;
distThresh = 0.1;
FrameRatePerSecond = 60;
%%%%%%%%%%

goal = [goal; 1];
numRuns = timeLimit * FrameRatePerSecond;
n = size(lengths, 2);
figure;
H = scatter(0, 0);
hold on;
scatter(goal(1), goal(2), 'g', 'filled');
distVals = zeros(numRuns, 1);


writerObj = VideoWriter('moving_arm13', 'MPEG-4');
writerObj.FrameRate = FrameRatePerSecond;
open(writerObj);

levVals = zeros(numRuns, 1);
for i = 1:numRuns
    X_w = forwardKinematic(angles, lengths);
    J = forwardKinematicDerv(angles, lengths);
    invJ = (J'*J)\J';
    
    dX = zeros(3, n);
    dX(:, end) = goal - X_w(:, end);
    
    dTheta = invJ * dX(:);
    dTheta = dTheta / (sqrt(dTheta'*dTheta));
    
    angles = angles + speed * dt * dTheta';
    
    H = scatter(X_w(1, :), X_w(2, :), 'r');
    scatter(0, 0, 'r', 'filled');
    plot(sum(lengths), sum(lengths), 'r');
    plot(-sum(lengths), -sum(lengths), 'r');
    
    plots(1) = plot([0, X_w(1, 1)], [0, X_w(2, 1)], 'r', 'LineWidth', 3);
    for j = 2:n
        plots(j) = plot([X_w(1, j-1), X_w(1, j)], ...
            [X_w(2, j-1), X_w(2, j)], 'r', 'LineWidth', 3);
    end
    
    levVals(i) = sqrt((X_w(1,1)-X_w(1,2))^2 + (X_w(2,1)-X_w(2,2))^2);
    
    hold on;
    scatter(X_w(1, end), X_w(2, end), 1, 'b', 'filled');
    scatter(goal(1), goal(2), 'g', 'filled');
    
    d = sqrt((goal - X_w(:, end))'*(goal - X_w(:, end)));
    distVals(i) = d;
    if (d < distThresh)
        fprintf('Distance reached within tolerance.\n');
        break;
    end
    
    drawnow
    F = getframe(gcf) ;
    writeVideo(writerObj, F);
    
    delete(H);
    for j = 1:n
        delete(plots(j));
    end
end

close(writerObj);
fprintf('Sucessfully generated the video.\n')

%%
figure;
plot((1:i), distVals(1:i));
legend('distance from goal')
xlabel('steps');








function X_w = forwardKinematic(angles, lengths)
n = size(lengths, 2);
X_local = [lengths; zeros(1, n)];
X_w = zeros(3, n);

Dw_ord = [
    0 1 0
    1 0 0
    0 0 1
    ];

Dw = [[rotateMat(angles(1)); zeros(1, 2)], [zeros(2, 1); 1]];
X_w(:, 1) = Dw * [X_local(:, 1); 1];
for i = 2:n
    Ds = [[rotateMat(angles(i)); zeros(1, 2)], [X_local(:, i-1); 1]];
    Dw = Dw * Ds;
    X_w(:, i) = Dw * [X_local(:, i); 1];
end
end

function result = forwardKinematicDerv(angles, lengths)
d = 1e-3;
n = size(lengths, 2);
ds = cell(n, 1);
result = [];
for i = 1:n
    msk = zeros(1, n);
    msk(i) = d;
    xp = forwardKinematic(angles+msk, lengths);
    xn = forwardKinematic(angles-msk, lengths);
    ds{i} = (xp-xn)/(2*d);
    
    result = [result, ds{i}(:)];
end
end

function result = rotateMat(theta)
result = [
    cos(theta), -sin(theta)
    sin(theta), cos(theta)
    ];
end