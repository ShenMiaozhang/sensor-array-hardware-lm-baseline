function run_baseline_solver()
% Baseline solver (magnetic dipole model + LM optimization).
% One-click usage:
%   1) Open MATLAB in this Software folder
%   2) Run: run_baseline_solver

cfg = defaultConfig();
[sensorPos, sensorN] = buildSensorLayout(cfg);

init = cfg.init;
result = init;
elapsedTimeS = NaN;
drawTrace = false;
trace = zeros(cfg.maxLoopTimes, 6);
traceCnt = 0;

fig = figure('Name', 'Baseline LM Solver', 'Color', 'w');
set(fig, 'CurrentCharacter', char(0));

ser = serialport(cfg.comPort, cfg.baudRate);
cleanupObj = onCleanup(@() cleanupSerial(ser));

disp('Collecting compensation frame...');
blCompensate = readSensorFrame(ser, sensorN, cfg.serialMaxRetry);

for cnt = 1:cfg.maxLoopTimes
    if ~ishandle(fig)
        break
    end

    bl = readSensorFrame(ser, sensorN, cfg.serialMaxRetry);

    % Geomagnetic compensation (Z axis is reversed w.r.t. X/Y in this setup).
    bl(:, 1) = -bl(:, 1) + blCompensate(:, 1);
    bl(:, 2) =  bl(:, 2) - blCompensate(:, 2);
    bl(:, 3) =  bl(:, 3) - blCompensate(:, 3);

    magStrength = vecnorm(bl, 2, 2);
    strongIdx = magStrength >= cfg.magThreshold;
    strongSensorPos = sensorPos(strongIdx, :);
    strongBl = bl(strongIdx, :);
    strongSensorCount = nnz(strongIdx);

    if strongSensorCount >= cfg.minStrongSensorCount
        tStart = tic;
        result = lsqnonlin( ...
            @(param) solverResidual(param, strongSensorPos, strongBl, cfg.NT), ...
            init, cfg.lb, cfg.ub, [], [], [], [], @unitVectorConstraint, cfg.lsqOptions);
        elapsedTimeS = toc(tStart);
        init = result;

        [alphaDeg, betaDeg, gammaDeg] = directionAnglesDeg(result(1:3));
        fprintf('solve_time = %.6f s\n', elapsedTimeS);
        fprintf('x: %.4f mm\ty: %.4f mm\tz: %.4f mm\tm: %.4f\tn: %.4f\tp: %.4f\tcnt: %d\n', ...
            result(4) * 1000, result(5) * 1000, result(6) * 1000, ...
            result(1), result(2), result(3), cnt);
        fprintf('alpha: %.1f deg\tbeta: %.1f deg\tgamma: %.1f deg\n', alphaDeg, betaDeg, gammaDeg);
    else
        fprintf('Not enough valid sensor data to perform optimization.\n');
    end

    renderFrame(sensorPos, bl, result, drawTrace, trace, traceCnt);

    currentKey = lower(get(fig, 'CurrentCharacter'));
    set(fig, 'CurrentCharacter', char(0));

    switch currentKey
        case 'q'
            close(fig);
            break
        case 't'
            drawTrace = true;
        case 'c'
            drawTrace = false;
            traceCnt = 0;
        case 'r'
            [alphaDeg, betaDeg, gammaDeg] = directionAnglesDeg(result(1:3));
            appendResultRow(cfg.outputFile, result, alphaDeg, betaDeg, gammaDeg, elapsedTimeS);
            disp('Magnet data has been saved to the Excel file.');
    end

    if drawTrace
        traceCnt = traceCnt + 1;
        if traceCnt <= size(trace, 1)
            trace(traceCnt, :) = result;
        end
    end

    drawnow limitrate;
end

clear cleanupObj;
end

function cfg = defaultConfig()
cfg.comPort = 'COM3';
cfg.baudRate = 115200;
cfg.sensorN = 5;
cfg.sensorIntervalX = 40.14 / 1000;
cfg.sensorIntervalY = 41.94 / 1000;
cfg.magThreshold = 1;
cfg.minStrongSensorCount = 3;
cfg.maxLoopTimes = 10000;
cfg.serialMaxRetry = 10;

% Magnetic moment magnitude (N_T). Default for baseline magnet is 5.8.
% Different magnets require different N_T; update this value before running.
cfg.NT = 5.8;

cfg.init = [0, 0, 1, 0, 0, 0.1];
cfg.lb = [-1, -1, -1, -0.1, -0.1, -0.1];
cfg.ub = [1, 1, 1, 0.15, 0.15, 0.15];
cfg.outputFile = 'magnet_data.xlsx';

cfg.lsqOptions = optimoptions('lsqnonlin', ...
    'Display', 'off', ...
    'FunctionTolerance', 1e-10, ...
    'StepTolerance', 1e-10);
end

function [sensorPos, sensorN] = buildSensorLayout(cfg)
sensorN = cfg.sensorN^2;
sensorPos = zeros(sensorN, 3);
for i = 1:cfg.sensorN
    for j = 1:cfg.sensorN
        idx = (i - 1) * cfg.sensorN + j;
        sensorPos(idx, 1) = ((j - 1) - (cfg.sensorN - 1) / 2) * cfg.sensorIntervalX;
        sensorPos(idx, 2) = ((i - 1) - (cfg.sensorN - 1) / 2) * cfg.sensorIntervalY;
        sensorPos(idx, 3) = 0;
    end
end
end

function bl = readSensorFrame(ser, N, maxRetry)
% Serial frame protocol:
%   Line 1: "U<space-separated N values>"
%   Line 2: "V<space-separated N values>"
%   Line 3: "W<space-separated N values>"
% The three lines are parsed into an N-by-3 matrix:
%   column 1 -> U axis, column 2 -> V axis, column 3 -> W axis.
bl = [];
errCnt = 0;

while isempty(bl)
    flush(ser);
    lineU = strtrim(readline(ser));
    lineV = strtrim(readline(ser));
    lineW = strtrim(readline(ser));

    U = parseAxisPayload(lineU, 'U');
    V = parseAxisPayload(lineV, 'V');
    W = parseAxisPayload(lineW, 'W');

    if numel(U) == N && numel(V) == N && numel(W) == N
        bl = [U(:), V(:), W(:)];
        break
    end

    errCnt = errCnt + 1;
    if errCnt >= maxRetry
        error('Fail to parse serial frame too many times.');
    end
end
end

function values = parseAxisPayload(lineText, axisPrefix)
values = [];
if startsWith(lineText, axisPrefix)
    values = sscanf(extractAfter(lineText, 1), '%f');
end
end

function F = solverResidual(param, sensorPos, bl, NT)
m = param(1);
n = param(2);
p = param(3);
a = param(4);
b = param(5);
c = param(6);

sensorCount = size(sensorPos, 1);
F = zeros(sensorCount * 3, 1);

for i = 1:sensorCount
    [bx, by, bz] = magneticFluxDensityComponent( ...
        NT, m, n, p, a, b, c, sensorPos(i, 1), sensorPos(i, 2), sensorPos(i, 3), 0);
    idx = (i - 1) * 3 + 1;
    F(idx:idx + 2) = [bx - bl(i, 1); by - bl(i, 2); bz - bl(i, 3)];
end
end

function [c, ceq] = unitVectorConstraint(param)
c = [];
ceq = param(1)^2 + param(2)^2 + param(3)^2 - 1;
end

function [blx, bly, blz] = magneticFluxDensityComponent(NT, m, n, p, a, b, c, x, y, z, noise)
dx = x - a;
dy = y - b;
dz = z - c;
Rl = sqrt(dx * dx + dy * dy + dz * dz);
if Rl < eps
    Rl = eps;
end

comm = 3 * (m * dx + n * dy + p * dz);
blx = NT * ((comm * dx) / (Rl^5) - m / (Rl^3)) + noise * randn;
bly = NT * ((comm * dy) / (Rl^5) - n / (Rl^3)) + noise * randn;
blz = NT * ((comm * dz) / (Rl^5) - p / (Rl^3)) + noise * randn;
end

function [alphaDeg, betaDeg, gammaDeg] = directionAnglesDeg(directionVec)
directionVec = max(min(directionVec(:), 1), -1);
alphaDeg = rad2deg(acos(directionVec(1)));
betaDeg = rad2deg(acos(directionVec(2)));
gammaDeg = rad2deg(acos(directionVec(3)));
end

function renderFrame(sensorPos, bl, result, drawTrace, trace, traceCnt)
quiver3(sensorPos(:, 1), sensorPos(:, 2), sensorPos(:, 3), ...
    bl(:, 1), bl(:, 2), bl(:, 3), '-b', 'LineWidth', 1.2);
hold on;
scatter3(sensorPos(:, 1), sensorPos(:, 2), sensorPos(:, 3), 12, 'r', 'filled');

xlim([-0.1, 0.1]);
ylim([-0.1, 0.1]);
zlim([-0.1, 0.2]);
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
set(gca, 'YDir', 'reverse');
grid on;

if drawTrace && traceCnt > 1
    validCnt = min(traceCnt, size(trace, 1));
    plot3(trace(1:validCnt, 4), trace(1:validCnt, 5), trace(1:validCnt, 6), '-r', 'LineWidth', 1);
    scatter3(trace(1:validCnt, 4), trace(1:validCnt, 5), trace(1:validCnt, 6), 18, 'b', 'filled');
end

scatter3(result(4), result(5), result(6), 80, 'g', 'filled');
quiver3(result(4), result(5), result(6), result(1), result(2), result(3), ...
    0.02, 'Color', 'm', 'LineWidth', 2, 'MaxHeadSize', 1);
hold off;
end

function appendResultRow(outputFile, result, alphaDeg, betaDeg, gammaDeg, solveTimeS)
if exist(outputFile, 'file') == 0
    header = {'X (mm)', 'Y (mm)', 'Z (mm)', 'Alpha (deg)', 'Beta (deg)', 'Gamma (deg)', 'Solve Time (s)'};
    writecell(header, outputFile);
end

xMm = result(4) * 1000;
yMm = result(5) * 1000;
zMm = result(6) * 1000;
row = {xMm, yMm, zMm, alphaDeg, betaDeg, gammaDeg, solveTimeS};
writecell(row, outputFile, 'WriteMode', 'append');
end

function cleanupSerial(ser)
if isa(ser, 'serialport')
    delete(ser);
end
end
