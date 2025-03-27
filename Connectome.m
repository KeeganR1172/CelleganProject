%% Load data
clc, clear, close all

load(['Connectome', filesep, 'GHermChem.mat'])
load(['Connectome', filesep, 'GHermElec_sym.mat'])

%%
t = tiledlayout(1,2);
fig1 = nexttile(t, 1);
fig1 = HermChemNodeConns(["dBWML10", "dBWMR10", "vBWML10", "vBWMR10"],1);

%%
fig2 = nexttile(t, 2);
fig2 = HermElecNodeConns(["dBWML10", "dBWMR10", "vBWML10", "vBWMR10"],1);

%% Maps all predecessors of Target

function nodefig = HermElecNodeConns(StartNode, depth)
    load(['Connectome', filesep, 'GHermElec_sym.mat'])
    
    PredNode = StartNode;
    NewPredNode = [];
    FiltNode = StartNode;
    FiltNeuron = StartNode;
    d_iterator = 0;
    
    while(not(isempty(PredNode)) && d_iterator < depth)
        numNodes = size(PredNode,2);
        d_iterator = d_iterator + 1;
        for i = 1:numNodes
            NewPredNode = [NewPredNode ; neighbors(GHermElec_Sym, PredNode(i))];
        end
        NewPredNode = unique(NewPredNode);
        PredNode = [];
        for i = 1:size(NewPredNode,1)
            if(sum(matches(FiltNode, NewPredNode(i))) == 0)
                PredNode = [PredNode, NewPredNode(i)];
            end
        end
        FiltNode = [FiltNode, PredNode];
    end
    
    for i = 1:size(FiltNode,2)
        if(sum(matches([InterNeurons, MotorNeurons, SensoryNeurons], FiltNode(i))) > 0)
            FiltNeuron = unique([FiltNeuron, FiltNode(i)]);
        end
    end
    GHermMusc = subgraph(GHermElec_Sym,FiltNode);
    GHermOnlyNeuron = subgraph(GHermMusc, FiltNeuron);
    
    GHermPlot = GHermMusc;
    
    NodeColor = [];
    for i = 1:size(GHermPlot.Nodes,1)
        if(matches(GHermPlot.Nodes.Group(i), 'BodyWallMuscles'))
            NodeColor(i,:) = [0 1 0];
        elseif(matches(GHermPlot.Nodes.Group(i), 'InterNeurons'))
            NodeColor(i,:) = [1 1 0];
        elseif(matches(GHermPlot.Nodes.Group(i), 'MotorNeurons'))
            NodeColor(i,:) = [0 0 1];
        elseif(matches(GHermPlot.Nodes.Group(i), 'SensoryNeurons'))
            NodeColor(i,:) = [1 0 0];
        end
    end
    LWidths = 5 * GHermPlot.Edges.Weight / max(GHermElec_Sym.Edges.Weight);
    nodefig = plot(GHermPlot, 'NodeColor', NodeColor, 'LineWidth', LWidths);
end

function nodefig = HermChemNodeConns(StartNode, depth)
    load(['Connectome', filesep, 'GHermChem.mat'])

    PredNode = StartNode;
    NewPredNode = [];
    FiltNode = StartNode;
    FiltNeuron = StartNode;
    d_iterator = 0;
    
    while(not(isempty(PredNode)) && d_iterator < depth)
        numNodes = size(PredNode,2);
        d_iterator = d_iterator + 1;
        for i = 1:numNodes
            NewPredNode = [NewPredNode ; predecessors(GHermChem, PredNode(i))];
        end
        NewPredNode = unique(NewPredNode);
        PredNode = [];
        for i = 1:size(NewPredNode,1)
            if(sum(matches(FiltNode, NewPredNode(i))) == 0)
                PredNode = [PredNode, NewPredNode(i)];
            end
        end
        FiltNode = [FiltNode, PredNode];
    end
    
    for i = 1:size(FiltNode,2)
        if(sum(matches([InterNeurons, MotorNeurons, SensoryNeurons], FiltNode(i))) > 0)
            FiltNeuron = unique([FiltNeuron, FiltNode(i)]);
        end
    end
    GHermMuscChem = subgraph(GHermChem,FiltNode);
    GHermOnlyNeuron = subgraph(GHermMuscChem, FiltNeuron);
    
    GHermPlot = GHermOnlyNeuron;
    
    NodeColor = [];
    for i = 1:size(GHermPlot.Nodes,1)
        if(matches(GHermPlot.Nodes.Group(i), 'BodyWallMuscles'))
            NodeColor(i,:) = [0 1 0];
        elseif(matches(GHermPlot.Nodes.Group(i), 'InterNeurons'))
            NodeColor(i,:) = [1 1 0];
        elseif(matches(GHermPlot.Nodes.Group(i), 'MotorNeurons'))
            NodeColor(i,:) = [0 0 1];
        elseif(matches(GHermPlot.Nodes.Group(i), 'SensoryNeurons'))
            NodeColor(i,:) = [1 0 0];
        end
    end
    LWidths = 5 * GHermPlot.Edges.Weight / max(GHermChem.Edges.Weight);
    nodefig = plot(GHermPlot, 'NodeColor', NodeColor, 'LineWidth', LWidths);
end

%% Test successors of BWMs
test = 0;
for i = 1:95
    test = [test; successors(GHermChem, BodyWallMuscles(i))];
end


% === Step 1: Build Motor Neuron → Muscle Map ===
motorToMuscleMap = containers.Map();  % key: motor neuron, value: list of muscles

for i = 1:length(BodyWallMuscles)
    muscle = BodyWallMuscles{i};
    
    % Get upstream (predecessor) nodes to this muscle
    preds = predecessors(GHermChem, muscle);
    
    % Filter for motor neurons only
    motorInputs = preds(ismember(preds, MotorNeurons));
    
    for j = 1:length(motorInputs)
        mNeuron = motorInputs{j};
        if isKey(motorToMuscleMap, mNeuron)
            motorToMuscleMap(mNeuron) = [motorToMuscleMap(mNeuron), {muscle}];
        else
            motorToMuscleMap(mNeuron) = {muscle};
        end
    end
end

% === Step 2: Display Preview of Mappings ===
disp("Motor neurons connected to muscles:");
motorKeys = keys(motorToMuscleMap);
for i = 1:min(5, length(motorKeys))
    fprintf("%s → %s\n", motorKeys{i}, strjoin(motorToMuscleMap(motorKeys{i}), ', '));
end

% === Step 3: Assign Fake Activation Values to Motor Neurons ===
motorActivation = containers.Map();
for i = 1:length(motorKeys)
    motorActivation(motorKeys{i}) = rand();  % Random value between 0–1
end

% === Step 4: Create Subgraph of Motor Neurons + Muscles ===
muscleLists = values(motorToMuscleMap);
flatMuscleList = [muscleLists{:}];  % flatten nested cells
allNodes = unique([motorKeys, flatMuscleList]);

Gsub = subgraph(GHermChem, allNodes);

% === Step 5: Assign Node Colors Based on Activation ===
colorMap = zeros(length(Gsub.Nodes.Name), 3);  % RGB

for i = 1:length(Gsub.Nodes.Name)
    n = Gsub.Nodes.Name{i};
    
    if ismember(n, MotorNeurons)
        val = motorActivation(n);  % 0–1
        % Blue (low) to white (high)
        colorMap(i,:) = [1 - val, 1 - val, 1];
        
    elseif ismember(n, BodyWallMuscles)
        % Muscles → green
        colorMap(i,:) = [0.3, 1, 0.3];
        
    else
        % Default (if other types sneak in)
        colorMap(i,:) = [0.8, 0.8, 0.8];
    end
end

% === Step 6: Visualize the Graph ===
figure;
p = plot(Gsub, ...
    'Layout', 'layered', ...
    'Direction', 'right', ...
    'NodeColor', colorMap, ...
    'MarkerSize', 6, ...
    'EdgeAlpha', 0.4);

title('Motor Neuron → Muscle Map with Activation Levels');

%%


% Make sure this variable exists
if ~exist('motorToMuscleMap', 'var')
    error('Please run generateMotorNeuronToMuscleMap.m first to generate motorToMuscleMap.');
end

% === Step 1: Choose 5 Motor Neurons (spaced along body) ===
selectedMotorNeurons = {'AS01', 'AS02', 'AS03', 'AS04', 'AS05'};
% Check that all are in your map
missing = ~isKey(motorToMuscleMap, selectedMotorNeurons);
if any(missing)
    error('Some selected motor neurons not found in motorToMuscleMap. Try different ones.');
end

% Get corresponding muscles
selectedMuscles = cell(1, length(selectedMotorNeurons));
for i = 1:length(selectedMotorNeurons)
    muscles = motorToMuscleMap(selectedMotorNeurons{i});
    selectedMuscles{i} = muscles{1};  % Just use first muscle mapped for each motor neuron
end

fprintf("Selected Neuron → Muscle Mappings:\n");
for i = 1:length(selectedMotorNeurons)
    fprintf("  %s → %s\n", selectedMotorNeurons{i}, selectedMuscles{i});
end

% === Step 2: Simulate Traveling Activation Wave ===
T = 100;  % time steps
numSegments = length(selectedMotorNeurons);
activationMatrix = zeros(numSegments, T);

for t = 1:T
    for seg = 1:numSegments
        phaseShift = (seg-1) * pi/4;  % stagger each segment
        activationMatrix(seg, t) = 0.5 + 0.5 * sin(2 * pi * 0.05 * t - phaseShift);
    end
end

% === Step 3: Animate Muscle Activation as Bars ===
figure('Color', 'w');
x = 1:numSegments;
h = bar(x, activationMatrix(:,1), 'FaceColor', 'flat');
ylim([0 1]); grid on;
xlabel('Muscle Segment'); ylabel('Activation');
title('Worm-like Muscle Activation');
xticklabels(selectedMuscles);

% Animation loop
for t = 1:T
    h.YData = activationMatrix(:,t);
    h.CData = [zeros(numSegments,1), activationMatrix(:,t), ones(numSegments,1)];  % blue → cyan
    title(['Time Step: ', num2str(t)]);
    pause(0.03);
    drawnow;
end

%%


% ==== Parameters ====
numSegments = 5;           % Number of worm body segments
segmentLength = 1;         % Fixed length (no contraction)
T = 100;                   % Time steps
waveSpeed = 0.05;          % Speed of wave across time
phaseOffset = pi/4;        % Phase shift per segment

% ==== Generate activation wave (like bar graph) ====
activationMatrix = zeros(numSegments, T);
for t = 1:T
    for seg = 1:numSegments
        activationMatrix(seg, t) = 0.4 * sin(2 * pi * waveSpeed * t - (seg-1) * phaseOffset);
    end
end

% ==== Setup Figure ====
figure('Color', 'w');
wormLine = plot(nan, nan, '-o', 'LineWidth', 3, 'Color', [0.2 0.6 1]);
axis equal;
xlim([-numSegments, numSegments]);
ylim([-numSegments, numSegments]);
xlabel('X'); ylabel('Y');
title('Worm Bending via Motor Neuron Activation');
grid on;

% ==== Animate Worm ====
for t = 1:T
    angles = activationMatrix(:, t);  % bend angles from motor output
    x = zeros(numSegments+1, 1);
    y = zeros(numSegments+1, 1);
    theta = 0;

    for i = 1:numSegments
        theta = theta + angles(i);  % cumulative angle (for curvature)
        x(i+1) = x(i) + segmentLength * cos(theta);
        y(i+1) = y(i) + segmentLength * sin(theta);
    end

    % Update plot
    set(wormLine, 'XData', x, 'YData', y);
    title(sprintf('Motor Neuron Driven Worm — Frame %d', t));
    pause(0.03);
    drawnow;
end

%%

T = 100;
numSegments = 5;
segmentLength = 1.0;
motor_activity = zeros(numSegments, T);

% Simulate noisy neural output for each segment
for i = 1:numSegments
    noise = 0.1 * randn(1, T);
    motor_activity(i, :) = 0.5 + 0.5 * tanh(cumsum(noise));  % smooth random activation
end

% Setup figure
figure('Color', 'w');
wormLine = plot(nan, nan, '-o', 'LineWidth', 3, 'Color', [0.9 0.4 0.2]);
axis equal;
xlim([-numSegments, numSegments]);
ylim([-numSegments, numSegments]);
xlabel('X'); ylabel('Y');
title('Worm Bending — Realistic Neural Input');
grid on;

% Animate
for t = 1:T
    angles = motor_activity(:, t);  % activation = bend angle
    x = zeros(numSegments+1, 1);
    y = zeros(numSegments+1, 1);
    theta = 0;

    for i = 1:numSegments
        theta = theta + (angles(i) - 0.5) * pi;  % map [0,1] to [-π/2, π/2]
        x(i+1) = x(i) + segmentLength * cos(theta);
        y(i+1) = y(i) + segmentLength * sin(theta);
    end

    set(wormLine, 'XData', x, 'YData', y);
    title(sprintf('Realistic Motor Neuron Bending — Frame %d', t));
    pause(0.03);
    drawnow;
end