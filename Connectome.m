%% Load data
%load("Connectome\GHermChem.mat")
%load("Connectome\GHermElec_sym.mat")
%%
t = tiledlayout(1,2);
fig1 = nexttile(t, 1);
fig1 = HermChemNodeConns(["dBWML10", "dBWMR10", "vBWML10", "vBWMR10"],1);

%%
fig2 = nexttile(t, 2);
fig2 = HermElecNodeConns(["dBWML10", "dBWMR10", "vBWML10", "vBWMR10"],1);

%% Maps all predecessors of Target

%StartNode = GHermChem.Nodes.Name(GHermChem.findnode("dBWML1"));
function  nodefig = HermElecNodeConns(StartNode, depth)
    load("Connectome\GHermElec_sym.mat")
    
    %StartNode = BodyWallMuscles(6:10);
    PredNode = StartNode;
    NewPredNode = [];
    FiltNode = StartNode;
    FiltNeuron = StartNode;
    %depth = 1;
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
    GHermLastFirst = subgraph(GHermElec_Sym,unique([StartNode, NewPredNode.']));
    GHermOnlyNeuron = subgraph(GHermMusc, FiltNeuron);
    
    numInterNeurons = sum(matches(FiltNode, InterNeurons));
    numMotorNeurons = sum(matches(FiltNode, MotorNeurons));
    numSensoryNeurons = sum(matches(FiltNode, SensoryNeurons));
    numNodes = size(FiltNode,2);
    numNeurons = numInterNeurons + numMotorNeurons + numSensoryNeurons;
    % find predecessor nodes (PredNode)
    % save pred_nodes to filter nodes (FiltNode)
    % for each pred_node, find pred node (NewPredNode)
    % repeat until we get no new nodes
    %
    
    GHermPlot = GHermMusc; % reassign so we only need to change digraph once
    
    % green for bodywall muscles, y for interneurons, b for motor, r for
    % sensory
    NodeColor = [];
    for i = 1:size(GHermPlot.Nodes,1)
        if(matches(GHermPlot.Nodes.Group(i), 'BodyWallMuscles'))
            NodeColor(i,:) = 	[0 1 0];
        elseif(matches(GHermPlot.Nodes.Group(i), 'InterNeurons'))
            NodeColor(i,:) = [1 1 0];
        elseif(matches(GHermPlot.Nodes.Group(i), 'MotorNeurons'))
            NodeColor(i,:) = [0 0 1];
        elseif(matches(GHermPlot.Nodes.Group(i), 'SensoryNeurons'))
            NodeColor(i,:) = [1 0 0];
        end
    end
    LWidths = 5*GHermPlot.Edges.Weight/max(GHermElec_Sym.Edges.Weight);
    
    %nodefig = figure;
    nodefig = plot(GHermPlot, 'NodeColor', NodeColor, 'LineWidth',LWidths);
end
%% Maps all predecessors of Target

%StartNode = GHermChem.Nodes.Name(GHermChem.findnode("dBWML1"));
function  nodefig = HermChemNodeConns(StartNode, depth)
%    load("Connectome\GHermChem.mat")
    
    %StartNode = BodyWallMuscles(6:10);
    PredNode = StartNode;
    NewPredNode = [];
    FiltNode = StartNode;
    FiltNeuron = StartNode;
    %depth = 1;
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
    GHermLastFirst = subgraph(GHermChem,[StartNode, NewPredNode.']);
    GHermOnlyNeuron = subgraph(GHermMuscChem, FiltNeuron);
    
    numInterNeurons = sum(matches(FiltNode, InterNeurons));
    numMotorNeurons = sum(matches(FiltNode, MotorNeurons));
    numSensoryNeurons = sum(matches(FiltNode, SensoryNeurons));
    numNodes = size(FiltNode,2);
    numNeurons = numInterNeurons + numMotorNeurons + numSensoryNeurons;
    % find predecessor nodes (PredNode)
    % save pred_nodes to filter nodes (FiltNode)
    % for each pred_node, find pred node (NewPredNode)
    % repeat until we get no new nodes
    %
    
    GHermPlot = GHermOnlyNeuron; % reassign so we only need to change digraph once
    
    % green for bodywall muscles, y for interneurons, b for motor, r for
    % sensory
    NodeColor = [];
    for i = 1:size(GHermPlot.Nodes,1)
        if(matches(GHermPlot.Nodes.Group(i), 'BodyWallMuscles'))
            NodeColor(i,:) = 	[0 1 0];
        elseif(matches(GHermPlot.Nodes.Group(i), 'InterNeurons'))
            NodeColor(i,:) = [1 1 0];
        elseif(matches(GHermPlot.Nodes.Group(i), 'MotorNeurons'))
            NodeColor(i,:) = [0 0 1];
        elseif(matches(GHermPlot.Nodes.Group(i), 'SensoryNeurons'))
            NodeColor(i,:) = [1 0 0];
        end
    end
    LWidths = 5*GHermPlot.Edges.Weight/max(GHermChem.Edges.Weight);
    
    %nodefig = figure;
    nodefig = plot(GHermPlot, 'NodeColor', NodeColor, 'LineWidth',LWidths);
end

%% test for successors of BWMs
% test = 0;
% for i = 1:95
%     test = [test; successors(GHermChem, BodyWallMuscles(i))];
% end
