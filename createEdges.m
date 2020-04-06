function edges = createEdges(YieldNodes,nodes,nelem_x)
edges = []; 
 for i = 1 : size(YieldNodes, 1) - 1
    for j = i + 1 : size(YieldNodes, 1) 
        if sqrt((YieldNodes(i, 2) - YieldNodes(j, 2))^2 + (YieldNodes(i, 3) - YieldNodes(j, 3))^2)...
                <= 1.1*sqrt(nodes(2,2)^2 + nodes(nelem_x + 3,2)^2)
           edges = [edges; YieldNodes(i) YieldNodes(j) ...
               sqrt((YieldNodes(i, 2) - YieldNodes(j, 2))^2 + (YieldNodes(i, 3) - YieldNodes(j, 3))^2)];
        end    
    end
 end
end