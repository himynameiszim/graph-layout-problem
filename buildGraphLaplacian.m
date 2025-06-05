% buildGraphLaplacian.m
% build the graph Laplacian submatrix for unknown vertices,
% along with full graph degree and adjacency information.

function [L_mat, degrees_full, adjacency_full] = buildGraphLaplacian(total_num_vertices, edges_list, unknown_idx_list)

    degrees_full = zeros(total_num_vertices, 1); % degree matrix
    adjacency_full = cell(total_num_vertices, 1); % adjacency matrix
    m_unknown = length(unknown_idx_list); % number of unknown vertex
    L_mat = zeros(m_unknown, m_unknown); % graph Laplacian matrix

    % init adjacency matrix
    for k = 1:total_num_vertices
        adjacency_full{k} = [];
    end

    % construct adjacency matrix and degree matrix
    for k = 1:size(edges_list, 1)

        degrees_full(edges_list(k,1)) = degrees_full(edges_list(k,1)) + 1;
        degrees_full(edges_list(k,2)) = degrees_full(edges_list(k,2)) + 1;

        adjacency_full{edges_list(k,1)}(end+1) = edges_list(k,2);
        adjacency_full{edges_list(k,2)}(end+1) = edges_list(k,1);
    end

    % map unknown vertices to the graph Laplacian matrix L
    unknown_map = containers.Map('KeyType', 'double', 'ValueType', 'double');
    for k=1:m_unknown
       unknown_map(unknown_idx_list(k)) = k; 
    end

    % construct graph Laplacian matrix L
    for k=1:m_unknown
        L_mat(k, k) = degrees_full(unknown_idx_list(k));

        % off-diagonal elements
        % if an unknown vertex i is adjacent to another unknown vertex j 
        % -> L(i,j) = -1, otherwise it is already set to 0
        k_neighbours = adjacency_full{unknown_idx_list(k)};
        for i=1:length(k_neighbours)
            if isKey(unknown_map, k_neighbours(i))
                j = unknown_map(k_neighbours(i));
                L_mat(k, j) = -1;
            end
        end
    end
end