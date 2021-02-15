clear variables; close all; clc
addpath ../../

%% Read in the entire Exodus file
flag = 0;
% flag = 0 => everything gets read (default)
% flag = 1 => elemvars don't get read
% flag = 2 => elemvars and nodalvars don't get read
% flag = 3 => elemvars, nodalvars, and global vars don't get read
% flag = 4 => elemvars, nodalvars, global vars, and sidesets don't get read
% flag = 5 => elemvars, nodalvars, global vars, sidesets, and nodesets don't get read
exo_obj = exo_rd('example.exo',flag);
time = exo_obj.Time.'; % this might be frequency depending on the output file

%% Plot the AccY data for the all of the nodesets
nodalvar_2_plot = 'AccX';
nodalvar_index = exo_obj.NodalVars.isNodalVar(nodalvar_2_plot);

assert(nodalvar_index ~= 0); % if its 0 then nodalvar_2_plot isnt a NodalVar

num_nodesets = size(exo_obj.Nodesets, 2);

for ns_idx = 1:num_nodesets
    nodeset_id = exo_obj.Nodesets(ns_idx).ID;
    
    % get the row numbers of the data in the NodalVars array for each node
    % in the current nodeset
    row_indices_in_NodalVars = exo_obj.Nodesets.id2idx(nodeset_id);
    
    % Now we plot the data for each node in the current nodeset
    figure(ns_idx)
    for row_num = row_indices_in_NodalVars
        nodal_data = exo_obj.NodalVars(nodalvar_index).Data(row_num,:);
        plot(time, nodal_data, 'linewidth', 2)
        hold on
    end
    grid on
    title(sprintf('%s Data for Nodeset %d',nodalvar_2_plot, nodeset_id))
    xlabel('Time (s)')
    ylabel(nodalvar_2_plot)
    hold off
end

%% Get the element data VStressX for Block 4
block_id = 4;
elem_var_name = 'VStressX';
var_index = exo_obj.ElemVars.isElemVar(elem_var_name);

assert(var_index ~= 0); % otherwise elem_var_name is not an ElemVar
block_index = exo_obj.ElemVars.id2idx(block_id);

fprintf('Getting %s Data for Block Id %d\n',exo_obj.ElemVars(block_index, var_index).Name,...
    exo_obj.ElemVars(block_index, var_index).BlockID);

% This data is in the following format: row = element_index, column = time step index
elem_data_for_block_id = exo_obj.ElemVars(block_index, var_index).Data;

% If there is no element data for this block then we should give an error
%   Test this by making block_id = 2 or 3
if isempty(elem_data_for_block_id)
    error('There is no %s data for Block %d of type %s.',elem_var_name,...
        block_id, exo_obj.Blocks(block_index).ElementType)
end

% For fun, let's look at the maximum VStressX over time in this block
figure
plot(time, max(elem_data_for_block_id, [], 1), 'linewidth', 2)
grid on
xlabel('Time (s)')
ylabel(elem_var_name)
title(sprintf('Maximum %s in Block %d over time',elem_var_name, block_id))

%% Add a new NodalVar to the Exodus database and write out new Exodus file

% Lets scale AccX to be in g's rather than inches/second^2 and make it a new
% NodalVar
new_NodalVar_data = exo_obj.NodalVars(nodalvar_index).Data ./ 386.22;
new_NodalVar_name = 'AccX_in_Gs';

% This call assumes each column of data corresponds to the timestep in the
% Exodus file. Alternatively, you could write the data for a specific time
% by passing the time the column of data corresponds to into the 3rd arg
exo_obj.AddNodalVar(new_NodalVar_name, new_NodalVar_data);

new_exodus_file_name = 'example_new.exo';
% Write out the new Exodus file
exo_wr(new_exodus_file_name, exo_obj);