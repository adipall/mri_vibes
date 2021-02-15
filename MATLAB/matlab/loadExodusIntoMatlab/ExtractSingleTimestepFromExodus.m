function ExtractSingleTimestepFromExodus(filenamein, filenameout, time, timestep,varNameMap)

if exist(filenameout,'file')
    disp(sprintf('Output file: %s already exists. Overwrite?',filenameout));
    reply = input('(y/n) ','s');
    
    if strcmp(lower(reply),'y') || strcmp(lower(reply),'yes')
        delete(filenameout);
    end
end

exofile = exo_rd(filenamein);



if nargin < 4 || timestep < 1
    timeVals = exofile.Time;
    diffs = timeVals - time;
    [~,timestep] = min(abs(diffs));
end

% operate on times
exofile.Time = exofile.Time(timestep);

% operate on nodal variables
for rowInd = 1:size(exofile.NodalVars,1)
    for colInd = 1:size(exofile.NodalVars,2)
        if ~isempty(exofile.NodalVars(rowInd,colInd).Data)
            exofile.NodalVars(rowInd,colInd).Data = exofile.NodalVars(rowInd,colInd).Data(:,timestep);
            
            if nargin > 4 
                for k = 1:length(varNameMap(1,:))
                    if strcmp(exofile.NodalVars(rowInd,colInd).Name,varNameMap{1,k})
                        exofile.NodalVars(rowInd,colInd).Name = varNameMap{2,k};
                        disp(sprintf('Replaced variable name %s with %s',varNameMap{1,k},varNameMap{2,k}));
                    end
                end
            end
        end
    end
end

% operate on global variables
for rowInd = 1:size(exofile.GlobalVars,1)
    for colInd = 1:size(exofile.GlobalVars,2)
        if ~isempty(exofile.GlobalVars(rowInd,colInd).Data)
            exofile.GlobalVars(rowInd,colInd).Data = exofile.GlobalVars(rowInd,colInd).Data(:,timestep);
            if nargin > 4 
                for k = 1:length(varNameMap(1,:))
                    if strcmp(exofile.GlobalVars(rowInd,colInd).Name,varNameMap{1,k})
                        exofile.GlobalVars(rowInd,colInd).Name = varNameMap{2,k};
                        disp(sprintf('Replaced variable name %s with %s',varNameMap{1,k},varNameMap{2,k}));
                    end
                end
            end
        end
    end
end

% operate on elem variables
for rowInd = 1:size(exofile.ElemVars,1)
    for colInd = 1:size(exofile.ElemVars,2)
        if ~isempty(exofile.ElemVars(rowInd,colInd).Data)
            exofile.ElemVars(rowInd,colInd).Data = exofile.ElemVars(rowInd,colInd).Data(:,timestep);
            if nargin > 4 
                for k = 1:length(varNameMap(1,:))
                    if strcmp(exofile.ElemVars(rowInd,colInd).Name,varNameMap{1,k})
                        exofile.ElemVars(rowInd,colInd).Name = varNameMap{2,k};
                        disp(sprintf('Replaced variable name %s with %s',varNameMap{1,k},varNameMap{2,k}));
                    end
                end
            end
        end
    end
end

% operate on nodeset variables
for rowInd = 1:size(exofile.NodesetVars,1)
    for colInd = 1:size(exofile.NodesetVars,2)
        if ~isempty(exofile.NodesetVars(rowInd,colInd).Data)
            exofile.NodesetVars(rowInd,colInd).Data = exofile.NodesetVars(rowInd,colInd).Data(:,timestep);
            if nargin > 4 
                for k = 1:length(varNameMap(1,:))
                    if strcmp(exofile.NodesetVars(rowInd,colInd).Name,varNameMap{1,k})
                        exofile.NodesetVars(rowInd,colInd).Name = varNameMap{2,k};
                        disp(sprintf('Replaced variable name %s with %s',varNameMap{1,k},varNameMap{2,k}));
                    end
                end
            end
        end
    end
end


% operate on sideset variables
for rowInd = 1:size(exofile.SidesetVars,1)
    for colInd = 1:size(exofile.SidesetVars,2)
        if ~isempty(exofile.SidesetVars(rowInd,colInd).Data)
            exofile.SidesetVars(rowInd,colInd).Data = exofile.SidesetVars(rowInd,colInd).Data(:,timestep);
            if nargin > 4 
                for k = 1:length(varNameMap(1,:))
                    if strcmp(exofile.SidesetVars(rowInd,colInd).Name,varNameMap{1,k})
                        exofile.SidesetVars(rowInd,colInd).Name = varNameMap{2,k};
                        disp(sprintf('Replaced variable name %s with %s',varNameMap{1,k},varNameMap{2,k}));
                    end
                end
            end
        end
    end
end

exo_wr(filenameout,exofile);