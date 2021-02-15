function fexo=AddElemVar(fexo,varname,blkidx,data)
    %
    %  fexo=AddElemVar(fexo,varname,blkidx,data)
    %
    %  where
    %    fexo => the output of exo_get
    %    varname => name of variable described by the variable data (string)
    %    blkidx  => index of block that the variable data is affiliated with
    %    data => nelements x ntime matrix of element data
    %
    %  fexo.Time should be defined before calling this function
    %
    %
    %  This function will add one data to one element variable for one block.
    %    Call it repeatedly to add more data.  If the element variable has not
    %    been used before, this function will create it.  If it exists, this 
    %    function will modify ElemVars appropriately.
    %
    nblks=length(fexo.Blocks);
    if length(fexo.Time)~=size(data,2),
        error('Number of columns of data is not equal to the number of timesteps')
    end
    if size(fexo.Blocks(blkidx).Connectivity,1)~=size(data,1),
        error('Number of columns of data is not equal to the number of elements in block %d',fexo.Blocks(blkidx).ID);
    end
    
    if ~isfield(fexo,'ElemVars'),
        for i=1:nblks,
            fexo.ElemVars(i).BlockID=fexo.Blocks(i).ID;
            fexo.ElemVars(i).ElemVarData.Name=varname;
            if blkidx==i,
                fexo.ElemVars(i).ElemVarData.Data=data;
            else
                fexo.ElemVars(i).ElemVarData.Data=[];
            end
        end
    else
        names={fexo.ElemVars(1).ElemVarData.Name};
        nevar=length(names);
        idx=strmatch(varname,names,'exact');
        if ~isempty(idx),
            fexo.ElemVars(blkidx).ElemVarData(idx).Data=data;
        else
            for i=1:nblks,
                fexo.ElemVars(i).ElemVarData(nevar+1).Name=varname;
                if i==blkidx,
                    fexo.ElemVars(i).ElemVarData(nevar+1).Data=data;
                else
                    fexo.ElemVars(i).ElemVarData(nevar+1).Data=[];
                end
            end
        end
    end
end
