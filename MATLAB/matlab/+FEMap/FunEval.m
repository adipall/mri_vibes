classdef FunEval < FEMap.Map
    properties
        FunctionHandle
    end
    methods
        function obj=FunEval(fun,fexo)
            if nargin>0,
                obj.FunctionHandle=fun;
                obj.ModelObjTo=fexo;
            end
        end
        function obj=EvalNodal(obj,nname,time)
            
            nvar=obj.FunctionHandle(obj.ModelObjTo.Nodes.Coordinates, time);
            obj=obj.WriteNodal(nvar,nname,time);
        end
        function obj=EvalElem(obj,blkidx,gauss,ename,time)
            [enodexyz,nelem,gp_id,ngauss]=FEMap.Map.CalcElemCoords(obj.ModelObjTo,blkidx,gauss);
            evar=obj.FunctionHandle(enodexyz,time);
            obj=obj.WriteElem2(evar,ename,blkidx,nelem,ngauss,gp_id,time);
        end
    end
end