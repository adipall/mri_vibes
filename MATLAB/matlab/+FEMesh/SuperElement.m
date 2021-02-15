classdef SuperElement
    properties
        Kr                  = []
        Cr                  = []
        Mr                  = []
        CBMap               = []
        OTM                 = []
        OutMap              = []
        OTME                = []
        OutElemMap          = []
        NumInterfaceNodes   = []
        NumEig              = []
        NumConstraints      = []
    end
    methods
        function obj=SuperElement(Kr,Cr,Mr,CBMap,OTM,OutMap,OTME,OutElemMap,NumInterfaceNodes,NumEig,NumConstraints)
            if nargin>0,
                obj.Kr=Kr;
                obj.Cr=Cr;
                obj.Mr=Mr;
                obj.CBMap=CBMap;
                obj.OTM=OTM;
                obj.OutMap=OutMap;
                obj.OTME=OTME;
                obj.OutElemMap=OutElemMap;
                obj.NumInterfaceNodes=NumInterfaceNodes;
                obj.NumEig=NumEig;
                obj.NumConstraints=NumConstraints;
            end
        end
        
    end
end