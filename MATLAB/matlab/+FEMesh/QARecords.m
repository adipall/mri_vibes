classdef QARecords < handle
    properties
        CodeName
        CodeVersion
        Date
        Time
    end
    methods
        function obj=QARecords(code,ver)
            if nargin>0,
                if iscell(code), %% this is from a file read
                    for i=1:length(code),
                        obj(i).CodeName=code{i}{1};
                        obj(i).CodeVersion=code{i}{2};
                        obj(i).Date=code{i}{3};
                        obj(i).Time=code{i}{4};
                    end
                else
                    obj.CodeName=code;
                    obj.CodeVersion=ver;
                    obj.Date=datestr(now,23);
                    obj.Time=datestr(now,13);
                end
            end
        end
        function display(obj)
            disp(sprintf('CodeName     CodeVersion       Date       Time'))
            for i=1:length(obj),
                disp(sprintf('%s         %s       %s            %s',obj(i).CodeName,obj(i).CodeVersion, ...
                obj(i).Date,obj(i).Time))
            end 
        end
        function qacell=outputQA(obj)
            qacell=cell(1,length(obj));
            for i=1:length(obj),
                qacell{i}{1}=obj(i).CodeName;
                qacell{i}{2}=obj(i).CodeVersion;
                qacell{i}{3}=obj(i).Date;
                qacell{i}{4}=obj(i).Time;
            end
        end
    end
end
