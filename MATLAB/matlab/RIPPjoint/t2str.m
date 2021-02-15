function str = t2str(time)
% Converts a time in seconds to a string

days = floor(time/60/60/24) ;
hours = floor(time/60/60-days*24) ;
minutes = floor(time/60 - hours*60 - days*24*60) ;
seconds = floor(time - minutes*60 -hours*60*60 - days*24*60*60) ;

if days > 0
    str = [num2str(days) ':' num2str(hours,'%02d') ':' num2str(minutes,'%02d') ':' num2str(seconds,'%02d')] ;
elseif hours > 0
    str = [num2str(hours) ':' num2str(minutes,'%02d') ':' num2str(seconds,'%02d')] ;
elseif minutes > 0
    str = [num2str(minutes) ':' num2str(seconds,'%02d')] ;
else
    str = [num2str(seconds) ' seconds'] ;
end
