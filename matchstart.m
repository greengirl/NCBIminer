function b=matchstart(s,exps)
s=strtrim(s);
	[~,a]=regexp(s,strcat('^',exps),'once');
    if isempty(a)
        b=0;
    else
        b=1;
    end
end