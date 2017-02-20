function geneidices= featurefind(Gi,featuretypes,genenames,timeout)
% created by Xiaoting XU, xiaotingxu@gmail.com
% modified 2014/06/30
geneidices = [];
% fprintf('\nreading sequence from genbank. Gi: %s\n',Gi);
% site='http://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&sendto=on&log$=seqview&db=nuccore&dopt=genbank&sort=ORGN&query_key=1&qty=910&filter=all'
retrieveurls = ['https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore'...
            '&id=' char(Gi) '&rettype=gb&retmode=text'];
internet_wait=0;
seq=[];
while internet_wait<=timeout&&isempty(seq)
    try 
        gbtext=urlread(retrieveurls);
        seq = genbankread(gbtext);
    catch exception
        if strcmp(exception.identifier,'MATLAB:urlread:ConnectionFailed')||...
                strcmp(exception.message,'FILE is not a valid GenBank file or url.')||...
                strcmp(exception.message,'The IP address of "eutils.ncbi.nlm.nih.gov" could not be determined.')
            display(['canot read url: '  retrieveurls]);
            display('try again');
%             pause(5);
            internet_wait=internet_wait+5;
        else
            display(Gi)
            rethrow(exception);
        end
    end
end
if isempty(seq)&&internet_wait==timeout
    error('Time out. The service is not available now or your internet connection is broken');
end
n=length(seq);
for t=1:n
tt = featuresparse(seq(t));
m=1;
while isempty(geneidices)&&m<=numel(featuretypes)
    if isfield(tt, featuretypes{m})
        tagfeature=tt.(featuretypes{m});
        xnames = fieldnames(tagfeature);
        genenames=cellfun(@(x) x(~isspace(x)),genenames, 'uni',0);
        for i = 1:numel(tagfeature)
            for j = 3:numel(xnames)
                if ~isempty(tagfeature(i).(xnames{j})) && ischar(tagfeature(i).(xnames{j}))
                    note = tagfeature(i).(xnames{j});
                    x = cellfun(@(x) strfind(note(~isspace(note)),x),genenames,'UniformOutput',false);
                    if any(cell2mat(x))
                        geneidices = tagfeature(i).Indices;
                        return
                    end
                end
            end
        end
    end
    m=m+1;
end
end
