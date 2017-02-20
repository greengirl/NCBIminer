function [Organism,Accs,seq,specimen_voucher,country,lat_lon,collection_date,collected_by,...
    identified_by]=getncbiseqs(giID,timeout,varargin)
% created by Xiaoting Xu.
% modified in 2014/06/30
%start position, stop position and strand (1 original sequence and 2 reverse completment sequence)
% start position, stop position and strand (1 original sequence and 2 reverse completment sequence)

if numel(varargin)&&numel(varargin{1})>=2
    partialseq=varargin{1};
end
db = 'nucleotide';

retrieveurls = ['https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=' char(db)...
            '&id=' char(giID) '&rettype=gb&retmode=text'];
if exist('partialseq','var')
    retrieveurls = [retrieveurls '&seq_start=',num2str(partialseq(1)),'&seq_stop=',num2str(partialseq(2))];
    if size(partialseq,2)==3
            retrieveurls=[retrieveurls '&strand=',num2str(partialseq(3))];
    end
end

% fprintf(strcat('\nread sequence from genbank:',char(giID)))
interwait=0;
data=[];
while interwait<=timeout&&isempty(data)
try
    gbtext = urlread(retrieveurls);
    gbtext=regexprep(gbtext,'\t','');
    data=genbankread(gbtext);
catch exception
    if strcmp(exception.identifier,'MATLAB:urlread:ConnectionFailed')||...
            strcmp(exception.message,'FILE is not a valid GenBank file or url.')||...
         	strcmp(exception.message,'The IP address of "eutils.ncbi.nlm.nih.gov" could not be determined.')
        display(['cannot read url: '  retrieveurls]);
      	display('try again');
        pause(5);
       	interwait=interwait+5;
    elseif interwait==timeout
        error('time out!');
    else
        display(['can not read GI: ' Gi]);
        rethrow(exception);
    end
end
end
n=length(data);
Organism=cell(n,1);
Accs=cell(n,1);
seq=cell(n,1);
for i=1:n
    Organism{i}=strtrim(data(i).SourceOrganism(1,:));
    Accs{i}=strtrim(data(i).Version);
%     if strcmp(data(i).Version,giID)
    seq{i}=data(i).Sequence;
    if nargout==9
        if i==1
            specimen_voucher=cell(n,1);
            country=cell(n,1);
            lat_lon=cell(n,1);
            collection_date=cell(n,1);
            collected_by=cell(n,1);
            identified_by=cell(n,1);
        end
        features=data(i).Features;
        okargs={'/specimen_voucher=','/country=','/lat_lon=','/collection_date=',...
            '/collected_by=','/identified_by='};
        ln=2;
        tln=size(features,1);
        while strncmpi(features(ln,:),repmat(' ',1,10),10)
            tmp=strtrim(features(ln,:));
            id = regexp(tmp,'=');
            if ~isempty(id)
                pname = tmp(1:id(1));
                pval = tmp(id(1)+1:length(tmp));
                k = find(strncmpi(pname,okargs,numel(pname)));
                if any(k)
                    switch k
                        case 1
                            specimen_voucher{i}=pval;
                        case 2
                            country{i}=pval;
                        case 3
                            lat_lon{i}=pval;
                        case 4
                            collection_date{i}=pval;
                        case 5
                            collected_by{i}=pval;
                        case 6
                            identified_by{i}=pval;
                    end
                end
            end
            if tln==ln
                break;
            end
            ln=ln+1;
        end
    end
end

        
