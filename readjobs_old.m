function [featuretype,featurename,Refseq,ExpectValue,seqlen,Tax,timeout,Alignments,...
    Entrezs,fprefix,grpdist,location,extseq,database]=readjobs_old(varargin)
% varargin: input file name with full path
if numel(varargin)
%     fid=fopen(varargin{1});
    FileName=varargin{1};
else
    [FileName,PathName]=uigetfile('*.txt'); % open input file
     FileName=[PathName FileName];
end

fid=fopen(FileName);
n=0;
tline1=fgetl(fid);
if exist(tline1,'file')
    filename=cell(1,1);
    while ischar(tline1)
        n=n+1;
        filename{n}=tline1; 
        tline1=fgetl(fid);
    end
else
    filename{1}=FileName;
    n=1;
end
fclose(fid);
featurename=cell(n,1);
featuretype=cell(n,1);
Entrezs=cell(n,1);
seqlen=zeros(n,2);
timeout=zeros(n,1);
Alignments=zeros(n,1);
ExpectValue=zeros(n,1);
grpdist=zeros(n,1);
location=zeros(n,1);
extseq=zeros(n,2);
tax=cell(n,1);
Refseq=cell(n,1);
fprefix=cell(n,1);
database=cell(n,1);
keywords={'feature type:' 'feature name:' 'initial query:','ExpectValue:', 'sequence length:' 'taxon list:' 'database:',...
    'time out:','Alignments:','Entrezs:','output prefix:','maximum query sequences difference:','location:','extended length:'};
for i=1:n
    fid=fopen(filename{i});
    numref = 0;
    t = 1;
    featurename{i} = [];
    featuretype{i} = [];
    Entrezs{i}=[];
    seqlen(i,:)=[50,inf];
    timeout(i)=1800;
    Alignments(i)=5000;
    ExpectValue(i)=0.0000000001;
    grpdist(i)=0.3;
    refseq(1).Header='';
    refseq(1).Sequence='';
    database{i}={'nr'};
    [tline1, x] = readline(fid,keywords);
    tnum=1;
    while ischar(tline1)
        while sum(x)
        [tline2, y] = readline(fid,keywords); % read te
            while sum(y)~=1
                switch tline1
                    case 'feature type:'
                        featuretype{i} = [featuretype{i}; cellstr(tline2)];
                    case 'feature name:'
                        featurename{i} = [featurename{i}; cellstr(tline2)];
                    case 'initial query:'
                        if matchstart(tline2,'>')
                            numref=numref+1;
                            refseq(numref).Header=tline2;
                        elseif exist(tline2,'file')
                            refseq=fastaread(tline2);
                        elseif numref
                            refseq(numref).Sequence = strcat(refseq(numref).Sequence,tline2);         
                        end
                    case 'sequence length:'
                        tline2=regexp(tline2,',','split');
                        seqlen(i,1)=str2double(strtrim(tline2{1}));
                        seqlen(i,2)=str2double(strtrim(tline2{2}));
                    case 'taxon list:'
                        taxlist=regexp(tline2,'\','split');
                        if t==1
                            tax={[taxlist{1} '[ORGN]'] ''};
                        elseif ~strcmp(tax(tnum,1),[taxlist{1} '[ORGN]'])
                            tnum=tnum+1;
                            tax{tnum,1} = [taxlist{1} '[ORGN]'];
                            tax{tnum,2} = '';
                        end
                        tax{tnum,2} = [tax{tnum,2} taxlist{2} '[ORGN] OR '];
                        t=t+1;
                    case 'time out:'
                        timeout(i)=str2double(tline2);
                    case 'Alignments:'
                        Alignments(i)=str2double(tline2);
                    case 'ExpectValue:'
                        ExpectValue(i)=str2double(tline2);
                    case 'Entrezs:'
                        Entrezs{i}=tline2;
                    case 'output prefix:'
                        fprefix{i}=tline2;
                    case 'maximum query sequences difference:'
                        grpdist(i)=str2double(tline2);
                    case 'location:'
                        location(i)=str2double(tline2);
                    case 'extended length:'
                        tline2=regexp(tline2,',','split');
                        extseq(i,1)=str2double(strtrim(tline2{1}));
                        extseq(i,2)=str2double(strtrim(tline2{2}));
                     case 'database:'
                         database{i}=regexp(tline2,',','split');   
                end
                [tline2, y] = readline(fid,keywords); % read te
                if ~ischar(tline2)
                    break;
                end
            end
        x=y;
        tline1=tline2;
        end
    [tline1, x] = readline(fid,keywords);
    end
    fclose(fid);
    tax(1:tnum,2) = cellfun(@(x) cellstr(x(1:length(x)-3)),tax(1:tnum,2));
    Refseq{i}=refseq;
    Tax{i}=tax;
end
if any(cellfun(@isempty,featuretype))
    error('myApp:argChk','feature type is empty! Please check your input file');
end
if any(cellfun(@isempty,featurename))
    error('myApp:argChk','feature name is empty! Please check your input file');
end
if any(cellfun(@isempty,Refseq))
    error('myApp:argChk','initial query is empty! Please check your input file');
end
if any(cellfun(@isempty,Tax))
    error('myApp:argChk','taxon list is empty! Please check your input file');
end
if any(cellfun(@isempty,database))
    error('myApp:argChk','database is empty! Please check your input file');
end
fprintf('\nread jobs succefully\n');
end

function [tline, keyword] = readline(fid,keywords)
% igore blank line and comments
keyword=[];
tline = '';
while isempty(tline)
    tline = fgetl(fid);
    if ischar(tline)
        tline = strtrim(tline);
        [x, y] = regexp(tline,'#','split');
        if ~isempty(y) 
            tline = strtrim(x{1});
        end
        if ~isempty(tline)
            if strcmpi(tline(length(tline)),':')
                keyword = strcmpi(tline,keywords);
                if sum(keyword)==0
                    errmessage=[tline ' ' 'is not a valid parameter!'];
                    error('myApp:argChk',errmessage);
                end
            end
        end
    end
end
end