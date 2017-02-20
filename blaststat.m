function data = blaststat(rid,desc,timeout)
% modified by Xiaoting Xu from the MATLAB function "getblast". blaststat
% only read tabulated BLAST results without alignments

% Child function blast_sta parse the tabulated results into structure with
% filds 'Gi', 'alignlen'(alignment length), 'mismatches',
% 'q_start'(q. start), 'q_end'(q. end), 's_start'(s. start), 's_end'(s. end),
% 'Expect'(expection value), 'Score', 'Strand'(1 means positive oritation while 
% %  2 means reversed). 
% 

%   parameters:
%   desc, specified number of descriptions in the report. Acceptable values are 1-100000.

%   Modified from fucntion "getblast" in MATLAB bioinformatic toolbox.


% verify java is available
if (~usejava('jvm'))
    error('No Java. Please install it')
end

% Set default parameters to retrieve the BLAST report
if desc<1&&desc>100000
    error('desc is a number between 1-100000');
end
site = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=';
paras=['&ALIGNMENT_VIEW=Tabular&FORMAT_TYPE=Text&DESCRIPTIONS=' num2str(desc)];

% Check for a valid RID type and format
if iscellstr(rid)
    rid=char(rid);
end
%if regexpi(rid,'\d+-\d+-\d+\.BLASTQ\d')
%if regexpi(rid, '\<(\w{11})\>') % change in RID format April 16 2007
if regexpi(rid, '\<(\w+)\>') % change in RID format April 16 2008
    texturl = [site rid paras];
else
    error('invalid rid');
end

% polls the NCBI BLAST SERVER for results
searchurl = [site rid '&FORMAT_OBJECT=SearchInfo'];
wait_timer = 0;
blaststa=[];
data=[];
while wait_timer <= timeout||~isempty(blaststa)
	searchInfo = urlread(searchurl);
	start = regexp(searchInfo, 'QBlastInfoBegin');
	stop = regexp(searchInfo, 'QBlastInfoEnd');
    if isempty(start) || isempty(stop)
        error('');
    end
	status = regexpi(searchInfo(start(1):stop(1)), '(?<=Status=)\w+', 'match', 'once');
    if strcmp(status, 'READY')
        disp 'Reading Blast results. Please wait ...';
        blaststa = urlread(texturl);
        blaststa = strrep(blaststa,char(0),'a'); % Correct invalid character in returned text
        if isempty(regexpi(searchInfo(start(2):stop(2)), 'ThereAreHits=yes'))
            disp (['No singificant hits for' rid]);
        else
            data=blast_sta(blaststa);%read tabulate BLAST results. new function written by Xiaoting Xu
        end
        return;
    elseif strcmp(status, 'WAITING')
        disp 'Blast results are not available yet. Please wait ...';
    elseif strcmp(status, 'FAILED')
        disp ['Search failed, try again in 5 seconds'];
    elseif strcmp(status, 'UNKNOWN')
        disp ['Search may be expired, try again in 5 seconds'];
    else
        disp(['Unknow status for search ' rid]);
    end
    pause(8);
    wait_timer = wait_timer + 8;
end

if wait_timer >= timeout
    error('exec-timeout!')
end 
end

function data=blast_sta(blaststa1)
    data=[];
    id1 = regexp(blaststa1,' hits found');
    id2 = regexp(blaststa1,'#');
    numhits = str2double(blaststa1(id2(numel(id2))+2:id1-1)); 
    id2=regexp(blaststa1,'</PRE>');
    blaststa=blaststa1(id1+length(' hits found'):id2-1);%Delete header £¨7lines£© and blank lines
    hitbegin = regexp(blaststa,'\d\n');
    if numhits~=numel(hitbegin)
        error('blast report is not a tublar');
    end
%     Fields={'Name', 'identity', 'alignment length', 'mismatches',...
%         'gap opens', 'q. start', 'q. end', 's. start', 's. end', 'Expect', 'Score'};
%     m=numel(Fields);
    blaststa=regexprep(blaststa,'\n','	');% replacy newline with table seprator
    blaststa=regexprep(blaststa,'\n','');%delete extra newline
    tags1 = find(blaststa == '	')+1;
    tags2 = [tags1(2:numel(tags1))-1 numel(blaststa)];
    j=1;
    for i = 1:numhits
        line=(i-1)*14+1;
        Gi=blaststa(tags1(line+1):tags2(line+1));
        x = regexp(Gi,';','split');
        Gi=regexp(x{1},'([\d]*)','tokens','once');
        data.stat(j).Gi=Gi{1};
        % modified on 17/11/2015
%         Accession=regexp(x{1},'\|(\w*\d*\.\d)\|','tokens','once');
        Accession=regexp(x{1},'\|(\w*\d*\.\d*)\|','tokens','once');
        data.stat(j).Accession=Accession{1};
        data.stat(j).identity=str2double(blaststa(tags1(line+4):tags2(line+4)));
        data.stat(j).alignlen=str2double(blaststa(tags1(line+5):tags2(line+5)));
        data.stat(j).mismatches=str2double(blaststa(tags1(line+6):tags2(line+6)));
        data.stat(j).q_start=str2double(blaststa(tags1(line+7):tags2(line+7)));
        data.stat(j).q_end=str2double(blaststa(tags1(line+9):tags2(line+9)));
        data.stat(j).s_start=str2double(blaststa(tags1(line+10):tags2(line+10)));
        data.stat(j).s_end=str2double(blaststa(tags1(line+11):tags2(line+11)));
        data.stat(j).Expect=str2double(blaststa(tags1(line+12):tags2(line+12)));
        data.stat(j).Score=str2double(blaststa(tags1(line+13):tags2(line+13)));
        if data.stat(j).s_start>data.stat(j).s_end
            data.stat(j).Strand=2;
        else
            data.stat(j).Strand=1;
        end
        a=1;
        while numel(x)>a
            data.stat(j+a)=data.stat(j);
            Gi=regexp(x{a+1},'([\d]*)','tokens','once');
            data.stat(j+a).Gi=Gi{1};
            Accession=regexp(x{a+1},'\|(\w*\d*\.\d)\|','tokens','once');
            data.stat(j+a).Accession=Accession{1};
            a=a+1;
        end    
        j=j+a;
    end
    data=data.stat;
end
