function NCBIminer
% main function to run NCBIminer
% Created by Xiaoting Xu 2014/06/29 (xiaotingxu@gmail.ocm)
% for linux 2014/07/14: file input through command line.
% for MAC 2014/07/14: file input through command line or by dobule click
% the application
%  Citation: Xu, X., Dimitrov, D., Rahbek, C., Wang, Z. 2014. NCBIminer: Sequences harvest from Genbank. 
%  Ecography (ver. 0). doi: 10.1111/ecog.01055
% Lasted updation: 2015/12/02

fprintf('\nProgramm start\n');
% mafftdir=dirmafft;
try
    if ismac||ispc
        fprintf('please select an input file\n');
        [FileName,PathName]=uigetfile('*.txt');% select an input file
    elseif isunix
        prompt = 'Where is the input file, e.g. ''/home/NCBIminer/'':';
        PathName = input(prompt);
        prompt = 'input file name, e.g. ''Demo1.txt'':';
        FileName = input(prompt);
        if ~strcmp(PathName(length(PathName)),'/')
            PathName=strcat(PathName,'/');
        end
    else
        Error('Unknow operating system'); % read jobs
    end
    cd(PathName);  % output folder
    [featuretype,featurename,refseq,ExpectValue,seqlen,taxlist_txt,timeout,Alignments,...
        Entrezs,fprefix,grpdist,location,extseq,database]=readjobs_old([PathName,FileName]); % read jobs
tic;
% if matlabpool('size')<=0
%     matlabpool open
% end

BlastProgram='blastn';
for fileid=1:numel(featuretype) % how many files
    filename=[fprefix{fileid} '_log.txt'];
%     filename1=[PathName fprefix{fileid} '.fas'];
%     if exist(filename,'file')||exist(filename1,'file')
%         fprintf(['\n"%s" or ''%s'' exsit, go to the next. If you want to renew it,\n ', ...
%             'please delete it first\n'],filename, filename1);
%         continue;
%     end
    fid=fopen(filename,'a+');
    refname=[PathName 'ref_' fprefix{fileid} '.fas'];
    tax_lev1i=taxlist_txt{fileid};
    for i = 1:size(tax_lev1i,1) % how many unique big groups (taxon before '\' in input files)
        tstart=tic;
        m=0;n=0;blastres2=[];
        taxname=tax_lev1i{i,1}(1:length(tax_lev1i{i,1})-6); 
        Entrez=['(' tax_lev1i{i,1} ') ' Entrezs{fileid}];
        % iterated BLAST for finding new queries
        [newrefseq, blastres, nref]=blastref(refseq{fileid},BlastProgram,Entrez,...
            Alignments(fileid),ExpectValue(fileid),timeout(fileid),featuretype{fileid},...
            featurename{fileid},grpdist(fileid),database(fileid));
       	% final blast
        if ~isempty(newrefseq)
            m=numel(newrefseq);
            fastawrite(refname,newrefseq);
           	Entrez=['(' tax_lev1i{i,2} ') '  Entrezs{fileid}];
            if strcmp(tax_lev1i{i,2},tax_lev1i{i,1})
                newrefid=setdiff(1:m,nref);
                if isempty(newrefid)
                    blastres2{1}=blastres;
                else
                    blastres2=blast_final(newrefseq(newrefid),BlastProgram,Entrez,Alignments(fileid),ExpectValue(fileid),timeout(fileid),database(fileid));
                    blastres2{length(blastres2)+1}=blastres;
                end
            else
                blastres2=blast_final(newrefseq,BlastProgram,Entrez,Alignments(fileid),ExpectValue(fileid),timeout(fileid),database(fileid));
            end
        end
        fprintf(fid,'%s: %d refseqs; ',taxname, m);
        fprintf('%s: %d refseqs; ',taxname, m);
        % combine blast results and read from genbank
        if ~isempty(blastres2)
          	[Header,Sequences,colltable]=readblastres(blastres2,seqlen(fileid,:),location(fileid),timeout(fileid),extseq(fileid,:));
           	n=numel(Header);
          	if n
                fastawrite([PathName fprefix{fileid} '.fas'], Header, Sequences);
                colltable(:,length(colltable(1,:))+1)=cellstr(repmat(taxname,n,1));
                colltable=colltable';
                if location(fileid)
                    if ~exist([PathName fprefix{fileid} '_table.txt'],'file')
                        exportid=fopen([PathName fprefix{fileid} '_table.txt'],'a+');
                        fprintf(exportid,['Accession_number\tOrganism\tspecimen_voucher\tcountry\tlat_lon\t'...
                            'collection_date\tcollected_by\tidentified_by\tstart\tstop\tstrand\ttax_lev1i\n']);
                        % updated on 2/12/2015 NCBIminer v1.02
                        fprintf(exportid,[repmat('%s\t',1,8),'%d\t%d\t%d\t','%s\n'],colltable{:,:});
                    else
                        exportid=fopen([PathName fprefix{fileid} '_table.txt'],'a+');
                        fprintf(exportid,[repmat('%s\t',1,8),'%d\t%d\t%d\t','%s\n'],colltable{:,:});
                    end
                else
                    if ~exist([PathName fprefix{fileid} '_table.txt'],'file')
                        exportid=fopen([PathName fprefix{fileid} '_table.txt'],'a+');
                        fprintf(exportid,'Accession_number\tOrganism\tstart\tstop\tstrand\ttax_lev1i\n');
                         % updated on 2/12/2015 NCBIminer v1.02
                        fprintf(exportid,[repmat('%s\t',1,2),'%d\t%d\t%d\t','%s\n'],colltable{:,:});
                    else
                        exportid=fopen([PathName fprefix{fileid} '_table.txt'],'a+');
                        fprintf(exportid,[repmat('%s\t',1,2),'%d\t%d\t%d\t','%s\n'],colltable{:,:});
                    end
                end
                fclose(exportid);
            end
        end
       	a = toc(tstart); % total time
     	fprintf('%d sequences were found for %s\n',n,taxname);
     	fprintf(fid,'\n%d sequences were found for %s. Calculation time: %10.5f s\n',n,featurename{fileid}{1},a);
        clear refseqs taxlisti m n Entrez taxoni taxlistij blastres2 blastres newrefseq
    end
    fclose(fid);
end
% fprintf(fid,'Calculation time: %10.5f s\nPress any key to exit the program\n',toc);
fprintf('Calculation time: %10.5f s\nPress any key to exit the program\n',toc);
catch err
    n=length(err.stack);
    display(err.message)
    if ~exist('fid','var')
        filename='NCBIminer_log.txt';
        fid=fopen(filename,'a+');
    end
    if exist('fileid','var')
        fprintf('Error in %s (line %d)\n',FileName,fileid);
        fprintf(fid,'Error in %s (line %d)\n',FileName,fileid);
    end
    fprintf(fid,'\n%s\n',err.message); %2015-03-06
    for i=1:n
        fprintf('Error in %s (line %d)\n',err.stack(i).file,err.stack(i).line);
        fprintf(fid,'Error in %s (line %d)\n',err.stack(i).file,err.stack(i).line);
    end
    fclose(fid);
end
if exist('err','var')
    fprintf('Job finish but there is an error, please check the log fiel: %s\n', filename);
    display('Press any key to exit the program');
end
% if matlabpool('size')>0
%     matlabpool close
% end
pause

