function blastres = blast_new(refseqs, BlastProgram,Entrez,Alignments,ExpectValue,waiting,database)
m=numel(database);
blastnum=0;
dbnum=1;
tinterval=50;
while dbnum<=m
    internet_wait=0;
    rid=[];
    while internet_wait<=waiting&&isempty(rid)
        try
            fprintf('\nsubmiting blast, please wait.....\n');
            if tinterval<30
                pt=30-tinterval;
                pause(pt);
            	fprintf('two blasts interval less than 30s, please waiting for %ds\n',pt);
            end
            t1=tic;
            rid = blastncbi(refseqs, BlastProgram,'Database',database{1}{dbnum},'Entrez',Entrez,...
                'Description',Alignments,'Expect', ExpectValue);break;
        catch exception
            if strcmp(exception.identifier,'MATLAB:urlread:ConnectionFailed')||...
                strcmp(exception.message,'The IP address of "blast.ncbi.nlm.nih.gov" could not be determined.')
                disp('Wainting for internet connection, try in 10 seconds');
                pause(10);internet_wait=internet_wait+10;
            else 
                rethrow(exception);
            end
        end
    end
    if internet_wait>waiting&&isempty(rid)
        error('myApp:argChk','time out!');
    end
    internet_wait=0;
    while internet_wait<=waiting
        try
            blastrestmp = blaststat(rid,Alignments,waiting);break;
        catch theException
            if strcmp(theException.identifier,'MATLAB:urlread:ConnectionFailed')
%                 disp(theException.message);
                disp('Internet connection failed. Wainting for internet connection, try in 5 seconds');
                pause(5);internet_wait=internet_wait+5;
            % modified on 17/11/2015
            elseif strcmp(theException.message,'Search may be expired, try again in 5 seconds')
               	disp ['Search may be expired, try again in 5 seconds'];
              	pause(5);internet_wait=internet_wait+5;
            elseif  strcmp(theException.message,['Unknow status for search ' rid])
               	disp ['Unknow status for search ' rid];
                pause(5);internet_wait=internet_wait+5;
            else
                rethrow(theException)
%                 disp([theException.message, ', try in 5 seconds']);
%                 pause(5);internet_wait=internet_wait+5;
            end
        end
    end
    if internet_wait>waiting&&isempty(blastrestmp)
        error('myApp:argChk','time out!');
    end
    blastnumi=numel(blastrestmp);
    if blastnumi>0
        blastres(blastnum+1:blastnum+blastnumi) = blastrestmp(1:blastnumi);
        blastnum=numel(blastres);
    end
    dbnum=dbnum+1;
    tinterval=toc(t1);
end
if ~blastnum
    blastres=[];
end



