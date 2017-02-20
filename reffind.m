function [refseqs,n]=reffind(refseq,blastres,featuretype,featurename,refdist,timeout)
% Created by Xiaoting Xu-2014/06/21 (xiaotingxu@gmail.com)
% Modified 2014/06/30
% input, refseq: query sequence; blastres: BLAST results get from function
% blast_new; featuretype: gene feature type; featurename: gene names and
% synonyms; refdist: maximum base pair differences in match positions (either start or stop position)
% of two blast hits
refseqs=[];
n=0;
GI={blastres.Gi};
Qstart_stop=[blastres.q_start; blastres.q_end]';
[uniGi, cSstart_stop, cQstart_stop]=blastcombine(blastres);
 reflen=length(refseq.Sequence);
if length(GI)>1
    groupdist=sqrt(2)*reflen*refdist;
	disEu = pdist(Qstart_stop,'euclidean'); % Euclidean distance
	clustTreeEu = linkage(disEu,'complete');
% 	[h,nodes] = dendrogram(clustTreeEu);
% 	maximum difference is coverage*length bps\
	hidx = cluster(clustTreeEu,'criterion','distance','cutoff',groupdist);
%       fprintf('\n show group plot...\n');
% 	ptsymb={'r*','bo','gs','go','m*','mo','rs','ms','bs','g*','y*'};
% 	for i = 1:max(hidx)
%       clust = find(hidx==i);
%       scatter(Qstart_stop(clust,1),Qstart_stop(clust,2),ptsymb{i});
%       hold on
% 	end
% 	hold off
% 	xlabel('Qstart'); ylabel('Qstop');
% 	grid on

else
	hidx=1;
end
indecies=[];
seqlen=abs(cSstart_stop(:,1)-cSstart_stop(:,2));
finishgi=zeros(size(cSstart_stop,1),1);
display('reference-sequences-identification, please wait......');
ngroups=max(hidx);
Headers=cell(ngroups,1);
Sequences=cell(ngroups,1);
for i = 1:ngroups
  	clust = hidx==i;
	GIs=GI(clust);
    id=find(ismember(uniGi,GIs));
  	clustlen=seqlen(id);
    [clustlen, sortid] = sort(clustlen,'descend');
    GIs=uniGi(id(sortid));
   	clustSIndices=cSstart_stop(id(sortid),:);
    clustQIndices=cQstart_stop(id(sortid),:);
    for j=1:length(clustlen)
        interwait=0;
        finishid=strcmp(uniGi,GIs{j});
        if finishgi(finishid), continue;end
     	finishgi(finishid)=1;
        if j>20,break;end
        tryurl=1;
        while interwait<=timeout&&tryurl
            try
                indecies=featurefind(GIs{j},featuretype,featurename,timeout);
                tryurl=0;
            catch err
                if matchstart(err.message,'Time out')
                    display('Time out, turn to the next\n');
                   	indecies=[]; 
               	elseif strcmp(err.identifier,'MATLAB:urlread:ConnectionFailed\n')
                   	pause(5), interwait=interwait+5;
               	elseif matchstart(err.message,'When extracting sequences for the feature\n');
                   	indecies=[];
                else
                  	rethrow(err);
                end
            end
        end
        if ~isempty(indecies)
            indecies=sort(indecies);
             if any(cell2mat(strfind(featuretype,'misc')))
                Sstartend=sort(clustSIndices(j,1:2));
                Sremstart=Sstartend(1)-indecies(1);
                Sremend=indecies(2)-Sstartend(2);
                Qremstart=min(clustQIndices(j,1:2))-1;
                Qremend=reflen-max(clustQIndices(j,1:2));
                if Sremstart<=0
                    newstart=Sstartend(1);
                elseif Sremstart<=Qremstart
                    newstart=indecies(1);
                else
                    newstart=Sstartend(1)-Qremstart;
                end
                if Sremend<=0
                    newend=Sstartend(2);
                elseif Sremend<=Qremend
                    newend=indecies(2);
                else
                    newend=Sstartend(2)+Qremend;
                end
                indecies=[newstart,newend];
             end
             [Organism,GIj,Sequence]=getncbiseqs(GIs{j},timeout,[indecies clustSIndices(j,3)]);
             id2=strcmp(GIj,GIs{j});
             if sum(id2)
                GIj=GIj(id2);
                Organism=Organism(id2);
                Sequence=Sequence(id2);
             end
             if clustSIndices(j,3)==1
                Headers{i}=[Organism{1},'|GI:',GIj{1},':' num2str(indecies(1)),'-', num2str(indecies(2))];
             else
                Headers{i}=[Organism{1},'|GI:',GIj{1},':' num2str(indecies(2)),'-', num2str(indecies(1))];
             end
             Sequences{i}=Sequence{1};
             lenseq=numel(strfind(Sequence,'n'));
             if  lenseq<0.3*reflen&&range(indecies)<=1.5*reflen
                 Sequences{i}=Sequence{1};break;
             else
                 Sequences{i}=[];
                 Headers{i}=[];
             end
        end
    end
    if isempty(Headers{i})&&clustlen(1)>=0.3*reflen&&clustlen(1)<=1.5*reflen
        [Organism,GIj,Sequence]=getncbiseqs(GIs{1},timeout,clustSIndices(1,:));
        id2=strcmp(GIj,GIs{j});
       	if sum(id2)
            GIj=GIj(id2);
            Organism=Organism(id2);
            Sequence=Sequence(id2);
        end
        if clustSIndices(1,3)==1
          	Headers{i}=[Organism{1},'|',GIj{1},':' num2str(clustSIndices(1,1)),'-', num2str(clustSIndices(1,2))];
        else
        	Headers{i}=[Organism{1},'|',GIj{1},':' num2str(clustSIndices(1,2)),'-', num2str(clustSIndices(1,1))];
        end
        Sequences{i}=Sequence{1};
    end
end
id=~cellfun(@isempty,Headers);
if sum(id)
    Headers=Headers(id);
    Sequence=Sequences(id);
    [Header, id]=unique(Headers,'First');
    refseqs=struct('Header',Header,'Sequence',Sequence(id));
    id=strcmpi({refseqs.Sequence},refseq.Sequence); 
    if any(id)
        n=find(id);
    else
        [~,~,Sequence]=getncbiseqs(uniGi{1},timeout,cSstart_stop(1,:));
        if strcmpi(Sequence,refseq.Sequence)
            n=sum(~id)+1;
            refseqs(n)=refseq;
        end
    end
end
