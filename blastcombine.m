function [Gi, Sstart_stop, Qstart_stop]=blastcombine(blastres)
% combine BLAST segments
% input: BLAST statitical result in structrue format with field Accession, 
% s_start (subject start position), s_end (subject end position), Strand, 
% q_start (Query Start position), q_end (Query end position), 
% alignlen (Aligned length). 
% Created by Xiaoting Xu-2014/06/21 (xiaotingxu@gmail.com)
% Modified 2014/06/30

[Giidx, Gi]=grp2idx({blastres.Gi});
numgi=max(Giidx);
uniidx=unique(Giidx);
num_gi=arrayfun(@(x) sum(Giidx==x), uniidx);

Sstart_stop=zeros(numgi,3);
Qstart_stop=zeros(numgi,2);

id=num_gi==1;
if sum(id)
    idx=ismember(Giidx,find(id));
  	Sstart_stop(id,:)=[cell2mat({blastres(idx).s_start})', cell2mat({blastres(idx).s_end})',...
       	cell2mat({blastres(idx).Strand})'];
    Qstart_stop(id,:)=[cell2mat({blastres(idx).q_start})', cell2mat({blastres(idx).q_end})'];
    leftid=find(~id);
    m=length(leftid);
else
    leftid=1:numgi;
    m=numgi;
end

for i =1:m
    idi=leftid(i);
    id=find(Giidx==idi);
	alignleni=cell2mat({blastres(id).alignlen});
    [~,idx]=sort(alignleni,'descend');
 	id=id(idx);
   	SIndicesi=[cell2mat({blastres(id).s_start})', cell2mat({blastres(id).s_end})',...
    	cell2mat({blastres(id).Strand})'];
 	QIndicesi=[cell2mat({blastres(id).q_start})', cell2mat({blastres(id).q_end})'];
    strand=SIndicesi(1,3);
    id=SIndicesi(:,3)==strand;
    Sstart_stop(idi,1:3)=[min(min(SIndicesi(id,1:2))),max(max(SIndicesi(id,1:2))),SIndicesi(1,3)];
    Qstart_stop(idi,1:2)=[min(min(QIndicesi(id,1:2))),max(max(QIndicesi(id,1:2)))];
end



    


