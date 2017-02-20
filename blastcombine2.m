function [Gi, Sstart_stop, Qstart_stop]=blastcombine2(blastres)
a=cellfun(@(x,y,z) strcat(x,'-',num2str(y),'-',num2str(z)),{blastres.Gi}, {blastres.s_start}, {blastres.s_end},'UniformOutput',0);
[~,id]=unique(a);
blastres=blastres(id);
SIndices=[cell2mat({blastres.s_start})', cell2mat({blastres.s_end})',...
	cell2mat({blastres.Strand})'];
maxlen=max(abs(SIndices(:,1)-SIndices(:,2)));
Gi_total={blastres.Gi};

[Gi_total, idx] = sort(Gi_total);
SIndices=SIndices(idx,:);
[Giidx, Gis]=grp2idx(Gi_total);
maxGiidx=max(Giidx);
uni_Gi=unique(Giidx);

num_gi=arrayfun(@(x) sum(Giidx==x), uni_Gi);
Gisingle=uni_Gi(num_gi==1);
Gisingle_id=ismember(Giidx,Gisingle);

Gi=cell(maxGiidx,1);
Sstart_stop=zeros(maxGiidx,3);


num_single=numel(Gisingle);
Gi(1:num_single)=Gi_total(Gisingle_id);
Sstart_stop(1:numel(Gisingle),:)=SIndices(Gisingle_id,:);

Gimultiple =find(num_gi>1); 
num_multiple=numel(Gimultiple);
if num_multiple+num_single~=maxGiidx
    error('Missing Gis');
end
Gimultiple_id=ismember(Giidx,Gimultiple);
Giidx=Giidx(Gimultiple_id);
uni_Gi=uni_Gi(Gimultiple);
Gi(num_single+1:maxGiidx)=Gis(Gimultiple);
SIndices=SIndices(Gimultiple_id,:);
QIndices=[[blastres.q_start]', [blastres.q_end]'];
QIndices=QIndices(idx,:);
Qstart_stop=zeros(maxGiidx,2);
Qstart_stop(1:numel(Gisingle),:)=QIndices(Gisingle_id,:);
QIndices=QIndices(Gimultiple_id,:);

for i=1:num_multiple
    idx= Giidx==uni_Gi(i); 
    SIndicesi=SIndices(idx,:);
   	seqlen=abs(SIndicesi(:,2)-SIndicesi(:,1));
    [~,id]=sort(seqlen,'descend');
    SIndicesi=SIndicesi(id,:);
    QIndicesi=QIndices(idx,:);
    QIndicesi=QIndicesi(id,:);
    minS=min(SIndicesi(1,1:2));
    maxS=max(SIndicesi(1,1:2));
	minQ=min(QIndicesi(1,1:2));
	maxQ=max(QIndicesi(1,1:2));
    j=1;
    id=SIndicesi(:,3)==SIndicesi(1,3)&abs(SIndicesi(:,1)-SIndicesi(:,2))>10;
    SIndicesi=SIndicesi(id,:);
    QIndicesi=QIndicesi(id,:);
    while j<=sum(id)-1
        Qstart=min(QIndicesi(j+1,1:2));
%         Qend=max(QIndicesi(j+1,1:2));
       	Sstart=min(SIndicesi(j+1,1:2));
      	Send=max(SIndicesi(j+1,1:2)); 
      	if (minS<=Sstart&&minQ<=Qstart&&(Sstart-maxS)<=100)||...
          	(minS>=Sstart&&minQ>=Qstart&&(minS-Send)<=100)    
           	minS=min(min(SIndicesi(1:j+1,1:2)));
           	maxS=max(max(SIndicesi(1:j+1,1:2)));
           	minQ=min(min(QIndicesi(1:j+1,1:2)));
         	maxQ=max(max(QIndicesi(1:j+1,1:2))); 
        else
          	minS=min(min(SIndicesi(1:j,1:2)));
          	maxS=max(max(SIndicesi(1:j,1:2)));
          	minQ=min(min(QIndicesi(1:j,1:2)));
           	maxQ=max(max(QIndicesi(1:j,1:2)));           
           	break;
        end
        j=j+1;
    end
    Sstart_stop(num_single+i,1:3)=[minS,maxS,SIndicesi(1,3)];
   	Qstart_stop(num_single+i,:)=[minQ,maxQ];
end
end

