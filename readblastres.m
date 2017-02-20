function [Header, Sequences, colltable]=readblastres(blastres2,lenlimit,location,timeout,extseq)
% modified by xiaoting xu 2014/06/25
blastnum=0;
for k=1:numel(blastres2)
    blastnumi=numel(blastres2{k});
    blastres(blastnum+1:blastnum+blastnumi) = blastres2{k};
    blastnum=numel(blastres);
end
Header=[];
Sequences=[];
colltable=[];
% [Accessions, Sindecies]=combineall(Accessions',Sstart_stops);
[GIs, Sindecies]=blastcombine2(blastres);
a=Sindecies(:,1)>Sindecies(:,2);
tmp=Sindecies(a,1);
Sindecies(a,1)=Sindecies(a,2);
Sindecies(a,2)=tmp;
seqlen=abs(Sindecies(:,1)-Sindecies(:,2));
id=(seqlen>lenlimit(1)&seqlen<=lenlimit(2));
GIs=GIs(id);
Sindecies=Sindecies(id,:);
a=Sindecies(:,3)==2;
if isinf(extseq(1))
    Sindecies(~a,1)=0;
    Sindecies(a,2)=0;
end
if isinf(extseq(2))
    Sindecies(a,1)=0;
    Sindecies(~a,2)=0;
end

Sindecies(a,1)=Sindecies(a,1)-extseq(2);
Sindecies(~a,1)=Sindecies(~a,1)-extseq(1);
Sindecies(a,2)=Sindecies(a,2)+extseq(1);
Sindecies(~a,2)=Sindecies(~a,2)+extseq(2);

id=Sindecies<0;
Sindecies(id)=1;
id=(Sindecies(:,2)-Sindecies(:,1))<=0;
Sindecies(id,:)=[];
GIs(id)=[];
m=numel(GIs);
if m
	Organism=cell(m,1);
    Accs=cell(m,1);
	Sequences=cell(m,1);
	startbp=Sindecies(:,1);
	endbp=Sindecies(:,2);
 	strand=Sindecies(:,3);
    display('read sequence from genbank, please wait');
    if location
        specimen_voucher=cell(m,1);
        country=cell(m,1);
        lat_lon=cell(m,1);
        collection_date=cell(m,1);
        collected_by=cell(m,1);
        identified_by=cell(m,1);
        for i=1:m
            [Organism{i},Accs{i},Sequences{i},specimen_voucher{i},country{i},lat_lon{i},collection_date{i},collected_by{i},...
                identified_by{i}]=getncbiseqs(GIs{i},timeout, [startbp(i) endbp(i) strand(i)]);
%          	Header{i}=[Organism{i} '|' Accessions{i} ':' num2str(startbp(i)) '-' num2str(endbp(i))];
        end
        % modified on 24th,July
        len=cellfun(@length,Accs);
        if sum(len)>m
            id1=find(len>1);
            for i=1:numel(id1)
                id2=strcmp(Accs{id1(i)},GIs(id1(i)));
                Organism{id1(i)}=Organism{id1(i)}(id2);
                Sequences{id1(i)}=Sequences{id1(i)}(id2);
                specimen_voucher{id1(i)}=specimen_voucher{id1(i)}(id2);
                lat_lon{id1(i)}=lat_lon{id1(i)}(id2);
                collection_date{id1(i)}=collection_date{id1(i)}(id2);
                collected_by{id1(i)}=collected_by{id1(i)}(id2);
                identified_by{id1(i)}=identified_by{id1(i)}(id2);
                country{id1(i)}=country{id1(i)}(id2);
            end
        end
        % modified on 24th,July
        colltable=[GIs,[Organism{:}]',[specimen_voucher{:}]',[country{:}]',...
            [lat_lon{:}]',[collection_date{:}]',[collected_by{:}]',[identified_by{:}]',num2cell(Sindecies)];
    else
        for i=1:m
            display(i);
            [Organism{i},Accs{i},Sequences{i}]=getncbiseqs(GIs{i},timeout,[startbp(i) endbp(i) strand(i)]);
%             Header{i}=[Organism{i} '|' Accessions{i} ':' num2str(startbp(i)) '-' num2str(endbp(i))];
        end
        % modified on 24th,July
        len=cellfun(@length,Accs);
        if sum(len)>m
            id1=find(len>1);
            for i=1:sum(id1)
                id2=strcmp(Accs{id1(i)},GIs(id1(i)));
                Organism{id1(i)}=Organism{id1(i)}(id2);
                Sequences{id1(i)}=Sequences{id1(i)}(id2);
            end
        end
        % modified on 24th,July
        colltable=[GIs,[Organism{:}]',num2cell(Sindecies)];
    end
    Header = cell(m,1);
    a=Sindecies(:,3)==2;
    Sequences=[Sequences{:}]';
    seqlen=cellfun(@length,Sequences);
    if sum(a)
        Header(a)=cellfun(@(x,y,z) [x, '|' y, ':', num2str(z(1)+z(2)-1) '-' num2str(z(1))],[Organism{a}]',GIs(a),...
            mat2cell([startbp(a),seqlen(a)],ones(sum(a),1),2),'uni',0);
    end
    if sum(~a)
        Header(~a)=cellfun(@(x,y,z) [x, '|' y, ':', num2str(z(1)) '-' num2str(z(1)+z(2)-1)],[Organism{~a}]',GIs(~a),...
            mat2cell([startbp(~a),seqlen(~a)],ones(sum(~a),1),2),'uni',0);
    end
end
