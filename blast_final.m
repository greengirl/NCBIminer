function blastres = blast_final(refseqs,BlastProgram,Entrez,Alignments,ExpectValue,waiting,database)
% Created by Xiaoting Xu-2014/06/21 (xiaotingxu@gmail.com)
% Modified 2014/06/30
% Modified 2014/07/24, one more parameter: database
m=numel(refseqs);
blastres2=cell(m,1);
for i = 1:m
    tstart=tic;
    blastres2{i} = blast_new(refseqs(i), BlastProgram,Entrez,Alignments,ExpectValue,waiting,database);
    if toc(tstart)<=30&&i<m
       pt=30-toc(tstart);
       fprintf('two blasts interval less than 30s, please waiting for %ds',pt);
       pause(pt);
    end
end
idnull=cellfun(@isstruct,blastres2);
if sum(idnull)
    blastres=blastres2(idnull);
else
    blastres=[];
end 
end

