function [newrefseq, blastres, n] = blastref(refseqs,BlastProgram,Entrez,Alignments,...
    ExpectValue,waiting,featuretype,featurename,refdist,database)
% Created by Xiaoting Xu-2014/06/21 (xiaotingxu@gmail.com)
% Modified 2014/06/30

m=numel(refseqs);
blastres=[];
newrefseq=[];
n=[];
i=1;
tinterval=30;
Entrez=[Entrez, 'AND 50:1000000[sequence length]'];
    while isempty(newrefseq)&&i<=m
        j=1; 
        refseq=refseqs(i);
        while j<=2
            if tinterval<30
                pt=30-tinterval;
                fprintf('two blasts interval less than 30s, please waiting for %ds\n',pt);
                pause(pt);
            end
            tstart = tic;
            blastres = blast_new(refseq, BlastProgram,Entrez,Alignments,ExpectValue,waiting,database);
            if isstruct(blastres)
                [newrefseq,n]=reffind(refseq,blastres,featuretype,featurename,refdist,waiting);
%             refseq=reffind_new_notcombine(refseqs(i),blastres,featuretype,featurename,refdist,waiting);
%                 newrefseq=reffind_new(refseqs(i),blastres,featuretype,featurename,refdist,waiting);
                refcount=numel(newrefseq);
            %iterated once
                if refcount
                    if n
                        break;
                    else
                        seqlen1=cellfun(@length, {newrefseq(:).Sequence});
                        [~,id1]=max(seqlen1);
                        refseq=newrefseq(id1);
                    end
                else
                    break;
                end
            end
            j=j+1;
            tinterval=toc(tstart);
        end
        i=i+1;
    end