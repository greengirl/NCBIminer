function [RID RTOE] = blastncbi(seq, program, varargin)
%BLASTNCBI generates a remote BLAST request.
%
%   BLASTNCBI(SEQ, PROGRAM) sends a BLAST request against sequence SEQ to
%   NCBI using the specified PROGRAM. The basic local alignment search tool
%   (BLAST) offers a fast and powerful comparative analysis of interesting
%   protein and nucleotide sequences against known structures in existing
%   online databases. With no output arguments, BLASTNCBI returns a command
%   window link to the actual NCBI Report. A call with one output argument
%   returns the report ID (RID), while two output arguments return both the
%   RID and the request time of execution (RTOE) (an estimate of the time
%   until completion). SEQ can be a GenBank or RefSeq accession number, a
%   GI, a FASTA file, a URL, a string, a character array, or a MATLAB
%   structure that contains a sequence. PROGRAM can be 'blastn', 'blastp',
%   'psiblast', 'blastx', 'tblastn', 'tblastx' or 'megablast'. BLASTNCBI
%   uses the NCBI default values for the optional arguments: 'nr' for the
%   database, 'L' for the filter, and '10' for the expectation threshold.
%   The default values for the remaining optional arguments depend on which
%   program is used. For help in selecting a BLAST program, visit
%   https://www.ncbi.nlm.nih.gov/BLAST/producttable.shtml.
%
%   Information on all the optional parameters can be found at
%   https://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/blastcgihelp_new.html.
%
%   BLASTNCBI(SEQ, PROGRAM, ... , 'DATABASE', dbase) selects the database
%   for the alignment search. Compatible databases depend on the type of
%   sequence submitted. For nucleotide sequences, choices are 'nr', 'refseq_rna',
%   'refseq_genomic', 'est', 'est_human', 'est_mouse', 'est_others', 'gss',
%   'htgs', 'pat', 'pdb', 'month', 'alu_repeats', 'dbsts', 'chromosome',
%   'env_nt', or 'wgs'. For peptides, choices are 'nr', 'env_nr', 'refseq_protein',
%   'swissprot', 'pat', 'pdb', or 'month'. The nonredundant database, 'nr',
%   is the default database for both nucleotide and peptide sequences.
%
%   BLASTNCBI(SEQ, PROGRAM, ... , 'DESCRIPTIONS', quantity) sets the number
%   of short descriptions returned to the quantity specified. The default
%   limit is 100 descriptions. For PSI-BLAST default is 500.
%
%   BLASTNCBI(SEQ, PROGRAM, ... , 'ALIGNMENTS', quantity) specifies the
%   number of sequences for which high-scoring segment pairs (HSPs) are
%   reported. The default is 100 for all programs except PSI-BLAST, which is
%   500.
%
%   BLASTNCBI(SEQ, PROGRAM, ..., 'FILTER', filter) sets the filter to be
%   applied to the query sequence. The possible values are 'L' for low
%   complexity filter, 'R' for human repeats, 'm' for mask for lookup
%   table, and 'lcase' to turn on the lowercase mask. The default setting
%   is 'L' for a low-complexity
%   filter.
%
%   BLASTNCBI(SEQ, PROGRAM, ... , 'EXPECT', threshold) provides a
%   statistical significance threshold for matches against database
%   sequences. The default value for the threshold is 10. Learn more about
%   the statistics of local sequence comparison at
%   https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html#head2.
%
%   BLASTNCBI(SEQ, PROGRAM, ... , 'WORD', size) defines the word size. For
%   peptides, W can be 2 or 3 (3 is the default), and for nucleotides, 7,
%   11, or 15 (11 is the default). MEGABLAST accepts values of 11, 12,
%   16, 20, 24, 28, 32, 48, or 64, with a default of 28.
%
%   BLASTNCBI(SEQ, PROGRAM, ... , 'MATRIX', matrix) (protein sequences
%   only) specifies the substitution matrix that assigns the score for a
%   possible alignment of two amino acid residues. Allowed values include
%   'PAM30', 'PAM70', 'BLOSUM80', 'BLOSUM62', and 'BLOSUM45'. 'BLOSUM62' is
%   the default value.
%
%   BLASTNCBI(SEQ, PROGRAM, ... , 'GAPOPEN', penalty) (protein sequences
%   only) defines the penalty for opening a gap in the alignment.
%
%   BLASTNCBI(SEQ, PROGRAM, ... , 'EXTENDGAP', penalty) (protein sequences
%   only) defines the penalty for extending a gap in the alignment.
%
%   BLASTNCBI(SEQ, PROGRAM, ... , 'GAPCOSTS', penalty) (protein sequences
%   only) defines the penalties for opening and extending a gap in the
%   alignment. The penalty must be a two-element array.
%
%   Allowed values for gap costs vary by the actual substitution matrix
%   used. For the default 'BLOSUM62' substitution matrix, allowed values
%   are [9 2], [8 2], [7 2], [12 1], [11 1] (default), and [10 1]. Allowed
%   values for 'PAM30' are [7 2], [6 2], [5 2], [10 1], [9 1] (default),
%   [8 1]. Allowed values for 'PAM70' are [8 2], [7 2], [6 2], [11 1], 
%   [10 1] (default), [9 1].  Allowed values for 'BLOSUM80' are [8 2], [7 2], 
%   [6 2], [11 1], [10 1] (default), [9 1]. Allowed values for 'BLOSUM45' are
%   [13 3], [12 3], [11 3], [10 3], [15 2] (default), [14 2], [13 2], [12 2],
%   [19 1], [18 1], [17 1], [16 1]. 
%
%   BLASTNCBI(SEQ, PROGRAM, ... , 'INCLUSION', threshold) (PSI-BLAST only)
%   defines the statistical significance threshold for including a sequence
%   in the Position Specific Score Matrix (PSSM) created by PSI-BLAST for
%   the subsequent iteration. The default value is 0.005.
%
%   BLASTNCBI(SEQ, PROGRAM, ... , 'PCT', percent) (MEGABLAST only) is the
%   percent identity and the corresponding match and mismatch score for
%   matching existing sequences in public databases. Allowed values are
%   None, 99, 98, 95, 90, 85, 80, 75, and 60 (match and mismatch scores are
%   automatically applied). The default is 99 (99, 1, -3).
%   
%   BLASTNCBI(SEQ, PROGRAM, ..., 'ENTREZ', entrezquery) specifies the
%   Entrez query syntax to search a subset of the selected database. This
%   can be helpful to limit searches to molecule types, sequence lengths or
%   to exclude organisms. For more information about Entrez query syntax see
%   https://www.ncbi.nlm.nih.gov/books/bv.fcgi?rid=helpentrez.chapter.EntrezHelp.
%
%   Examples:
%
%       % Get a sequence structure from the Protein Data Bank.
%       S = getpdb('1CIV')
%
%       % Then use the structure as input for a PSI-BLAST search with an
%       % inclusion threshold of 1e-3.
%
%       blastncbi(S,'psiblast','inclusion',1e-3)
%
%       % Now click the URL link (NCBI BLAST Report) to view the first
%       % iteration of the report and to run subsequent searches.
%
%       % You can also perform a typical protein BLAST with an accession
%       % number and an alternative scoring matrix.
%       RID = blastncbi('AAA59174','blastp','matrix','PAM70','expect',1e-10)
%
%       % You can pass the RID to GETBLAST to download the report and parse it
%       % into a MATLAB structure.
%       getblast(RID)
%
%   See also BLASTFORMAT, BLASTLOCAL, BLASTREAD, BLASTREADLOCAL, GETBLAST.

%   Copyright 2003-2012 The MathWorks, Inc.


%   Reference: Altschul, Stephen F., Thomas L. Madden, Alejandro A.
%   Schaffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman
%   (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein
%   database search programs",  Nucleic Acids Res. 25:3389-3402.

if ~usejava('jvm')
    error(message('bioinfo:blastncbi:NeedJVM', mfilename));
end

%=== determine the format of the sequence
if isstruct(seq)
    seqtext = bioinfoprivate.seqfromstruct(seq);
elseif regexpi(seq,'^\d{6,8}$|^[a-zA-Z]{1,3}_?[a-zA-Z]{0,4}\d{5,8}(.\d)?$')
    seqtext = seq;  %true for a GI or a GenBank or RefSeq Accession number
elseif regexpi(deblank(seq),'^[a-zA-Z\*\-]+$')
    seqtext = deblank(seq);  %true for a sequence string
else
    try %try as a FASTA file
        fastaseq = fastaread(seq);
    catch theErr
        error(message('bioinfo:blastncbi:InvalidSeq'));
    end
    if size(fastaseq,1) > 1
        error(message('bioinfo:blastncbi:MultiSeq'));
    else
        seqtext = fastaseq.Sequence;
    end
end

%=== check sequence length
if length(seqtext) > 8000
    disp('Go to <a href="https://www.ncbi.nlm.nih.gov/BLAST/">NCBI</a> to BLAST sequences in excess of 8000 characters.');
    RID = 'Sequence exceeds maximum length.';
    RTOE = [];
    return;
else
    query = ['&QUERY=' seqtext];
end

%=== Set default values and acceptable parameters for a BLASTP search

% Request parameters
%site = 'https://www.ncbi.nlm.nih.gov/blast/';
site = 'https://blast.ncbi.nlm.nih.gov/'; % new URL as of May 2008
cbs = '&COMPOSITION_BASED_STATISTICS=yes';
dbase = '&DATABASE=nr';
cdd = '&CDD_SEARCH=off'; % changed default to 'off' due to stricter parameter checking at NCBI
filter = '&FILTER=L';
threshold = '&EXPECT=10';
word = '&WORD_SIZE=3';
matrix = '&MATRIX_NAME=BLOSUM62';
gapopen = '';
gapextend = '';
gapcosts = '';
gap = '&GAPCOSTS=11%201';
pct = '';
entrez = '';
doBlast = true;

% Report parameters
desc = '&DESCRIPTIONS=100';
algnmt = '&ALIGNMENTS=100';
ithresh = '';
format_out = '&NCBI_GI=on&SHOW_LINKOUT=on&SHOW_OVERVIEW=on';

% Acceptable options
okdb = {'nr','refseq_protein','swissprot','pat','pdb','month','env_nr'};
okfilter = {'L','m','lcase'};
okword = [2 3];
okmatrix = {'PAM30','PAM70','BLOSUM80','BLOSUM62','BLOSUM45'};
okgap = [9 2;8 2;7 2;12 1;11 1;10 1];
okpct = {};
okaaprog = {'blastp','blastx','psiblast','tblastn','tblastx'};
okntprog = {'blastn','megablast'};

if (ischar(program))
    p = strmatch(lower(program),[okaaprog okntprog]); 
    if isempty(p)
        error(message('bioinfo:blastncbi:BadProgram'));
    else
        program = ['&PROGRAM=' program];
        if (p == 3) % true if it is a psi-blast search
            program = '&PROGRAM=blastp&RUN_PSIBLAST=on&HITLIST_SIZE=500';
            desc = '&DESCRIPTIONS=500';
            algnmt = '&ALIGNMENTS=500';
            ithresh = '&I_THRESH=0.005';
            format_out = [format_out '&RUN_PSIBLAST=on&STEP_NUMBER=1'];
        end
        if (p > 3) % true if it is a translated or nucleotide search
            okdb = {'nr','est','est_human','est_mouse','est_others','gss','htgs',...
                'pat','pdb','month','alu_repeats','dbsts','chromosome','wgs',...
                'refseq_rna','refseq_genomic','env_nt'};
        end
        if (p >= 6) % true if it is a nucleotide search
            word = '&WORD_SIZE=11';
            matrix = '';
            cbs = '';
            cdd = '';
            gap = '';
            okfilter = {'L','m','lcase','R'};
            okmatrix = {};
            okgap = {};
            okword = [7 11 15];
        end
        if (p == 7) % true if it is a megablast search
            program = '&MEGABLAST=on&PROGRAM=blastn&PAGE=MegaBlast';
            word = '&WORD_SIZE=28';
            pct = '&PERC_IDENT=99,+1,+-3';
            okword = [11 12 16 20 24 28 32 48 64];
            okpct = {...
                'none,+1,+-2',...
                '99,+1,+-3',...
                '98,+1,+-3',...
                '95,+1,+-3',...
                '90,+1,+-2',...
                '85,+1,+-2',...
                '80,+2,+-3',...
                '75,+4,+-5',...
                '60,+1,+-1'};
        end
    end
end

%=== parse options
if nargin > 2
    if rem(nargin,2) == 1
        error(message('bioinfo:blastncbi:IncorrectNumberOfArguments', mfilename));
    end
    % noblast option just creates the URL but does not send to BLAST.
    okargs = {'database','descriptions','alignments','filter','expect',...
        'word','matrix','gapopen','extendgap','gapcosts','inclusion',...
        'pct','entrez','noblast'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:blastncbi:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:blastncbi:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % database
                    if (ischar(pval) && any(strcmpi(pval,okdb)))
                        dbase = ['&DATABASE=' pval];
                    else
                        error(message('bioinfo:blastncbi:BadDatabaseParam'));
                    end
                case 2  % descriptions
                    if isnumeric(pval)
                        desc = sprintf('&DESCRIPTIONS=%d',pval);
                    else
                        error(message('bioinfo:blastncbi:BadDesctParam'));
                    end
                case 3  % alignments
                    if isnumeric(pval)
                        algnmt = sprintf('&ALIGNMENTS=%d',pval);
                    else
                        error(message('bioinfo:blastncbi:BadAlignParam'));
                    end
                case 4  % filter
                    if (ischar(pval) && any(strcmpi(pval,okfilter)))
                        if strcmpi(pval,'lcase')
                            filter = '&LCASE_MASK=true';
                        else
                            filter = ['&FILTER=' pval];
                        end
                    else
                        error(message('bioinfo:blastncbi:BadFilterParam'));
                    end
                case 5  % expect
                    if isnumeric(pval)
                        threshold = sprintf('&EXPECT=%e',pval);
                    else
                        error(message('bioinfo:blastncbi:BadExpectParam'));
                    end
                case 6  % word
                    if (isnumeric(pval) && (ismember(pval,okword)))
                        word = sprintf('&WORD_SIZE=%d',pval);
                    else
                        error(message('bioinfo:blastncbi:BadWordParam'));
                    end
                case 7  % matrix
                    val = strmatch(lower(pval),lower(okmatrix)); 
                    if length(val) == 1
                        pval = okmatrix{val};
                        matrix = ['&MATRIX_NAME=' pval];
                        switch(val)
                            case 1 % PAM30
                                okgap = [7 2;6 2;5 2;10 1;9 1;8 1];
                                gap = '&GAPCOSTS=9%201';
                            case {2,3} % PAM70 or BLOSUM80
                                okgap = [8 2;7 2;6 2;11 1;10 1;9 1];
                                gap = '&GAPCOSTS=10%201';
                            case 5 % BLOSUM45
                                okgap = [13 3;12 3;11 3;10 3;15 2;14 2;13 2;12 2;19 1;18 1;17 1;16 1];
                                gap = '&GAPCOSTS=15%202';
                        end
                    else
                        error(message('bioinfo:blastncbi:BadMatrixParam'));
                    end
                case 8  % gapopen
                    gapopen = pval;
                case 9  % extendgap
                    gapextend = pval;
                case 10 % gap costs
                    if numel(pval) ~= 2 || ~isnumeric(pval)
                        error(message('bioinfo:blastncbi:BadGapCostsParam'));
                    else
                        gapcosts = pval;
                    end
                case 11  % inclusion
                    if isnumeric(pval)
                        ithresh = sprintf('&I_THRESH=%f',pval);
                    else
                        error(message('bioinfo:blastncbi:BadInclusionParam'));
                    end
                case 12  % percent identity
                    if isnumeric(pval)
                        pval = sprintf('%d',pval);
                    end
                    pct = strmatch(pval,okpct); 
                    if ~isempty(pct)
                        pct = ['&PERC_IDENT=' char(okpct(pct))];
                    else
                        error(message('bioinfo:blastncbi:BadPercentParam'));
                    end
                case 13 % entrez
                    pval = strrep(pval, ' ', '%20');
                    entrez = ['&ENTREZ_QUERY=' pval];
                case 14  % noblast
                    doBlast = ~bioinfoprivate.opttf(pval);
                    if isempty(doBlast)
                        error(message('bioinfo:blastncbi:InputOptionNotLogical', upper( char( okargs( k ) ) )));
                    end
            end
        end
    end
end

%=== if using gap, verify valid penalties
if ~isempty(gapcosts)
    val = strmatch(gapcosts,okgap);
    if ~isempty(val)
        gap = sprintf('&GAPCOSTS=%d%%20%d',okgap(val,1),okgap(val,2));
        if ~isempty(gapextend) || ~isempty(gapopen)
            warning(message('bioinfo:blastncbi:RedundantGapCostParam', okgap( val, 1 ), okgap( val, 2 ))); 
        end
    else
        error(message('bioinfo:blastncbi:BadGapParam'));
    end
else
    if ~isempty(gapopen)
        val =  strmatch(gapopen,okgap); 
        if ~isempty(val)
            gap = sprintf('&GAPCOSTS=%d%%20%d',okgap(val,1),okgap(val,2));
            if ~isempty(gapextend) && (gapextend ~= okgap(val,2))
                warning(message('bioinfo:blastncbi:BadExtendGapParam', okgap( val, 2 )));
            end
            
        else
            error(message('bioinfo:blastncbi:BadGapParam'));
        end
    else
        if ~isempty(gapextend)
            error(message('bioinfo:blastncbi:MissingGapOpenParam'));
        end
    end
end

%=== call NCBI BLAST, get RID and RTOE
blasturl = [site 'Blast.cgi?CMD=Put' program cbs dbase cdd filter threshold word matrix gap pct query algnmt desc ithresh entrez];

%=== call the URL
if doBlast
    s = urlread(blasturl);
    start = regexp(s,'QBlastInfoBegin\n'); % Find first line of the QBlast Info

    if isempty(start)
        RID = regexpi(s,'(?<=Request ID</td><td>\s*<b>)(\w+)', 'match'); % psi-blast searches
        RTOE = [];
        if isempty(RID) % search for error message
            start=regexp(s,'ERROR: ');
            if isempty(start)
                error(message('bioinfo:blastncbi:BlastSubmissionUnknown'));
            else
                [dummy, endOfFile] = regexp(s,'</form>.*?\n'); 
                blasterror = regexprep(s(start:endOfFile),'<.*>(\s)*','');
                error(message('bioinfo:blastncbi:BlastSubmission', blasterror));
            end

        end
    else
        [dummy, endOfFile] = regexp(s,'\nQBlastInfoEnd.*?\n'); 
        s = s(start+16:endOfFile-14);
        %RID = regexpi(s, '(?<=RID\s*=\s*)\w{11}', 'match', 'once'); % change in RID format April 16 2007
        %RID = regexpi(s, '(?<=RID\s*=\s*)\w+', 'match', 'once'); % change in RID format April 16 2008
        RID = regexpi(s, '(?<=RID\s*=\s*)\w+(?=\s*RTOE)', 'match', 'once'); 
        
        if isempty(RID)
            error(message('bioinfo:blastncbi:InvalidSearch'));
        end
            
        if nargout > 1
            RTOE = str2double(regexpi(s,'(?<=RTOE\s*=\s*)\d+', 'match', 'once'));
        end
    end

    
    if (nargout == 0)
        format_out = [algnmt desc ithresh format_out];
        retrieveurl = [site 'Blast.cgi?CMD=Get&RID=' char(RID) format_out];
        
        %    disp(['<a href="matlab: web(''' retrieveurl ''')">NCBI BLAST Report</a>']);
        disp(['<a href="matlab: web(''' retrieveurl ''', ''-browser'')">NCBI BLAST Report</a>']);
        clear RID;
    end

else
    RID = blasturl;
    RTOE = [];
end
