function [Xfit,Xbnd,normalization,chip] = pf0_fromexpresswiggle(fname, jmp, DNA_len, genome_file, minerr)
%% [Xfit,Xbnd,normalization,chip] = pf0(fname, jump, DNA_len, genome_file, minerr)
%%
%% If you use Grizzly Peak, please cite:
%% Melissa M. Harrison*, Xiao-Yong Li*, Tommy Kaplan*, Michael R. Botchan and Michael B. Eisen.
%% Zelda binding in the early Drosophila melanogaster embryo marks regions subsequently
%% activated at the maternal-to-zygotic transition
%% PLoS Genetics, 2011
%%
%% Distributed under the GNU General Public License (GPL)
%% Tommy Kaplan, tomkap@berkeley.edu


% [Xfit,Xbnd,normalization,chip] = pf('ZLD.cyc8.dm3.bed', 10, 150, 'genome_dm3.mat', 34)
if nargin<1, error('Usage: pf <.bed file> <jump [10]> <DNA_len [150]> <genome.mat file [genome_dm3.mat]> <minerr [30]>'); end;
if nargin<2, jmp=50; elseif ~isnumeric(jmp), jmp=str2double(jmp); end;
if nargin<3, DNA_len=250; elseif ~isnumeric(DNA_len), DNA_len=str2double(DNA_len); end;
if nargin<4, genome_file = 'genome_mm9.mat'; end;
if nargin<5, minerr=30; elseif ~isnumeric(minerr), minerr=str2double(minerr); end;
mname = [fname(1:end-4) '.mat'];
global lens; lens = ReadLensFromGenome(genome_file);
global jump; jump = jmp; clear jmp;

global clr; CLRs.dm3 = '255,92,92'; CLRs.dy13 = '55,192,55'; CLRs.dp26 = '127,165,255'; CLRs.dv12 = '192,88,255';
tkn=regexp(genome_file,'genome_(.*)\.mat','tokens');
if ~isempty(tkn) && isfield(CLRs,tkn{1}), clr = CLRs.(tkn{1}{1}); else, clr = '192,0,255'; end
clear tkn CLRs;

if ~exist(mname,'file')
    global chip ab
    [chip, ab] = CreateChip(fname,jump,DNA_len);
    save(mname,'chip','ab','jump','DNA_len');
else
    global chip ab
    load(mname)
end
% if ~exist(mname,'file'),
%     global chip ab; [chip, ab] = ReadReads(fname, jump, DNA_len); save(mname,'chip','ab','jump'); clear chip ab; %#ok<NASGU,REDEF>
% end
% global chip ab; load(mname); %#ok<REDEF>
% if ~exist([ab '.wig'],'file') && ~exist([ab '.wig.bz2'],'file'), dump_wig(ab); end
% %TF = upper(regexp(ab,'^[A-Za-z]*','match')); TF=TF{1};

% normalize to 100e6 % this actually normalizes to 100e6/jump, or in our case, 10e6.

chrs = fieldnames(chip); %#ok<*NODEF>
s=0; for i=1:length(chrs), s=s+sum(chip.(chrs{i})); end;
normalization=(1e7/jump*DNA_len)/s; for i=1:length(chrs), chip.(chrs{i}) = normalization*chip.(chrs{i}); end; %#ok<STRNU>
global KK; KK = DNA_len/jump;
% [c,l]=crosscorr(chipF.chr2R,chipR.chr2R); [~,K]=max(c); bar(jump*l,c); global KK; KK = l(K); clear c l K;
warning('off','MATLAB:rankDeficientMatrix');

% original shape
global shape lshape; shape = PeakShape(DNA_len,jump); shape = shape(shape>0); lshape = length(shape);
% remove flanks
global shp lshp; shp=shape(shape>0); shp = shp/max(shp); lshp = length(shp);
% extrapolate
global shape2 lshape2; shape2=shape(shape>0); shape2 = shape2/max(shape2); lshape2 = length(shape2);
shape2 = interp1(1:jump:jump*lshape2, shape2, 1:jump*lshape2,'linear','extrap'); lshape2 = length(shape2);

% initialize data
thrs = logspace(log10(50*minerr),log10(2*minerr),25); ext=minerr;
global Xdat Xfit Xbnd; Xdat = cell(0); Xfit = struct; Xbnd = cell(0); %#ok<REDEF,NASGU>
for t=1:length(chrs),
    chr = chrs{t};
    Xdat.(chr) = chip.(chr); nd = length(Xdat.(chr)); %#ok<NODEF>
    Xfit.(chr) = zeros(1, nd, 'single'); %#ok<STRNU>
end; clear y chr nd

% Load global genome
LoadGenome(genome_file); global genome; %#ok<NUSED>

kk=1;
for i=1:length(thrs),
    thr = thrs(i); fprintf('==> Threshold set at %.2f <==\n', thr);
    for t=1:length(chrs),
        chr = chrs{t}; nd=floor(lens.(chr)/jump); %length(Xdat.(chr));
        if ~any(Xdat.(chr)>thr), continue; end;
        % find segments above threshold
        ff = 2*(Xdat.(chr)>=thr)-1; % convert to 1/-1
        AA = conv(ff,[1 1 -1]); AA = AA(2:end-1); st = find(AA==3);  if AA(1)==2, st=[1,st]; end;
        AA = conv(ff,[1 -1 -1]); AA = AA(2:end-1); en = find(AA==-3); if AA(end)==-2, en=[en,length(AA)]; end;
        segs = [st;en]'; segs(segs<1)=1; segs(segs>nd)=nd;
        % iterate over putative peaks
        for s=1:size(segs,1)
            if ~any(Xdat.(chr)(segs(s,1):segs(s,2))>thr), continue; end;
            % extend
            XX1 = max(1,segs(s,1)-10*KK):segs(s,1); i1=find(Xdat.(chr)(XX1)<ext,1,'last');
            if isempty(i1),i1=XX1(1); else i1=XX1(i1); end;
            XX2 = segs(s,2):min(segs(s,2)+10*KK,nd); i2=find(Xdat.(chr)(XX2)<ext,1,'first');
            if isempty(i2),i2=XX2(end); else i2=XX2(i2); end;
            XXX = i1:i2; clear XX1 XX2;
            signal = Xdat.(chr)(XXX);
            % fit by peaks;
            [peaks,err,fitted,res]=pfr(signal,minerr,[]); %#ok<NASGU>
            % update tracks
            Xdat.(chr)(XXX) = ext-1; %Xdat.(chr)(XXX) - fitted; % override fitted region
            Xfit.(chr)(XXX) = Xfit.(chr)(XXX) + fitted; %#ok<STRNU>
            fprintf('%d %s:%d-%d %.1f | %d peaks (%.1f)\n',...
                kk, chr,jump*[i1 i2], max(signal), size(peaks,1), err);
            x=struct;
            x.chr = chr;
            x.pos = sprintf('%s:%d-%d',chr, jump*[i1 i2]);
            x.coords = jump*[i1 i2];
            % Peak information
            x.peaks = jump*XXX(peaks(:,1));
            x.heights = peaks(:,2)';
            x.err = err;

            Xbnd{kk}=x; %#ok<AGROW,NASGU>
            kk=kk+1;
        end
    end
end
%
warning('on','MATLAB:rankDeficientMatrix');
fname2 = [fname(1:end-4) '.fit.mat'];
save(fname2, 'ab','chip','Xfit','Xdat','Xbnd','jump','shp','DNA_len','normalization');
xname = [fname(1:end-4) '.xls'];
dump_xls(xname);
dump_peaks(ab);

%% ----------------------------------------------------------------------
function [peaks,Err,fitted,resi]=pfr(dat,thr,start)
yy = [zeros(1,100) double(dat) zeros(1,100)];
n = length(yy);
X = 101:n-100;
if nargin<2 || isempty(thr),thr=prctile(yy,50)/10; end;
if nargin<3, start = []; end;
param = start; [err,~,heights,resi] = fitshape(param,yy);
opts = optimset('TolX',.001,'Display','off','LargeScale','off');
% optimizition
while 1,
    improved=0;
    [~,res] = max(resi); respar = fminunc(@fitshape,res,opts,resi);
    tparam = fminsearch(@fitshape,[param respar],opts,yy);
    [terr,fitted,theights,tresi] = fitshape(tparam,yy);
    % accept new peak?
    if terr<err, err=terr; heights=theights; param=tparam; resi=tresi; improved=1; end;
    if ~improved || terr<thr, break; end;
end
% remove close-by peaks
[~,J]=sort(heights,'descend'); param=param(J); clear J;
P=[]; for i=1:length(param), if ~any(abs(param(i)-P)<=5), P(end+1)=param(i); end; end; param=P; clear i P; %#ok<AGROW>
param = fminsearch(@fitshape,param,opts,yy);

% remove sunken and internal and skewed peaks
[~,fitted,heights,~] = fitshape(param,yy); mh=max(heights);
rp = round(0.5+param);
PP=[rp'-1,rp'+1];
df = interp1(1:n,fitted,PP,'linear'); 
df=(df(:,2)-df(:,1))'/2;
for i=length(param):-1:1,
    if heights(i)<yy(rp(i))/2 || heights(i)<mh/10 || (heights(i)<mh/5 && abs(df(i))>5) || rp(i)<=100 || rp(i)>n-100
        param(i)=[]; heights(i)=[];
    end;
end;
clear i rp df;
param = fminsearch(@fitshape,param,opts,yy);
[Err,fitted,heights,resi] = fitshape(param,yy); % score all

% prepare to leave
fitted=fitted(X);
resi=resi(X);
param = round(0.5+param'-100); param(param<1)=1;n=n-200;param(param>n)=n;
peaks=[param,heights];

%% ----------------------------------------------------------------------
function [err,fitted,heights,residuals] = fitshape(lambda,y)
global shp; lshp=length(shp);
A = zeros(length(y),length(lambda));
for j = 1:length(lambda),
    f = round(lambda(j)-lshp/2+1);
    if f<1 || f+lshp-1>length(A), err=1e6; fitted=[]; heights=[]; residuals=[]; return; end;
    A(f:f+lshp-1,j)=shp;
end
heights = abs(A\y');
fitted = (A*heights)';
residuals = y-fitted;
% err = sqrt(mean(residuals.^2));
n=length(y); X=101:n-100;
err = mean(abs(residuals(X)));
%% ----------------------------------------------------------------------
function err = fun(beta,A,y,alpha,P)
fit = (A*beta)';
res = y-fit;
err = sqrt(mean(res.^2)) + alpha * (-log2(P) * beta);
% fprintf('%.2f ', beta); fprintf('=> %.2f\n',err);
%% ----------------------------------------------------------------------
function [chp, abt] = CreateChip(fname,jump,DNA_len)
i = regexp(fname, '\.bedgraph|\.wig'); abt = fname(1:i-1); clear i;
w=DNA_len; w=ceil(w/jump);

global lens
chrs = fieldnames(lens);
for t=1:length(chrs), chr = chrs{t}; L = lens.(chr); chp.(chr) = zeros(1,ceil(L/jump)+w,'single'); end

fid = fopen(fname,'r');
i = 0;
while 1 
    tline = fgetl(fid); if ~ischar(tline), break, end
    if i == 0
        i = i + 1;
        continue
    end
    
    % bedgraph format
    [chrom,coord,count] = strread(tline, '%s%d%f');
    
    crd = floor(coord/jump);
    % wig format
    %[chrom,coord1,coord2,count] = strread(tline,'%s%d%d%f');
	%crd = floor(coord1/jump):floor(coord2/jump)-1;
	
	chp.(char(chrom))(crd) = count;
end
fclose(fid);clear fid
%% ----------------------------------------------------------------------
function [shape] = PeakShape(DNA_len,jump)
KK = floor(DNA_len/jump);
mu=DNA_len;r=2;p=r/(mu+r);X=nbinrnd(r,p,1,1e6);X2=ceil(X.*rand(1,1e6));
I=(X>60&X<170);L1=ceil((1-X2(I))/jump);L2=floor((X(I)-X2(I))/jump);
H1=hist(L1,(-1000/jump):0);H2=hist(L2,0:(1000/jump));H1=H1/sum(H1);H2=H2/sum(H2);clear X X2 I L1 L2;
dat1=[H1 zeros(1,(1000/jump))];c1=conv(dat1,ones(1,KK));c1=c1(1:end-KK+1);
dat2=[zeros(1,(1000/jump)) H2];c2=conv(dat2,ones(1,KK));c2=c2(KK:end);
c=c1+c2; c=smooth(c,5)'; c=c/max(c);
% shape.X=-1000:jump:1000; shape.Y=c; shape.X0=ceil(length(shape.X)/2);
shape = c;
clear mu r p H1 H2 dat1 dat2 H1 H2 c c1 c2;
%% ----------------------------------------------------------------------
function [chip, ab] = ReadReads(fname,jump,DNA_len)
i = regexp(fname, '\.bed|\.bw'); ab = fname(1:i-1); clear i;
w=DNA_len; w=ceil(w/jump);
global lens;
chrs = fieldnames(lens);
for t=1:length(chrs), chr = chrs{t}; L = lens.(chr); chip.(chr) = zeros(1,ceil(L/jump)+w,'single'); end
[~,res] = system(['wc -l < ', fname] ); nl=str2double(res); clear res;
fid=fopen(fname,'r'); i=0;
nreads=0;
while 1 
    i=i+1;
    tline = fgetl(fid); if ~ischar(tline), break, end
    % bw format
    [p1,strand,chr,start,subseq,phred,score,mut] = strread(tline, '%s%c%s%d%s%s%f%s', 'delimiter', '\t');
    chr = chr{1};
    if strcmp(chr,'chr2L_Benjamin')
        chr = 'chr2L';
    end
    start = ceil(start/jump); endo = start;
    % bed format
    %[chr,start,endo,~,~,strand] = strread(tline, '%s%d%d%s%s%c', 'delimiter', '\t'); chr = chr{1};
    %start = ceil(start/jump); endo = ceil(endo/jump);
    if ~isfield(lens,chr), continue; end;
    if strand=='+', pp = start:start+w-1; else pp = endo-w+1:endo; end
    if min(pp)<1 || max(pp)>length(chip.(chr)), continue; end;
    chip.(chr)(pp) = chip.(chr)(pp) + 1; nreads = nreads + 1;
    if mod(i,500)==0, fprintf('\r%10d tags read (of %d, %.0f%%)', i, nl, 100*i/nl); end;
end
clear i p1 strand chr start subseq phred mut;
% normalize to 10Mb reads
ratio =  nreads / 1e7;
for t=1:length(chrs), chr = chrs{t}; chip.(chr) = chip.(chr) / ratio; end
fclose(fid);
fprintf('\r%10d tags read (of %d, %.0f%%)\n', nl, nl, 100);
%% ----------------------------------------------------------------------
function [lens] = ReadLensFromGenome(genome_file)
LoadGenome(genome_file); global genome;
chrs = fieldnames(genome)';
for i=1:length(chrs), 
    if ~strcmp(chrs{i},'fname'), 
        lens.(chrs{i})=length(genome.(chrs{i}).Sequence); 
    end 
end
%% ----------------------------------------------------------------------
function [] = LoadGenome(genome_file)
global genome;
if isempty(genome),
    if exist(genome_file,'file'),
	load(genome_file);
    else
        fname = [genome_file(1:end-4) '.seq'];
        fid=fopen(fname, 'r');
        while 1
            name = textscan(fid, '%[^\t]', 1); if isempty(name{1}), break; end;
            name = name{1}{1}(2:end); disp(name);
            seq = textscan(fid, '%[^\n]', 1, 'bufsize', 1e8); genome.(name) = seq{1}{1};
        end
        fclose(fid);
        save(genome_file, 'genome');
    end
end
%% ----------------------------------------------------------------------
function [] = dump_xls(xname)
global Xbnd;
fid = fopen(xname,'w');
fprintf(fid,'pos\tchr\tfrom\tto\tpeaks\theights\terr\n');
for i=1:length(Xbnd),
    X=Xbnd{i};
    fprintf(fid, '%s\t%s\t%d\t%d\t', X.pos, X.chr, X.coords);
    fprintf(fid, '%s\t', sprintf('%d,',X.peaks));
    fprintf(fid, '%s\t', sprintf('%.2f,',X.heights));
    fprintf(fid, '%.2f\n', X.err);
end
fclose(fid);
%% ----------------------------------------------------------------------
function [] = dump_wig(fname)
global chip clr jump;
if isempty(clr), clr='192,0,255'; end;
fid = fopen([fname '.wig'], 'w');
fprintf(fid, ...
    ['track type=bedGraph name="%s" description="%s" visibility=full color=%s autoScale=on '...
    'maxHeightPixels=100:24:21 graphType=bar alwaysZero=on yLineOnOff=on gridDefault=on windowingFunction=mean\n'],...
    fname, fname, clr);
% segment
for chr = (fieldnames(chip)')
    lastpos = NaN; lastval = NaN;
    chr = chr{1};
    X = jump:jump:jump*(length(chip.(chr))-17);
    dat = chip.(chr);
    for i=1:length(X)
       d = abs(dat(i)-lastval);
       if ~isfinite(d) | d>1e-1,
          if isfinite(lastval),
              fprintf(fid, '%s %d %d %.1f\n', chr, lastpos, X(i), lastval);
          end;
          lastval = dat(i); lastpos = X(i);
       end
   end
end
fclose(fid);
%% ----------------------------------------------------------------------
function [] = dump_peaks(fname)
global Xbnd clr;
if isempty(clr), clr='192,0,255'; end;

fid = fopen([fname '.peaks.bed'], 'w');
fprintf(fid, 'track type=bed name="%s peaks" description="%s peaks" visibility=dense itemRgb=On\n', fname, fname);

for i=1:length(Xbnd),
    x=Xbnd{i};
    % region
    fprintf(fid, '%s %d %d %d 1 . 0 0 %s\n', x.chr, x.coords, i, clr);
    % peaks
    for j=1:length(x.peaks),
        fprintf(fid, '%s %d %d %.1f 1 . 0 0 255,192,55\n', x.chr, x.peaks(j)-50, x.peaks(j)+50, x.heights(j));
    end
    % sites
    if isfield(x,'sites'),
        c=155-ceil(155*x.p');
        for j=1:length(x.sites),
            st = x.sites(j)-2; en = st + length(x.seq{j});
            fprintf(fid, '%s %d %d %.2f 1 . 0 0 %d,%d,%d\n', x.chr, st, en, x.p(j), c(j), c(j), c(j));
        end
    end
end
fclose(fid);
