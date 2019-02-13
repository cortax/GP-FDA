function M = srm_b(trajs,K,knots,order,Ops)
%SRM_B  SRM with transformation y=[x+b]
%
%   Model = SRM_B(Trajs,K,knots,order,[Options])
%    - Trajs : 'Trajs' structure (see HELP CCToolbox)
%    - K     : number of clusters to find
%    - knots : vector of knot values (pass in [] if you want the algorithm
%              to select the knots "automatically"; Use Options.Kn if
%              you would like to specify the number of knots); note that
%              the resulting augmented knot sequence may have a longer
%              length (this is normal).
%    - order : order of polynomial regression
%
%   DefOps = SRM_B('options) returns the default options.
%
%   Options (some fields)
%    .Interval : range of allowed translations (e.g., [-2 2]), you  must
%                provide this if knots is empty when this functions is called.
%    .KnotMult : (default .25) multiplied by the complete length of data 
%                interval to get the number of knots (minimum of 1)

% Scott Gaffney   8 May 2003
% Department of Information and Computer Science
% University of California, Irvine
%
% Changes
% -----------------

PROGNAME = 'srm_b';
METHOD   = PROGNAME;
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


%%% Handle Argument Processing
%%
n = nargin;
%
% Check for calling convention
%
% ('options')
%
if (n==1 & strcmp(trajs,'options'))
  M = DefaultOptions([]);
  return;
end
%
% (Trajs,K,knots,order,[Options])
if (n<4), error([PROGNAME,': incorrect number of arguments.']); end
M.K = K;
M.knots = knots;
M.order = order;
Ops = cexist('Ops',[]);
%%
%%% End Argument Processing


% Preprocessing
Ops = DefaultOptions(Ops);
M.Options = Ops;
M.zero  = Ops.zero;
M.method = METHOD;

% build data matrices
[Y,x,Seq] = trajs2seq(trajs,M.zero,M.Options.MinLen);
if (isempty(M.knots))
  if (isempty(Ops.Interval))
    errorbox([PROGNAME,': you must provide either ', ...
        '''knots'' or ''Options.Interval''.']);
    return;
  end
  M.knots = SelectKnots(M,x,Y,Seq);
end

%% Handle graphics output
if (~isfield(Ops,'MsgHnd'))
  Ops.MsgHnd = -1;
end
CreatedMsgBar=0;
if (isempty(Ops.MsgHnd))
  Ops.MsgHnd = msgbar([],'');
  CreatedMsgBar=1;
end
if (Ops.ShowGraphics & isempty(trajs))
  trajs = seq2cell(Y,Seq);
elseif (~Ops.ShowGraphics)
  clear trajs;
end



%***************************************************************************
%   Begin Main Function
%***************************************************************************

%% Define some stuff
NumIter = 0;
Lhood = zeros(1,Ops.IterLimit);
N = Seq(end)-1;

% Initialize the algorithm
M = InitE(M,x,Y,Seq);
M = InitM(M,x,Y,Seq);


%%%%%%%%%%%%%%%%%%% E-M Algorithm
while (1)
  NumIter = NumIter + 1;
  if(Ops.MsgHnd>=0)
    msgbar(Ops.MsgHnd,sprintf('%sIteration %d',Ops.MsgPrefix,NumIter));
  end
  if (Ops.ShowGraphics), M.Options = showmodel(M,trajs); end

  %%% E-Step
  M = Estep(M,x,Y,Seq);
  [Lhood(NumIter),M] = CalcLike(M,N);
  [DoStop,Ops] = StoppingCondition(Lhood,NumIter,Ops);
  if (DoStop), break; end

  %%% M-Step
  M = Mstep(M,x,Y,Seq);
end
%%%%%%%%%%%%%%%%%%% E-M Algorithm


M = permuteModel(M); 
M.Lhood = Lhood(1:NumIter);
M.NumPoints = prod(size(Y));
[trash, M.C] = max(M.Pik,[],2);
M.TrainLhood = M.Lhood(end);
M.TrainLhood_ppt = M.Lhood(end)./M.NumPoints;

% Calculate number of independent parameters
[P,K,D] = size(M.Mu);
M.NumIndParams = (K-1) + K*P*D + K*D + K;  % alpha, mu, sigma, s
if (M.Options.Sigma.Share==1)
  M.NumIndParams = M.NumIndParams - K*(D-1);
end

if (CreatedMsgBar)
  delete(Ops.MsgHnd);
end


%***************************************************************************
%   End Main Function
%***************************************************************************




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% postval 
%%
%%
function val = postval(b,x,Y,Mu,knots,order,s,sigma)
Xhat = bsplinebasis(knots,order,x-b);
val = sum(sum((Y-Xhat*Mu).^2)./sigma) + b^2/s;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estep 
%%
%%
function M = Estep(M,x,Y,Seq)
[P,D,K] = size(M.Mu);
[N,D] = size(Y);
n = length(Seq)-1;
fun = @postval;
SearchOps = M.Options.SearchOps;

% set M.Eb to the mode of the posterior and evaluate M.Vb at that point
for k=1:K
  Mu    = M.Mu(:,:,k);
  sigma = M.Sigma(k,:);
  s     = M.S(k);
  for i=1:n
    indx = Seq(i):Seq(i+1)-1;
    ni = length(indx);
    b0 = M.Eb(i,k) * M.Options.PropStart;
    pt0 = [b0];
    M.Eb(i,k) = fminsearch(fun,pt0,SearchOps, ...
      x(indx),Y(indx,:),Mu,M.knots,M.order,s,sigma);
    
    % calculate the expected information at b
    Xhat = bsplinebasis(M.knots,M.order,x(indx)-M.Eb(i,k));
    DXi  = bsplinebasis(M.knots,M.order,x(indx)-M.Eb(i,k),1);
    D2Xi = bsplinebasis(M.knots,M.order,x(indx)-M.Eb(i,k),2);
    DXiMu = DXi*Mu;
    D2XiMu = D2Xi*Mu;
    YxMu = Y(indx,:)-Xhat*Mu;
    Ib = sum(sum(YxMu.*D2XiMu - DXiMu.*DXiMu)./sigma) -1/s;
    M.Vb(i,k) = -1/Ib;
  end
end

M = CalcPik(M,x,Y,Seq);
% M.Pik = ones(n,1);  % for alignment
% M.scale = 1;        % for alignment



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CalcPik 
%%
function M = CalcPik(M,x,Y,Seq)
% Numerical integration
NumSamps = 100;
MaxTries = 5;
[N,D] = size(Y);
n = length(Seq)-1;
K = M.K;
P = M.order+1;
mlen = max(diff(Seq));
Piid = zeros(N,D);
Piimk = zeros(N,NumSamps,K);
S = M.S;
TotalSamps = 0;
tries = 1;

while (1)
  TotalSamps = TotalSamps + NumSamps;
  M.Pik(:) = 0;

  % calculate the density at sampled points
  for k=1:K
    b = randn(NumSamps,1).*sqrt(S(k));
    for j=1:NumSamps
      Xhat = bsplinebasis(M.knots,M.order,x-b(j));
      for d=1:D
        Piid(:,d) = normpdf(Y(:,d),Xhat*M.Mu(:,d,k),M.Sigma(k,d));
      end
      Piimk(:,j,k) = prod(Piid,2);
    end
  end

  % now scale the data to avoid underflow with long curves
  % and sum across the sample integration points
  M.scale = mean(mean(mean(Piimk)));
  Piimk_scl = Piimk./M.scale;  % we don't scale across D; we got rid of D above.
  for k=1:K
    for j=1:TotalSamps
      M.Pik(:,k) = M.Pik(:,k) + sprod(Piimk_scl(:,j,k),Seq,mlen);
    end
  end
  clear Piimk_scl;
  M.Pik = (M.Pik./TotalSamps) .* (ones(n,1)*M.Alpha');
  if (all(sum(M.Pik,2))), break; end  % check for stable integration

  % we have detected some zeros, try again?
  if (tries==MaxTries)
    fprintf(['srm_tt_sh: Integration failed, using realmin*1e100 ',...
      'instead.\n']);
    zero = find(sum(M.Pik,2)==0);
    M.Pik(zero,:) = realmin*1e100*(ones(length(zero),1)*M.Alpha');
    break;
  else
    fprintf(['srm_tt_sh: Zero membership detected, trying ', ...
        'integration again: %d\n'],tries);
    tries = tries+1;
    S = 1.25*S;  % biased, but gets over some tricky integrations
    Piimk = [zeros(N,NumSamps,K) Piimk]; % save current values
  end
end



  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CalcLike
%%
function [Lhood,M] = CalcLike(M,N)
K = M.K;
s = sum(M.Pik,2);
Lhood = sum(log(s)) + N*log(M.scale);
M.Pik = M.Pik ./ (s*ones(1,K));  % normalize the memberships






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% StoppingCondition
%%
function [dostop,Ops] = StoppingCondition(Lhood,NumIter,Ops)
dostop=0;

if (NumIter >= Ops.IterLimit)
  dostop = 1;
  return;
end
% return;  % for alignment

if (NumIter ~=1)
  if (isnan(Lhood(NumIter)))
    fprintf('the log-likelihood is equal to NaN.\n');
    dostop = 1;
  end
  if (Lhood(NumIter) < Lhood(NumIter-1))
    Ops.NumDec = Ops.NumDec+1;
    %fprintf(['the log-likelihood wabbled down', ...
    % ' on this iteration from ',num2str(Lhood(NumIter-1)),' to ', ...
    % num2str(Lhood(NumIter)),'.\n']);
    if (Ops.NumDec>=Ops.MaxDec)
      dostop = 1;
    end
  else
    abs_change = Lhood(NumIter)-Lhood(1);
    if (abs_change==0)
      dostop = 1;
    else
      delta = (Lhood(NumIter)-Lhood(NumIter-1)) / abs_change;
      if (abs(delta) < Ops.stopval)
        dostop = 1;
      end
    end
  end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mstep
%%
function M = Mstep(M,x,Y,Seq)
[N,D] = size(Y);
K = M.K;
P = length(M.knots)-M.order;
n = length(Seq)-1;
lens = diff(Seq)';

% Alpha
Wk = sum(M.Pik);
ff = find(Wk==0);
if (~isempty(ff))
  fprintf('Some Wk are zero\n');
  Wk(ff) = 0.01*sum(Wk);
end
M.Alpha = Wk'./n;

for k=1:K
  Piik = copy(M.Pik(:,k),lens);

  %% s
  M.S(k) = sum(M.Pik(:,k).*(M.Eb(:,k).^2 + M.Vb(:,k)))./Wk(k);

  % now make Xhat and finish up the M-Step
  Xhat = bsplinebasis(M.knots,M.order,x-copy(M.Eb(:,k),lens));
  PikXhat = Piik*ones(1,P).*Xhat;
  
  % check valid knot interval
  mass = sum(PikXhat);
  invalid = find(mass<1);
  valid = (1:P)';  valid(invalid)=[];
  
  % specialized error code: shouldn't be needed under normal circumstances
  if (isempty(valid))  % this only happens when there isn't enough data...
    fprintf('This cluster has no support.\n');  % ...and in this case, you...
    if (P<=2)          % ...shouldn't be fitting this model, but this...
      valid = [1:P]';  % ...section will get rid of crashing issues in any case
      invalid = [];
    else
      valid = [2:P-1]';  % often this will get over the problem, but the...
      invalid = [1 P]';  % ...resulting Mu won't be accurate
    end
  end
  
  % remove unwanted spline coefs (those with no support)
  PikXhat(:,invalid) = [];
  
  %% Mu
  M.Mu(:,:,k) = (PikXhat'*Xhat) \ (PikXhat'*Y);

  % fill-in invalid coefficients w/ closest valid ones
  f = find(invalid<valid(1));
  M.Mu(invalid(f),:,k) = ones(length(f),1) * M.Mu(valid(1),:,k);
  f = find(invalid>valid(end));
  M.Mu(invalid(f),:,k) = ones(length(f),1) * M.Mu(valid(end),:,k);

  %% sigma
  M.Sigma(k,:) = Piik'*((Y-Xhat*M.Mu(:,:,k)).^2)./ sum(lens.*M.Pik(:,k));
end
if (M.Options.Sigma.Share)
  M.Sigma = (sum(M.Sigma,2)./D)*ones(1,D);
end

% check for strange circumstances
fs = find(M.S<=0);
fsig = find(M.Sigma<=0);
if (~isempty([fs(:); fsig(:)]))
  fprintf('Some variances are zero\n');
  M.S(fs) = M.Options.minvar;
  M.Sigma(fsig) = M.Options.minvar;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% InitM
%%
function M = InitM(M,x,Y,Seq)
M = Mstep(M,x,Y,Seq);  % default init

if (M.Options.InitializeWithExamples)
[P,D,K] = size(M.Mu);
X = bsplinebasis(M.knots,M.order,x);
rnd = randperm(length(Seq)-1);
i=0;
for k=1:K
  i=i+1;
  indx = Seq(rnd(i)):Seq(rnd(i)+1)-1;
  x = X(indx,:);
  
  % check valid knot interval
  mass = sum(x);
  invalid = find(mass<1);
  x(:,invalid) = [];
  valid = (1:P)';  valid(invalid)=[];
  for d=1:D
    M.Mu(valid,d,k) = wls(x,Y(indx,d));
  end
  
  % fill-in invalid coefficients w/ closest valid ones
  fd = find(invalid<valid(1));   
  M.Mu(invalid(fd),:,k) = ones(length(fd),1)*M.Mu(valid(1),:,k);
  fd = find(invalid>valid(end)); 
  M.Mu(invalid(fd),:,k) = ones(length(fd),1)*M.Mu(valid(end),:,k);
end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% InitE
%%
function M = InitE(M,x,Y,Seq)
[N,D] = size(Y);
K = M.K;
n = length(Seq)-1;

% E-step vars
M.Eb  = zeros(n,K);
M.Vb  = zeros(n,K);

M.Eb  = randn(n,K)*1;
M.Vb  = rand(n,K)*2;
M.Pik = exprnd(.5,n,K);
M.Pik = M.Pik ./ (sum(M.Pik,2)*ones(1,K));

% set the default maximum translation interval
if (isempty(M.Options.Interval))
  min_start = min(x(Seq(1:end-1)));
  max_end = max(x(Seq(2:end)-1));
  M.Options.Interval(2) = min_start-M.knots(1);
  M.Options.Interval(1) = max_end-M.knots(end);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SelectKnots
%%
function knots = SelectKnots(M,x,Y,Seq)
Start = min(x(Seq(1:end-1)));
End = max(x(Seq(2:end)-1));
Start = Start - M.Options.Interval(2);  % adjust for translation interval
End   = End   - M.Options.Interval(1);
len = max(diff(Seq)) + (M.Options.Interval(2)-M.Options.Interval(1));

knotmult = .25;
if (isfield(M.Options,'KnotMult'))
  knotmult = M.Options.KnotMult;
end

% check for number of knots
if (isfield(M.Options,'Kn') & ~isempty(M.Options.Kn))
  Kn = M.Options.Kn;
else
  Kn = max(ceil(len*knotmult),1);
end
knots = linspace(Start,End,Kn);  % uniform spacing
knots = addendpts(knots,M.order);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper Functions
%%

function Ops = DefaultOptions(Ops);
Ops = SetFieldDef(Ops,'zero','nozero');
Ops = SetFieldDef(Ops,'stopval',1e-7);
Ops = SetFieldDef(Ops,'IterLimit',25);
Ops = SetFieldDef(Ops,'MaxDec',4);
Ops = SetFieldDef(Ops,'NumDec',0);
Ops = SetFieldDef(Ops,'NumEMStarts',1);
Ops = SetFieldDef(Ops,'MsgPrefix','');
Ops = SetFieldDef(Ops,'ShowGraphics',0);
Ops = SetFieldDef(Ops,'MinLen',[]);
Ops = SetFieldDef(Ops,'Interval',[]); % see InitE for default value
Ops = SetFieldDef(Ops,'Sigma',[]);
Ops.Sigma = SetFieldDef(Ops.Sigma,'Share',0);
Ops = SetFieldDef(Ops,'PropStart',0.8);
Ops = SetFieldDef(Ops,'InitializeWithExamples',1);
%
Ops.SearchOps = optimset('fminsearch');
Ops.SearchOps.Display = 'off';
Ops.SearchOps.MaxFunEvals = 200;
Ops.SearchOps.MaxIter = 200;

function Options = showmodel(M,trajs)
M.Mu = permute(M.Mu, [1 3 2]);  % make p-K-D
[trash, M.C] = max(M.Pik,[],2);
Options = viewmodel(M,trajs,M.Options);

function M = permuteModel(M)
M.Mu  = permute(M.Mu,[1 3 2]);
