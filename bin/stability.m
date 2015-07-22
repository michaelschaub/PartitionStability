function [S, N, VI, C] = stability(G, T, varargin)
%STABILITY    Graph partitioning optimizing stability with the Louvain
%             algorithm
%
%   [S, N, VI, C] = STABILITY(G, T) finds the optimal partitions of the 
%   graph G by optimizing the stability at each Markov time in vector T. G 
%   can either be the list of the edges in the graph (in the form [node i, 
%   node j, weight of link i-j; node k, node l, weight of link k-l;...] if 
%   the graph is weighted, or [node i, node j; node k, node l;...] if the 
%   graph is unweighted) or the adjacendy matrix of the graph. S, N, VI and 
%   C contain respectively the stability, the number of communities, the 
%   variation of information, and the optimal partition for each 
%   Markov time contained in T. If T is not specified, the modularity
%   (equivalent to stability for T=1) is calculated. Ideally, Markov time
%   should be sampled exponentially (e.g.: T = 10.^[-2:0.01:2]).
%
%   [S, N, VI, C] = STABILITY(G, T,'PARAM',VALUE) accepts one or more
%   comma-separated parameter name/value pairs. For a list of parameters 
%   and values, see "Parameter Options."
%
%    
%   Parameter Options:
% 
%        Parameter      Value                                 Default
%        ---------      -----                                 -------
%        L              Number of optimisations of the          100
%                       Louvain algorithm to be done at 
%                       each Markov time.
%
%        M              The top M partitions among the L        100
%                       given at each Markov time by the
%                       L louvain optimisations will be
%                       used to compute the variation of
%                       information.    
%
%        laplacian      Allows to choose which type of     'normalised'
%                       laplacian should be used to 
%                       calculate the stability. It can
%                       either be 'combinatorial', or 
%                       'normalised'.
%
%        directed       activate stability for directed         none
%                       graphs. Note that transition matrices
%                       are defined for left multiplications 
%                       here, i.e. A_ij is the link from i to j.
%
%        teleport_tau   teleportation probability               0.15
%                       (only active if directed == true)
%
%
%        noVI           Disables the calculation of the         none
%                       robustness of the partitions.
%                       Disabling this can significantly 
%                       speed up the calculations.
%
%        out            Enables saving step by step the         ''
%                       partitions found in a .mat file 
%                       located in the current folder, 
%                       along with the number of 
%                       communities (N), the value of 
%                       Stability (S), and the variation 
%                       of information (VI) for the 
%                       partitions obtained at each 
%                       Markov Time.
%
%       full            Enforces the calculation of the         none
%                       full stability.
%
%       linearised      Enforces the calculation of the         none
%                       linearised stability.
%
%       nocheck         Disables the checks for the             none
%                       encoding of the graph. This can 
%                       save computational time but can 
%                       also lead to serious errors if 
%                       the graph has not been properly
%                       encoded.
%
%       prec            Precision: defines a threshold for      1e-9  
%                       the range of weights allowed in  
%                       the laplacian exponential matrix 
%                       of the full stability.
%
%       plot            Plots the plots the results of          none
%                       the stability, number of 
%                       communities and variation of 
%                       information as a function of the
%                       Markov time.           
%
%       v               Verbose mode.                           none
%
%       p               Parallel mode.                          none
%
%       t               Output as text files:                   none
%                       Enables saving step by step the         
%                       partitions found at each Markov 
%                       time in a text file located in a 
%                       folder named 'Partitions', as 
%                       well as the outputs in a text 
%                       file output_'...'.stdout matlab 
%                       The option 'out' must be on.
%
%



% Unparsed default parameters
Graph = [];                                     % List of edges of the graph to be partitioned
Time = 1;                                       % Markov times at which the graph should be partitioned
flag_matlabpool = false;


%$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%
%$                                          $%
%$          Arguments parsing               $%
%$                                          $%
%$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%



[StabilityFunction, OutputFile, prefix, Full, Sanity, plotStability, verbose, TextOutput, PARAMS] = parseinput(length(varargin),G,varargin);



% Argument 1: G

if nargin > 0
    if size(G,1) == size(G,2) && ~issparse(G)
            G=sparse(G);
    end
    % Check if the graph is correctly encoded
    if Sanity
        G=check(G, verbose, PARAMS);
    end
    % If the full stability is to be computed, Graph should be the
    % adjacency matrix.
    if Full
        if size(G,1) ~= size(G,2)
            if size(G,2)==3
                Graph=sparse(G(:,1)+1,G(:,2)+1,G(:,3));
            elseif size(G,2)==2
                Graph=sparse(G(:,1)+1,G(:,2)+1,ones(length(G(:,1)),1));
            else
                error('Wrong size for G: G should be a graph saved either as a list of edges (size(G)=[N,3] if weighted, size(G)=[N,2] if unweighted) or as an adjacency matrix (size(G)=[N,N])');
            end
        else
            Graph = sparse(G);
        end
        PARAMS.NbNodes = size(Graph,2);
    % if the linearised stability is to be computed, Graph should be the
    % list of edges.
    else
        if size(G,1) == size(G,2)
            [rows,cols,vals] = find(G');
            if sum(vals)==length(vals)
                Graph=[cols-1, rows-1 ones(size(cols))];
            else
                Graph=[cols-1, rows-1, vals];
            end
            clear rows cols vals;
        elseif size(G,2)==3
            Graph=G;
        else
            error('Wrong size for G: G should be a graph saved either as a list of edges (size(G)=[N,3] if weighted, size(G)=[N,2] if unweighted) or as an adjacency matrix (size(G)=[N,N])');
        end
        PARAMS.NbNodes = max(Graph(:,1))+1;
    end
else
    error('Please provide at least the graph to be partitioned. Type "help stability" for more information.');
end    

% Argument 2: T

if nargin > 1
    if (isvector(T) && isnumeric(T))
        Time=T;
    else
        error('The second argument should be a numerical vector. Type "help stability" for more information.');
    end
end

% Parallel computation: Initialize the number of cores if matlabpool is not
% yet running.
if PARAMS.ComputeParallel && (matlabpool('size') == 0)
    flag_matlabpool = true;
    matlabpool
end

%$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%
%$                                          $%
%$      Computation of the stability        $%
%$                                          $%
%$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%

if verbose
    c = clock;
    disp(' ');
    disp(['   Stability will be computed with the following parameters:']);
    disp(' ');
    if Full; disp('      - Full Stability'); else disp('      - Linearised Stability'); end
    if strfind(func2str(StabilityFunction), 'N'); disp(['      - Normalised laplacian']); else disp(['      - Combinatorial laplacian']); end
    if PARAMS.directed; disp(['      - DIRECTED graph; teleportation parameter set to ' num2str(PARAMS.teleport_tau)]); end
    if PARAMS.ComputeVI; disp(['      - Computation of the variation of information (VI): Yes']); else disp(['      - Computation of the variation of information: No']); end
    if OutputFile; disp(['      - Save the results in a file: Yes ' prefix]); else disp(['      - Save the results in a file: No']); end
    if OutputFile; disp(['      - Saved file prefix: ' prefix]); end
    if Sanity; disp(['      - Check the input graph: Yes']); else disp(['      - Check the input graph: No']); end
    if plotStability; disp(['      - Plot the results: Yes']); else disp(['      - Plot the results: No']); end
    if verbose; disp(['      - Verbose mode: Yes']);    else disp(['      - Verbose mode: No']);    end
    if PARAMS.ComputeParallel; disp(['      - Parallel computation: Yes']); else disp(['      - Parallel computation: No']); end
    disp(['      - Number of Louvain iterations: ' int2str(PARAMS.NbLouvain)]);
    disp(['      - Number of Louvain iterations used for the computation of VI: ' int2str(PARAMS.M)]);
    disp(['      - Precision used: ' num2str(PARAMS.Precision)]); 
    disp(' ');
    tstart=tic;
end


% Initialisation
S = zeros(1, length(Time));
N = zeros(1, length(Time));
VI = zeros(1, length(Time));
C = zeros(PARAMS.NbNodes, length(Time));

if TextOutput
    mkdir(['Partitions_' prefix]);
end

if OutputFile
    save(['Stability_' prefix '.mat'],'Time','S','N','VI','C');
end

if plotStability
    figure_handle = figure;
end

if verbose
    step_prec=0;
end


% Loop over all Markov times
for t=1:length(Time)

    if verbose
        disp(['   Partitioning for Markov time = ' num2str(Time(t),'%10.6f') '...']);
    end
    
    
    
    [S(t), N(t), C(:,t), VI(t) VAROUT] = StabilityFunction(Graph, Time(t), PARAMS);
    if isfield(VAROUT,'precomputed')
        PARAMS.precomputed = VAROUT.precomputed;
        PARAMS.pi = VAROUT.pi;
        PARAMS.P = VAROUT.P;
    end
    
    if plotStability && t>1
        stability_plot(Time,t,S,N,VI,PARAMS.ComputeVI,figure_handle);
    end
    
    if TextOutput
        cd(['Partitions_' prefix]);
        dlmwrite(['Partition_' prefix '_' num2str(Time(t),'%10.6f') '.dat'],[[1:PARAMS.NbNodes]',C(:,t)],'delimiter','\t');
        cd ..;        
        dlmwrite(['Stability_' prefix '.stdout'],[Time(t), S(t), N(t), VI(t)],'-append', 'delimiter','\t')
    end   
    
    if OutputFile
        save(['Stability_' prefix '.mat'],'Time','S','N','VI','C','-append');
    end

    
    if verbose && 100*t/length(Time) >= step_prec+10        
        disp(' ');
        disp(['   Completed: ' num2str(round(100*t/length(Time)),10) '%']);
        remaining_time=toc(tstart)*(1-t/length(Time))/(t/length(Time));
        nb_hours = floor(remaining_time/3600);
        nb_min = floor((remaining_time - nb_hours*3600)/60);
        nb_sec = round(remaining_time - nb_hours*3600 - nb_min*60);
        disp(['   Estimated time remaining: ' datestr([2011  1 1 nb_hours nb_min nb_sec], 'HH:MM:SS')]);%num2str(nb_hours) ':' num2str(nb_min) ':' num2str(nb_sec)]);
        disp(' ');
        step_prec = floor(100*t/length(Time));
    end

end

if verbose
    c = clock;
    disp(' ');
    disp(['   Partitioning of the graph finished at ' datestr([2011 1 1 c(4) c(5) c(6)], 'HH:MM:SS')]);
    remaining_time=toc(tstart);
    nb_hours = floor(remaining_time/3600);
    nb_min = floor((remaining_time - nb_hours*3600)/60);
    nb_sec = round(remaining_time - nb_hours*3600 - nb_min*60);
    disp(['   Total time needed: ' datestr([2011 1 1 nb_hours nb_min nb_sec], 'HH:MM:SS')]);%num2str(nb_hours) ':' num2str(nb_min) ':' num2str(nb_sec)]);
end


if OutputFile
    save(['Stability_' prefix '.mat'],'Time','S','N','VI','C','-append');
end

if flag_matlabpool
    matlabpool close;
end



end

%------------------------------------------------------------------------------
function [StabilityFunction, OutputFile, prefix, Full, Sanity, plotStability, verbose, TextOutput,PARAMS] = parseinput(options,G,varargin)
% Parse the options from the command line

%$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%
%$                                          $%
%$          Default parameters              $%
%$                                          $%
%$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Global" options relevant for output and control flow

stability_type_specified=false;
threshold_nnodes_full = 1000;
threshold_nedges_full = 5000;

OutputFile = false;                             % No output file by default.

Laplacian = 'Normalised';                       % Default Laplacian
Full = false;                                    % If true, performs the full stability
Sanity = true;                                  % If true, performs the graph sanity checks
plotStability = false;                          % If true, plots the results of the stability, number of communities and variation of information vs Markov time.           
verbose = false;                                % Toggle verbose mode
prefix = '';                                    % Output prefix
TextOutput = false;                             % Toggles the text output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options stored in struct relevant for optimization etc.
PARAMS = struct;                                % create empty structure for storing parameters
PARAMS.precomputed = false;                     % Flag for precomputed transition matrix + stationary distribution
PARAMS.directed = false;                        % enables dealing with directed graphs
PARAMS.ComputeVI = true;                        % True if the variation of information should be computed
PARAMS.ComputeES = false;                       % True if edge statistics should be computed
PARAMS.ComputeParallel = false;                 % Toggles the computation in parallel
PARAMS.NbLouvain = 100;                         % Number of louvain optimisations at each Markov time
PARAMS.NbNodes = 0;                             % Total number of nodes;
PARAMS.Precision = 1e-9;                        % Threshold for stability and edges weigths
PARAMS.M = 100;                                 % Top M partitions among the L found by louvain are used to compute the variation of information
PARAMS.K = NaN;                                 % K stabilities value, only relevant for Ruelle random walk and k stabilities. K =-1 corresponds to the normalised Laplacian
PARAMS.teleport_tau = 0.15;                     % teleportation probability (only relevant for directed graphs)

% actual parsing begins here
attributes={'novi', 'l', 'm', 'out', 'full','linearised', 'nocheck', 'laplacian', 'prec', 'plot','v','t','p','k','directed','teleport','precomputed'};

if options > 0
    
    varargin = varargin{:}; % extract cell array input from varargin
    
    % test whether attribute-value pairs are specified, or fixed parameter order
    stringoptions = lower(varargin(cellfun('isclass',varargin,'char')));
    attributeindexesinoptionlist = ismember(stringoptions,attributes);
    newinputform = any(attributeindexesinoptionlist);
    if newinputform
        % parse values to functions parameters
        i = 1;
        while (i <= length(varargin))
            if strcmpi(varargin{i},'full')
                stability_type_specified=true;
                Full = true;
                i = i+1;
            elseif strcmpi(varargin{i},'linearised')
                if stability_type_specified
                    warning('The program can run either the linearised or the full stability, not both simultaneously. Please choose only one of them. Linearised stability will be used here.');
                end
                stability_type_specified=true;
                Full = false;
                i = i+1;
            elseif strcmpi(varargin{i},'directed')
                PARAMS.directed = true;
                i = i+1;
            elseif strcmpi(varargin{i},'novi')
                PARAMS.ComputeVI = false;
                i = i+1;
            elseif strcmpi(varargin{i},'nocheck')
                Sanity = false;
                i = i+1;
            elseif strcmpi(varargin{i},'plot')
                plotStability = true;
                i = i+1;
            elseif strcmpi(varargin{i},'v')
                verbose = true;
                i = i+1;
            elseif strcmpi(varargin{i},'precomputed')
                PARAMS.precomputed = true;
                i = i+1;
            elseif strcmpi(varargin{i},'p')
                if exist('matlabpool','file')
                    PARAMS.ComputeParallel = true;
                else
                    PARAMS.ComputeParallel = false;
                    warning('The Parallel Computing Toolbox of Matlab does not appear to be installed. Defaulting to single node computation...');
                end
                i = i+1;
            elseif strcmpi(varargin{i},'t')
                TextOutput = true;
                i = i+1;
            else
                %Check to make sure that there is a pair to go with
                %this argument.
                if length(varargin) < i + 1
                    error('MATLAB:stability:AttributeList', ...
                        'Attribute %s requires a matching value', varargin{i});
                elseif strcmpi(varargin{i},'laplacian')
                    if ischar(varargin{i+1})
                        Laplacian = varargin{i+1};
                    else
                        error('MATLAB:stability:laplacian',...
                            'Please provide a matching value for attribute laplacian. It must either be ''normalised'' or ''combinatorial''.');
                    end
                elseif strcmpi(varargin{i},'l')
                    if isnumeric(varargin{i+1})
                        PARAMS.NbLouvain = round(varargin{i+1});
                        PARAMS.M = round(varargin{i+1});
                    end
                elseif strcmpi(varargin{i},'prec')
                    if isnumeric(varargin{i+1})
                        PARAMS.Precision = varargin{i+1};
                    end
                elseif strcmpi(varargin{i},'m')
                    if isnumeric(varargin{i+1})
                        PARAMS.M = varargin{i+1};
                    end
                elseif strcmpi(varargin{i},'k')
                    if isnumeric(varargin{i+1})
                        PARAMS.K = varargin{i+1};
                    end
                elseif strcmpi(varargin{i},'teleport')
                    if isnumeric(varargin{i+1})
                        PARAMS.teleport_tau = varargin{i+1};
                    end
                elseif strcmpi(varargin{i},'out')
                    if ischar(varargin{i+1})
                        OutputFile = true;
                        prefix = varargin{i+1};
                    else
                        error('MATLAB:stability:out',...
                            'Please provide a matching value for attribute out. It must be a string.');
                    end
                else
                    error('MATLAB:stability:Attribute',...
                        'Invalid attribute tag: %s', varargin{i});
                end
                i = i+2;
            end
        end
    else 
        if ischar(varargin{1})
            error('MATLAB:stability:Attribute',...
                            'Invalid attribute tag: %s', varargin{1});
        else
            error('MATLAB:stability:Attribute',...
                            'Invalid attribute tag: %d', varargin{1});
        end
    end
end

TextOutput = TextOutput & OutputFile;

if ~stability_type_specified
    if (size(G,1) == size(G,2) && nnz(G)<threshold_nedges_full && size(G,1)<threshold_nnodes_full) || (size(G,1) ~= size(G,2) && size(G,1)<threshold_nedges_full && max(max(G(:,1:2)))<threshold_nnodes_full)
        Full = true;
        disp(' ');
        disp('  ---------------------------------------------------------------');
        disp(['  You have not specified whether the full or the linearised      ';...
              '  version of Stability should be computed. The full Stability    ';...
              '  will be used here based on the size of your graph. If you want ';...
              '  to use the linearised Stability, please add ''linearised''       ';...
              '  in the arguments of stability.                                 ']);
        disp('  ---------------------------------------------------------------');

    else
        Full = false;
        disp(' ');
        disp('  ------------------------------------------------------------------');
        disp(['  You have not specified whether the full or the linearised        ';...
              '  version of Stability should be computed. The linearised Stability';...
              '  will be used here based on the size of your graph. If you want   ';...
              '  to use the full Stability, please add ''full'' in the arguments    ';...
              '  of stability.                                                    ']);
        disp('  ------------------------------------------------------------------');
    end
 
end

% Choose which type of stability is to be computed
if Full
    if strcmpi(Laplacian, 'combinatorial')
        StabilityFunction = @louvain_FCL;
    elseif strcmpi(Laplacian, 'normalised')
        StabilityFunction = @louvain_FNL;
    else
        error('Please provide a valid matching value for attribute laplacian. It must either be ''normalised'' or ''combinatorial''.');
    end
else
    if strcmpi(Laplacian, 'combinatorial')
        StabilityFunction = @louvain_LCL;
    elseif strcmpi(Laplacian, 'normalised')
        StabilityFunction = @louvain_LNL;
    else
        error('Please provide a valid matching value for attribute laplacian. It must either be ''normalised'' or ''combinatorial''.');
    end
end

end

%--------------------------------------------------------------------------
function shares = split_even(N,nr_threads)
%Function to compute an even split of N runs between t threads 
shares = ones(1,nr_threads)*floor(N/nr_threads);
shares(1) = shares(1)+ rem(N,nr_threads);
end

%------------------------------------------------------------------------------
function [S, N, C, VI, VAROUT] = louvain_FNL(Graph, time, PARAMS)
% Computes the full normalised stabilty

VAROUT =[]; % init varying outputs 

% "transition matrix" and pi unknown
if PARAMS.precomputed == false
    % directed case: M_ij >> from i to j
    if PARAMS.directed == true
        dout = sum(Graph,2);
        dangling = (dout==0);
        dout(dangling) = 1;
        Dout = sparse(diag(dout));
        clear dout;
        M = (1-PARAMS.teleport_tau)*(Dout\Graph); % deterministic part of transition
        % teleportation according to arXiv:0812.1770
        M =	M + diag(PARAMS.teleport_tau + dangling.*(1-PARAMS.teleport_tau))...
            * ones(PARAMS.NbNodes)/PARAMS.NbNodes;
        
        clear Dout dangling
        [v lambda_all] = eigs(M'); % largest eigenvalue of transition matrix corresponds to stat.distribution.
        lambda = max(diag(lambda_all));
        v = v(:,diag(lambda_all) == lambda);
        v = abs(v);              % make sure eigenvector is positive
        clear lambda;
        % store results for future use
        VAROUT.precomputed = true;
        VAROUT.pi = v/sum(v);
        VAROUT.P = M;      
        % now compute exponential transition matrix
        solution = diag(v/sum(v))*expm(time* (M-eye(size(M))) );
        clear M v;
        % symmetrize solution
        solution = (solution+solution')/2;
        
        
        % undirected case
    else
        % Generate the matrix exponential
        trans=sparse(diag(    (sum(Graph)).^(-1)     ) * Graph);  %(stochastic) transition matrix        
        Lap=sparse(trans-eye(PARAMS.NbNodes));
        % store results for future use
        VAROUT.precomputed = true;
        VAROUT.P = trans;
        
        clear trans;
        exponential=expm(time.*Lap);
        clear Lap;
        
        PI=sparse((diag(sum(Graph)))/sum(sum(Graph)));  %diag matrix with stat distr   
        VAROUT.pi = diag(PI);   % store results for future use
        
        solution=sparse(PI*exponential);
        clear exponential;
        clear PI;
        
        
    end
    
    % stationary distribution and "transition matrix" have been computed before
else
    if PARAMS.directed == true
        solution = diag(PARAMS.pi)*expm(time* (PARAMS.P - eye(size(PARAMS.P))) );
        solution = (solution +solution')/2; % symetrization needed for directed case
    else
        solution = diag(PARAMS.pi)*expm(time* (PARAMS.P - eye(size(PARAMS.P))) );
    end
end



% prune out weights that are too small as defined by precision
solution=max(max(solution))*PARAMS.Precision*round(solution/(max(max(solution))*PARAMS.Precision));
[row,col,val] = find(solution);
clear solution
graph=[col-1,row-1,val];
weighted = 'w';

% Optimize louvain NbLouvain times
lnk = zeros(PARAMS.NbNodes, PARAMS.NbLouvain);
lnkS = zeros(PARAMS.NbLouvain,1);
precision = PARAMS.Precision;
if PARAMS.ComputeParallel
    parfor l=1:PARAMS.NbLouvain
        [stability, nb_comm, communities] = stability_louvain_LNL(graph, 1, precision, weighted);
        lnk(:,l) = communities;
        lnkS(l) = stability;
    end
else
    for l=1:PARAMS.NbLouvain
        [stability, nb_comm, communities] = stability_louvain_LNL(graph, 1, precision, weighted);
        lnk(:,l) = communities;
        lnkS(l) = stability;
    end
end
[S,indexS]=max(lnkS);
C=lnk(:,indexS);
N=max(C)+1;

clear communities;
clear graph;

ComputeVI = PARAMS.ComputeVI;
NbLouvain = PARAMS.NbLouvain;
NbNodes = PARAMS.NbNodes;
M = PARAMS.M ;

if ComputeVI && nnz(max(lnk)==NbNodes-1)~=NbLouvain && nnz(max(lnk)==0)~=NbLouvain
     VI = computeRobustness(lnk, lnkS, M, PARAMS.ComputeParallel);
else
    VI = 0;
end

clear lnk;

end

%------------------------------------------------------------------------------
function [S, N, C, VI, VAROUT] = louvain_FCL(Graph, time, PARAMS)
% Computes the full combinatorial stability

VAROUT =[]; % init varying outputs 

ComputeVI = PARAMS.ComputeVI;
precision = PARAMS.Precision;
NbLouvain = PARAMS.NbLouvain;



% "transition matrix" and pi unknown
if PARAMS.precomputed == false
    % directed case: M_ij >> from i to j
    if PARAMS.directed == true
        % Directed part so far not implemented, different variants are possible
        error(['This parameter combination is not allowed! If you want to use directed stability, use the normalised Laplacian']);        
        % undirected case
    else
        % standard Laplacian and average degree
        av_degree = sum(sum(Graph))/PARAMS.NbNodes;
        Lap=  -(Graph-diag(sum(Graph)));
        % store results for future use
        VAROUT.precomputed = true;
        VAROUT.P = Lap/av_degree;
        
 
        % Generate the matrix exponential
        exponential=sparse(expm(-time.*Lap/av_degree));
        clear Lap;
        
        PI=sparse(eye(PARAMS.NbNodes)/PARAMS.NbNodes);  %diag matrix with stat distr   
        VAROUT.pi = diag(PI);   % store results for future use
        
        solution=sparse(PI*exponential);
        clear exponential;
        clear PI;
        
        
    end
    
    % stationary distribution and "transition matrix" have been computed before
else
    if PARAMS.directed == true
        solution = diag(PARAMS.pi)*expm(-time* PARAMS.P);
        solution = (solution +solution')/2; % symetrization needed for directed case
    else
        solution = diag(PARAMS.pi)*expm(-time* PARAMS.P);
    end
end



% prune out weights that are too small as defined by precision
solution=max(max(solution))*PARAMS.Precision*round(solution/(max(max(solution))*PARAMS.Precision));
[row,col,val] = find(solution);
clear solution
graph=[col-1,row-1,val];
weighted='w';

% Optimize louvain NbLouvain times
lnk = zeros(PARAMS.NbNodes, PARAMS.NbLouvain);
lnkS = zeros(PARAMS.NbLouvain,1);
precision = PARAMS.Precision;
if PARAMS.ComputeParallel
    parfor l=1:PARAMS.NbLouvain
        [stability, nb_comm, communities] = stability_louvain_LNL(graph, 1, precision, weighted);
        lnk(:,l) = communities;
        lnkS(l) = stability;
    end
else
    for l=1:PARAMS.NbLouvain
        [stability, nb_comm, communities] = stability_louvain_LNL(graph, 1, precision, weighted);
        lnk(:,l) = communities;
        lnkS(l) = stability;
    end
end
[S,indexS]=max(lnkS);
C=lnk(:,indexS);
N=max(C)+1;

clear communities;
clear graph;

ComputeVI = PARAMS.ComputeVI;
NbLouvain = PARAMS.NbLouvain;
NbNodes = PARAMS.NbNodes;
M = PARAMS.M ;

if ComputeVI && nnz(max(lnk)==NbNodes-1)~=NbLouvain && nnz(max(lnk)==0)~=NbLouvain
     VI = computeRobustness(lnk, lnkS, M, PARAMS.ComputeParallel);
else
    VI = 0;
end

clear lnk;

end

%------------------------------------------------------------------------------
function [S, N, C, VI, VAROUT] = louvain_LCL(Graph, time, PARAMS)
VAROUT =[]; % init varying outputs 

if PARAMS.directed == true
        % Directed part so far not implemented, different variants are possible
        error(['This parameter combination is not allowed! If you want to use directed stability, use the normalised Laplacian']);
end


ComputeVI = PARAMS.ComputeVI  ;
precision = PARAMS.Precision;
NbLouvain = PARAMS.NbLouvain;
M = PARAMS.M ;
NbNodes = PARAMS.NbNodes;
ComputeParallel = PARAMS.ComputeParallel;
weighted = 'w';

% Optimize louvain NbLouvain times
lnk = zeros(NbNodes, NbLouvain);
lnkS = zeros(NbLouvain,1);
if ComputeParallel
    parfor l=1:NbLouvain
        [stability, nb_comm, communities] = stability_louvain_LCL(Graph, time, precision, weighted);
        lnk(:,l) = communities;
        lnkS(l) = stability;
    end
else
    for l=1:NbLouvain
        [stability, nb_comm, communities] = stability_louvain_LCL(Graph, time, precision, weighted);
        lnk(:,l) = communities;
        lnkS(l) = stability;
    end
end
[S,indexS]=max(lnkS);
C=lnk(:,indexS);
N=max(C)+1;

clear communities;
clear Graph;

if ComputeVI && nnz(max(lnk)==NbNodes-1)~=NbLouvain && nnz(max(lnk)==0)~=NbLouvain
     VI = computeRobustness(lnk, lnkS, M,ComputeParallel);
else
    VI = 0;
end

clear lnk;

end
%------------------------------------------------------------------------------
function [S, N, C, VI, VAROUT] = louvain_LNL(Graph, time, PARAMS)
VAROUT =[]; % init varying outputs 

weighted = 'w';


% directed case: M_ij >> from i to j
if PARAMS.directed == true
    Graph = sparse(Graph(:,1)+1,Graph(:,2)+1,Graph(:,3),PARAMS.NbNodes,PARAMS.NbNodes);
    % "transition matrix" and pi unknown
    if PARAMS.precomputed == false
        dout = sum(Graph,2);
        dangling = (dout==0);
        dout(dangling) = 1;
        Dout = sparse(diag(dout));
        clear dout;
        M = (1-PARAMS.teleport_tau)*(Dout\Graph); % deterministic part of transition
        % teleportation according to arXiv:0812.1770
        M =	M + diag(PARAMS.teleport_tau + dangling.*(1-PARAMS.teleport_tau))...
            * ones(PARAMS.NbNodes)/PARAMS.NbNodes;
        
        clear Dout dangling
        [v lambda_all] = eigs(M'); % largest eigenvalue of transition matrix corresponds to stat.distribution.
        lambda = max(diag(lambda_all));
        v = v(:,diag(lambda_all) == lambda);
        v = abs(v);              % make sure eigenvector is positive
        clear lambda;
        % store results for future use
        VAROUT.precomputed = true;
        VAROUT.pi = v/sum(v);
        VAROUT.P = M;        
        %symmetrise due to trace optimsation
        solution = diag(VAROUT.pi)*M;
        solution = (solution +solution')/2;
    else
        solution = diag(PARAMS.pi)* PARAMS.P;
        solution = (solution +solution')/2;
    end
    [row,col,val] = find(solution);
    clear solution
    Graph = [col-1,row-1,val];   
end



% Optimize louvain NbLouvain times
lnk = zeros(PARAMS.NbNodes, PARAMS.NbLouvain);
lnkS = zeros(PARAMS.NbLouvain,1);
if PARAMS.ComputeParallel
    parfor l=1:PARAMS.NbLouvain
        [stability, ~, communities] = stability_louvain_LNL(Graph, time, PARAMS.Precision, weighted);
        lnk(:,l) = communities;
        lnkS(l) = stability;
    end
else
    for l=1:PARAMS.NbLouvain
        [stability, ~, communities] = stability_louvain_LNL(Graph, time, PARAMS.Precision, weighted);
        lnk(:,l) = communities;
        lnkS(l) = stability;
    end
end
[S,indexS]=max(lnkS);
C=lnk(:,indexS);
N=max(C)+1;

clear communities;
clear Graph;

if PARAMS.ComputeVI && nnz(max(lnk)==PARAMS.NbNodes-1)~=PARAMS.NbLouvain...
        && nnz(max(lnk)==0)~=PARAMS.NbLouvain
    VI = computeRobustness(lnk,lnkS, PARAMS.M,PARAMS.ComputeParallel);
else
    VI = 0;
end

clear lnk;

end
%------------------------------------------------------------------------------
function Graph = check(Graph, verbose, PARAMS)
% Check that the graph is properly encoded.
    if verbose
        disp(' ');
        disp('   Graph sanity check...');
    end
    
    % Initialisation of Graph properties
    edgelist = false;
    unweighted = false;
    
    if size(Graph,2) < 2 || size(Graph,1) < 3
        error(['The size of the graph is [' num2str(size(Graph,1)) ',' num2str(size(Graph,2)) ']. Please check that it has been correclty encoded.'])
    end       
    
    if size(Graph,2) ~= size(Graph,1)
        if size(Graph,2) ~=2 && size(Graph,2) ~=3
            error('Wrong size for G: G should be a graph saved either as a list of edges (size(G)=[N,3] if weighted, size(G)=[N,2] if unweighted) or as an adjacency matrix (size(G)=[N,N])');
        end
        edgelist = true;
        if size(Graph,2) == 2
            unweighted = true;
        end
    end
    
    % Check nodes numbering and convert edgelist into adjacency matrix
    if edgelist
        if min(min(Graph(:,1:2))) ~=0
            warning('The numbering of the nodes in the graph should always start with zero.');
            old_node_1 = min(min(Graph(:,1:2)));
            Graph(:,1)=Graph(:,1)-old_node_1;
            Graph(:,2)=Graph(:,2)-old_node_1;
        end
        if unweighted == false
            Graph=sparse(Graph(:,1)+1,Graph(:,2)+1,Graph(:,3));
        else
            Graph=sparse(Graph(:,1)+1,Graph(:,2)+1,ones(length(Graph(:,1)),1));
        end
    end

    % Check that graph contains just numbers
    if any(any(~isnumeric(Graph)))
	error('The graph provided contains elements which are not numbers (isnumeric == false). Please check your graph, and try again.');
    end
        
    % Check symmetry of the adjacency matrix if graph is not directed
    if PARAMS.directed == false
    	if size(Graph,1) ~= size(Graph,2)
        	error('The graph provided is a directed graph. Specify the correct options or correct your graph');
    	end
	if any(any(Graph~=Graph'))
		if nnz(triu(Graph,1))>0 && nnz(tril(Graph,-1))>0
		    error('The graph provided is a directed graph.');
		else
		    warning('Adjacency matrix A of provided graph is triangular -- symmetrizing A = A + A^T');
		    Graph=Graph+Graph';
		end
	end
    end
    
    % Check for isolated nodes
    if ( any( sum(abs(Graph))' == 0 & sum(abs(Graph),2) == 0 ) )
        warning('There are isolated nodes in the graph!?');
    end
    
    % Check for disconnected components
    if exist('graphconncomp','file') == 2
        nbcomp=graphconncomp(sparse(Graph),'WEAK',true);
        if nbcomp>1
            warning(['There are ' num2str(nbcomp) ' not strongly connected components in the graph. If your graph is directed please be aware of the teleportation settings.']);
        end
    end
    
    % Return Graph to its original form
    if edgelist
        [row, col, val] = find(Graph);
        if unweighted
            Graph=[col-1, row-1];
        else
            Graph=[col-1, row-1, val];
        end
    else
        Graph=sparse(Graph);
    end
end
%------------------------------------------------------------------------------
function VI = computeRobustness(lnk, lnkS, M,ComputeParallel)

% Parameter

[~,i] = sort(lnkS);
lnk=lnk(:,i);
lnk=lnk(:,end-M+1:end);
VI = varinfo(lnk',ComputeParallel);
clear i;
end

%------------------------------------------------------------------------------
function [] = stability_plot(Time,t,S,N,VI,ComputeVI,figure_handle)

set(0,'CurrentFigure',figure_handle);

if ComputeVI
    subplot(2,1,1), ax=plotyy(Time(1:t),N(1:t),Time(N>1),S(N>1));
else
    ax=plotyy(Time(1:t),N(1:t),Time(N>1),S(N>1));
end
xlabel('Markov time');
set(ax(1),'YScale','log');
set(ax(2),'YScale','log');
set(ax(1),'YTickMode','auto','YTickLabelMode','auto','YMinorGrid','on');
set(ax(2),'YTickMode','auto','YTickLabelMode','auto','YMinorGrid','on');
set(get(ax(1),'Ylabel'),'String','Number of communities');
set(get(ax(2),'Ylabel'),'String','Stability');
set(ax(1),'XLim', [10^floor(log10(Time(1))) 10^ceil(log10(Time(end)))], 'YLim', [1 10^ceil(log10(max(N)))], 'XScale','log','XMinorGrid','on');
set(ax(2),'XLim', [10^floor(log10(Time(1))) 10^ceil(log10(Time(end)))], 'YLim', [10^floor(log10(min(S(N>1)))), 1], 'XScale','log');
ylabel('Number of communities');
if ComputeVI 
    subplot(2,1,2), semilogx(Time(1:t),VI(1:t));
    set(gca, 'XLim', [10^floor(log10(Time(1))) 10^ceil(log10(Time(end)))], 'YMinorGrid','on','XMinorGrid','on');
    if max(VI)>0
        set(gca,'YLim', [0 max(VI)*1.1]);
    end
    xlabel('Markov time');
    ylabel('Variation of information');
end
drawnow;

end
