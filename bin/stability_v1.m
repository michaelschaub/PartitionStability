function [S, N, VI, C] = stability_v1(G, T, varargin)
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
%        laplacian      Allows to choose which type of      'normalized'
%                       laplacian should be used to 
%                       calculate the stability. It can
%                       either be 'combinatorial', or 
%                       'normalised'.
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



%$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%
%$                                          $%
%$          Default parameters              $%
%$                                          $%
%$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%

Graph = [];                                     % List of edges of the graph to be partitioned
Time = 1;                                       % Markov times at which the graph should be partitioned
StabilityFunction = @louvain_LNL;    	        % Linearised stability with normalised laplacian is used by default
ComputeVI = true;                               % True if the variation of information should be computed
OutputFile = false;                             % No output file by default.
NbLouvain = 100;                                % Number of louvain optimisations at each Markov time
NbNodes = 0;                                    % Total number of nodes;
Full = false;                                   % If true, performs the full stability
Sanity = true;                                  % If true, performs the graph sanity checks
Precision = 1e-10;                              % Threshold for stability and edges weigths
plotStability = false;                          % If true, plots the results of the stability, number of communities and variation of information vs Markov time.           
verbose = false;                                % Toggles verbose mode
weighted = 'u';                                 % 'u' for unweighted graph, 'w' for weighted graph
prefix = '';                                    % Output prefix
M = 100;                                        % Top M partitions among the L found by louvain are used to compute the variation of information
ComputeParallel = false;                        % Toggles the computation in parallel
TextOutput = false;                             % Toggles the computation in parallel
flag_matlabpool = false;


%$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%
%$                                          $%
%$          Arguments parsing               $%
%$                                          $%
%$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%

% Options

if nargin > 2
    [StabilityFunction, ComputeVI, OutputFile, prefix, NbLouvain, M, Full, Sanity, Precision, plotStability, verbose, ComputeParallel, TextOutput] = parseinput(length(varargin),G,varargin);
end

% Display starting time
if verbose
    c = clock;
    disp(' ');
    disp(['   Stability will be computed with the following parameters:']);
    disp(' ');
    if Full; disp('      - Full Stability'); else disp('      - Linearised Stability'); end
    if strfind(func2str(StabilityFunction), 'N'); disp(['      - Normalised laplacian']); else disp(['      - Combinatorial laplacian']); end
    if ComputeVI; disp(['      - Computation of the variation of information (VI): Yes']); else disp(['      - Computation of the variation of information: No']); end
    if OutputFile; disp(['      - Save the results in a file: Yes ' prefix]); else disp(['      - Save the results in a file: No']); end
    if OutputFile; disp(['      - Saved file prefix: ' prefix]); end
    if Sanity; disp(['      - Check the input graph: Yes']); else disp(['      - Check the input graph: No']); end
    if plotStability; disp(['      - Plot the results: Yes']); else disp(['      - Plot the results: No']); end
    if verbose; disp(['      - Verbose mode: Yes']);    else disp(['      - Verbose mode: No']);    end
    if ComputeParallel; disp(['      - Parallel computation: Yes']); else disp(['      - Parallel computation: No']); end
    disp(['      - Number of Louvain iterations: ' int2str(NbLouvain)]);
    disp(['      - Number of Louvain iterations used for the computation of VI: ' int2str(M)]);
    disp(['      - Precision used: ' num2str(Precision)]); 
    disp(' ');
    tstart=tic;
end


% Argument 1: G

if nargin > 0
    if size(G,1) == size(G,2) && ~issparse(G)
            G=sparse(G);
    end
    % Check if the graph is correctly encoded
    if Sanity
        G=check(G, verbose);
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
        if max(max(Graph)) == 1 && isinteger(Graph)
            weighted='u';
        else
            weighted='w';
        end
        NbNodes = size(Graph,2);
    % if the linearised stability is to be computed, Graph should be the
    % list of edges.
    else
        if size(G,1) == size(G,2)
            [rows,cols,vals] = find(G);
            if sum(vals)==length(vals)
                Graph=[cols-1, rows-1];
                weighted='u';
            else
                Graph=[cols-1, rows-1, vals];
                weighted='w';
            end
            clear rows cols vals;
        elseif size(G,2)==2 
            Graph=G;
            weighted='u';
        elseif size(G,2)==3
            Graph=G;
            weighted='w';
        else
            error('Wrong size for G: G should be a graph saved either as a list of edges (size(G)=[N,3] if weighted, size(G)=[N,2] if unweighted) or as an adjacency matrix (size(G)=[N,N])');
        end
        NbNodes = max(Graph(:,1))+1;
    end
else
    error('Please provide at least the graph to be partitioned. Type "help stability" for more information.');
end    

% Argument 2: T

if nargin > 1
    if (isvector(T) && isnumeric(T))
        Time=T;
    else
        error('The second argument should be a vector. Type "help stability" for more information.');
    end
end

% Parallel computation: Initialize the number of cores if matlabpool is not
% yet running.
if ComputeParallel && matlabpool('size') == 0
    flag_matlabpool = true;
    matlabpool
end

%$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%
%$                                          $%
%$      Computation of the stability        $%
%$                                          $%
%$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%

% Display starting time
if verbose
    c = clock;
    disp(' ');
    disp(['   Partitioning of the graph started at ' datestr([2011 1 1 c(4) c(5) c(6)], 'HH:MM:SS')]);
    disp(' ');
    tstart=tic;
end


% Initialisation
S = zeros(1, length(Time));
N = zeros(1, length(Time));
VI = zeros(1, length(Time));
C = zeros(NbNodes, length(Time));

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
    
    [S(t), N(t), C(:,t), VI(t)] = StabilityFunction(Graph, Time(t), Precision, weighted, ComputeVI, NbLouvain, M, NbNodes, ComputeParallel);
    
    if plotStability && t>1
        stability_plot(Time,t,S,N,VI,ComputeVI,figure_handle);
    end
    
    if TextOutput
        cd(['Partitions_' prefix]);
        dlmwrite(['Partition_' prefix '_' num2str(Time(t),'%10.6f') '.dat'],[[1:NbNodes]',C(:,t)],'delimiter','\t');
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
function [StabilityFunction, ComputeVI, OutputFile, prefix, NbLouvain, M, Full, Sanity, Precision, plotStability, verbose, ComputeParallel, TextOutput] = parseinput(options,G,varargin)
% Parse the options

% Initialise parameters

% if (size(G,1) == size(G,2) && nnz(G)<5000 && size(G,1)<1000) || (size(G,1) ~= size(G,2) && size(G,1)<5000 && max(max(G(:,1:2)))<1000)
%     StabilityFunction = @louvain_FNL;     % Linearised stability with normalised laplacian is used by default
% else
%     StabilityFunction = @louvain_LNL;     % Linearised stability with normalised laplacian is used by default
% end

stability_type_specified=false;
threshold_nnodes_full = 1000;
threshold_nedges_full = 5000;
ComputeVI = true;                               % True if the variation of information should be computed
OutputFile = false;                             % No output file by default.
NbLouvain = 100;                                % Number of louvain optimisations at each Markov time
Full = false;
Laplacian = 'Normalised';
Sanity = true;
Precision = 1e-10; 
plotStability = false; 
verbose = false;
M = 100;
prefix = '';
ComputeParallel = false;
TextOutput = false;
attributes={'novi', 'l', 'm', 'out', 'full', 'linearised', 'nocheck', 'laplacian', 'prec', 'plot','v','t','p'};

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
            elseif strcmpi(varargin{i},'novi')
                ComputeVI = false;
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
            elseif strcmpi(varargin{i},'p')
                if exist('matlabpool','file')
                    ComputeParallel = true;
                else
                    ComputeParallel = false;
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
                        NbLouvain = round(varargin{i+1});
                        M = round(varargin{i+1});
                    end
                elseif strcmpi(varargin{i},'prec')
                    if isnumeric(varargin{i+1})
                        Precision = varargin{i+1};
                    end
                elseif strcmpi(varargin{i},'m')
                    if isnumeric(varargin{i+1})
                        M = varargin{i+1};
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

%------------------------------------------------------------------------------
function [S, N, C, VI] = louvain_FNL(Graph, time, precision, weighted, ComputeVI, NbLouvain, M, NbNodes, ComputeParallel)
% Computes the full normalised stabilty

% Generate the matrix exponential
diagdeg=sparse((diag(sum(Graph)))/sum(sum(Graph)));  %diag matrix with stat distr
trans=sparse(diag(    (sum(Graph)).^(-1)     ) * Graph);  %(stochastic) transition matrix
clear Graph;
Lap=sparse(trans-eye(NbNodes));
clear trans;
exponential=sparse(expm(time.*Lap));
clear Lap;
solution=sparse(diagdeg*exponential);
clear exponential;
clear diagdeg;
solution=max(max(solution))*precision*round(solution/(max(max(solution))*precision));
clear exponential;
clear diagdeg;
[row,col,val] = find(solution);
clear solution
graph=[col-1,row-1,val];

% Optimize louvain NbLouvain times
lnk = zeros(NbNodes, NbLouvain);
lnkS = zeros(NbLouvain,1);
if ComputeParallel
    parfor l=1:NbLouvain
	[stability, nb_comm, communities] = stability_louvain_LNL(graph, 1, precision, weighted);
        lnk(:,l) = communities;
        lnkS(l) = stability;
    end
else
    for l=1:NbLouvain
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

if ComputeVI && nnz(max(lnk)==NbNodes-1)~=NbLouvain && nnz(max(lnk)==0)~=NbLouvain
     VI = computeRobustness(lnk, lnkS, M, ComputeParallel);
else
    VI=0;
end

clear lnk;

end


%------------------------------------------------------------------------------
function [S, N, C, VI] = louvain_FCL(Graph, time, precision, weighted, ComputeVI, NbLouvain, M, NbNodes, ComputeParallel)
% Computes the full combinatorial stability

% Generate the matrix exponential
diagdeg=sparse(1/NbNodes .* eye(NbNodes));
Lap=sparse(Graph-diag(sum(Graph)));
total=sum(sum(Graph));
clear Graph;
avdegree=total/NbNodes;
exponential=sparse(expm((time/avdegree).*Lap));
clear Lap;
solution=sparse(diagdeg*exponential);
solution=max(max(solution))*precision*round(solution/(max(max(solution))*precision));
clear exponential;
clear diagdeg;
[row,col,val] = find(solution);
clear solution
graph=[col-1,row-1,val];


% Optimize louvain NbLouvain times
lnk = zeros(NbNodes, NbLouvain);
lnkS = zeros(NbLouvain,1);
if ComputeParallel
    parfor l=1:NbLouvain
        [stability, nb_comm, communities] = stability_louvain_LNL(graph, 1, precision, weighted);
        lnk(:,l) = communities;
        lnkS(l) = stability;
    end
else
    for l=1:NbLouvain
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

if ComputeVI && nnz(max(lnk)==NbNodes-1)~=NbLouvain && nnz(max(lnk)==0)~=NbLouvain
     VI = computeRobustness(lnk, lnkS, M, ComputeParallel);
else
    VI = 0;
end

clear lnk;

end

%------------------------------------------------------------------------------
function [S, N, C, VI] = louvain_LCL(Graph, time, precision, weighted, ComputeVI, NbLouvain, M, NbNodes, ComputeParallel)


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
     VI = computeRobustness(lnk, lnkS, M, ComputeParallel);
else
    VI = 0;
end

clear lnk;

end
%------------------------------------------------------------------------------
function [S, N, C, VI] = louvain_LNL(Graph, time, precision, weighted, ComputeVI, NbLouvain, M, NbNodes, ComputeParallel)

% Optimize louvain NbLouvain times
lnk = zeros(NbNodes, NbLouvain);
lnkS = zeros(NbLouvain,1);
if ComputeParallel
    parfor l=1:NbLouvain
        [stability, nb_comm, communities] = stability_louvain_LNL(Graph, time, precision, weighted);
        lnk(:,l) = communities;
        lnkS(l) = stability;
    end
else
    for l=1:NbLouvain
        [stability, nb_comm, communities] = stability_louvain_LNL(Graph, time, precision, weighted);
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
    VI = computeRobustness(lnk,lnkS, M, ComputeParallel);
else
    VI = 0;
end

clear lnk;

end
%------------------------------------------------------------------------------
function Graph = check(Graph, verbose)
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

    % Check for NaN's
    if any(any(isnan(Graph)))
	error('The graph provided contains elements which are not numbers (NaN). Please check your graph, and try again.');
    end
        
    % Check symmetry of the adjacency matrix
    if size(Graph,1) ~= size(Graph,2)
        error('The graph provided is a directed graph. This program only deals with undirected graphs.');
    end
    if any(any(Graph~=Graph'))
        if nnz(triu(Graph,1))>0 && nnz(tril(Graph,-1))>0
            error('The graph provided is a directed graph. This program only deals with undirected graphs.');
        else
            warning('Adjacency matrix A of provided graph is triangular -- symmetrizing A = A + A^T');
            Graph=Graph+Graph';
        end
    end
    
    % Check for isolated nodes
    if nnz(sum(Graph))~=size(Graph,2)
        warning('There are isolated nodes in the graph');
    end
    
    % Check for disconnected components
    if exist('graphconncomp','file') == 2
        nbcomp=graphconncomp(sparse(Graph));
        if nbcomp>1
            warning(['There are ' num2str(nbcomp) ' disconnected components in the graph.']);
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
function VI = computeRobustness(lnk, lnkS, M, ComputeParallel)

[i,i] = sort(lnkS);
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
