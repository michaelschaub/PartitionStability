###############################################################################
Copyright (C) 2012 A. Delmotte, M. Schaub, S. Yaliraki, M. Barahona

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

###############################################################################

-----------------------------------------------------------------------------
Community Detection using the stability of a graph partition.
-----------------------------------------------------------------------------

The code implements the stability method as introduced in the article

(1) "Stability of graph communities across time scales" 
Delvenne, J.-C.; Yaliraki, S. N. & Barahona, M. 
arXiv:0812.1811 (2009)
and then published in 
Proceedings of the National Academy of Sciences, 2010, 107, 12755-12760;

and further expanded in:

(2) "Laplacian Dynamics and Multiscale Modular Structure in Networks"
Lambiotte, R.; Delvenne, J.-C. & Barahona, M.
arxiv:0812.1770 (2009)

and

(3) J.-C. Delvenne, M. T. Schaub, S. N. Yaliraki, and M. Barahona, 
"The stability of a graph partition: A dynamics-based framework for community 
detection" in Time Varying Dynamical Networks (N. Ganguly, A. Mukherjee, 
M. Choudhury, F. Peruani, and B. Mitra, eds.), Birkhauser, Springer,
2012. to be published.


To optimize the stability quality function, we use the Louvain algorithm as 
described in the publication:

(4) "Fast unfolding of communities in large networks",
Vincent D Blondel, Jean-Loup Guillaume, Renaud Lambiotte, Etienne Lefebvre,
Journal of Statistical Mechanics: Theory and Experiment 2008 (10), P10008

The folder  /demo/ contains three examples to demonstrate the main functionality 
of the program:

(i) demo.mat : A simple network with a multiscale community structure.  
(ii) Protein_Adk.mat : A graph representation of the protein Adenylate Kinase 
  (AdK) as described in

(5) "Protein multi-scale organization through graph partitioning and robustness 
analysis: application to the myosin–myosin light chain interaction"
Delmotte, A.; Tate, E. W.; Yaliraki, S. N. & Barahona, M. 
Physical Biology, 2011, 8, 055010

(ii) ring_of-rings.mat  : The "ring-of-rings" graph introduced in:

(6) "Markov dynamics as a zooming lens for multiscale community detection: 
non clique-like communities and the field-of-view limit"
Schaub, M. T.; Delvenne, J.-C.; Yaliraki, S. N. & Barahona, M. 
PLoS ONE, Public Library of Science, 2012, 7, e32210

Further example graphs are available on request. 

***If you make use of any part of this toolbox, please cite the 
respective articles.***

For detailed instructions on how to compile the code in MATLAB see below.
If you find a bug or have further comments, please send an email and if 
necessary the input file and the parameters that caused the error.

- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Authors   : A. Delmotte and M. Schaub  
Email     : antoine.delmotte09@imperial.ac.uk, michael.schaub09@imperial.ac.uk  
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

###############################################################################

-----------------------------------------------------------------------------
Contributions to the code
-----------------------------------------------------------------------------

The C++ code performing the stability optimization is based on the 
implementation of the Louvain method as available from 
http://sites.google.com/site/findcommunities/ 
(Authors: Jean-Loup Guillaume / Etienne Lefebvre)

The Louvain code was adapted and extended by R. Lambiotte 
(http://www.lambiotte.be) to allow for the optimization of the stability quality 
function and was subsequently refined further by Antoine Delmotte. 

The MATLAB frontend was initally added by Antoine Delmotte. 
Final adjustments and additions, testing, and maintenance are due to 
Antoine Delmotte and Michael Schaub.
Version 2.0 which updates the first release and includes the addition of 
directed stability was created by Michael Schaub (see code history for further
details). 

###############################################################################

-----------------------------------------------------------------------------
How to install the stability package
-----------------------------------------------------------------------------

1. Open Matlab

2. Make sure you have a C++ compiler installed
  * For Linux, you can find one here: 
    http://www.gnu.org/software/gcc/
  * For Windows, you can use Visual C++ express: 
    http://www.microsoft.com/express/Windows/

3. Make sure mex is properly configured in Matlab:
  * Type "mex -setup" in Matlab, and choose your compiler.

4. In Matlab, go into the directory of the Stability toolbox.

5. Type "Install_Stability" in the Matlab command window.
  * If you get an error message concerning the libstdc++.so file, 
    you may want to try the following manipulation:

        cd "Matlab_root_directory"/sys/os/glnx86/
        sudo mv libgcc_s.so.1 libgcc_s.so.1.back
        sudo mv libstdc++.so.6 libstdc++.so.6.back

6. You will get a messge asking whether the stability toolbox should 
   be added to your Matlab path. Answering yes will allow you to use 
   the stability toolbox functions as standard Matlab functions.
            
7. Type "help stability" in Matlab to discover how to use the code.

8. Try this example to check that everything is working:
    
        cd('demo');   % go into the demo directory (in the stability folder)
        load demo;    % load data and then run stability
        [S, N, VI, C] = stability(Graph,Time,'plot','v');

NOTES:

* The install script provides the option to add the bin folder to your 
Matlab path. This will enable you to use stability as a standard Matlab 
function from any directory. If you don't want this option any more,
just remove it from the path by going in File/Set Path.

* If you get a warning message concerning savepath, and you want the 
stability code to be in your path, go, after the installation, in 
File/Set Path, and choose "save". Then choose where you want pathdef.m
to be saved. If at the next matlab startup, you notice that stability is
not in your matlab path anymore, try editing/creating the "startup.m" file
from your matlab user folder (type userpath to know where it is located)
and add the following line: addpath(' path to bin folder of stability 
package '). Alternatively, if you are the only user on your machine, you
can start matlab as a superuser ("sudo matlab" in linux) and rerun the
"Install_Stability" script. This will permanently add the stability folder 
in the path for all users.

* To speed up the calculations, you might consider adding the
option 'noVI'. This disables the calculation of the variation of information, 
which is usually slow at small Markov times, when the number of 
communities found is big. 
Another option is to decrease the number of optimisations on which the variation 
of information is calculated. To do so, add the option 'M' and put a value
such that M < L (L is the number of louvain optimisations).  
Example:
 
    [S, N, VI, C] = stability(Graph,time,'plot','v', 'L', 100, 'M', 10);

