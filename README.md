## ERGMs for Brain Networks
Source code for the app **ERGMs for Brain Networks** in [Brainlife.io](https://brainlife.io/).
This app estimates the parameters of an Exponential Random Graph Model (ERGM) on brain network data. 
#### Authors
- Vito Dichio (vito.dichio@etu.sorbonne-universite.fr)
## Basics
#### ERGM in a nutshell
Given a graph G and a set of sufficient statistics $\bm{x}(G)\in\mathbb{R}^m$ and a set of parameters $\bm{\theta}\in\mathbb{R}^m$, an ERGM assigns is specified by the following p.d.f.: 


<img width="305" alt="Screenshot 2021-11-20 at 12 50 53" src="https://user-images.githubusercontent.com/79842912/142725299-befccceb-51af-42b5-bdaa-034718b4fba6.png">


the goal of an ERGM fit is to estimate the parameters $\bm{\theta}$, given an observed graph $G$. Current implementations are based on MCMC-MLE procedures, from which estimates can be drawn together with st.deviations. In addition, model assessment can be performed by simulating synthetic networks from the estimated model and resorting graphical Goodness of Fit (GoF) methods.

#### Implementation
The current app is based on the *ergm* package for *R* language [1].  We refer to the documentation for details on the estimation procedure; useful [tutorials](https://github.com/statnet/Workshops/wiki) are provided by the *statnet* community e.g. [here](http://statnet.org/Workshops/ergm_tutorial.html).
For examples of applications to neuroscience, see for instance [2-3]. The default model implemented here is the one selected in [3].
## Input/Output
### Input
- **Network**: The input brain network must be stored in the [v2 JSON graph format](https://github.com/jsongraph/json-graph-specification). The *edge* objects contain information about the elements of connectivity (or adjacency) matrix, together with any other user-specified edge covariate. In the *node* objects it is possible to store nodal attributes. (*NB*: An example of *R* code to transform raw data into the v2JSON format is provided in *dichio/bl-ERGM/RawTov2JSON/v2JSON-writer.R*)
- **ERGM formula**: String containing the ERGM formula according to the guidelines of the *ergm* package for *R* [1].
- **nsim_gof**: Integer number of simulated networks used for Goodness of Fit (GoF) model assessment.
- **unfiltered**: Boolean variable specifying if the input network is not already discretized into binary states. In the case of connectivity matrices (**unfiltered**=TRUE), the app implements the ECO filtering criterion ($k=3$) [4].
### Output
- **estimation.txt**: Result of the fitting procedure i.e. estimated parameters, st.errors, covariance matrix.
- **log-computation.txt**
- **gof.pdf**: Standard plots for the Goodness of Fit as implemented in the *ergm* package [1].
- **mcmc-diagnostic.pdf**: sanity check for the MCMC procedure. 



## References
[1] Hunter, D. R. *et al.* (2008). ergm: A package to fit, simulate and diagnose exponential-family models for networks. _Journal of statistical software_, _24_(3), nihpa54860.

[2] Simpson, S. L., Hayasaka, S., & Laurienti, P. J. (2011). Exponential random graph modeling for complex brain networks. _PloS one_, _6_(5), e20039.

[3] Obando, C., & De Vico Fallani, F. (2017). A statistical model for brain networks inferred from large-scale electrophysiological signals. _Journal of The Royal Society Interface_, _14_(128), 20160940.

[4] De Vico Fallani, F., Latora, V., & Chavez, M. (2017). A topological criterion for filtering information in complex brain networks. _PLoS computational biology_, _13_(1), e1005305.
