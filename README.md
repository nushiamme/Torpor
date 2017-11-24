**Paper authors: Anusha Shankar\*, Rebecca J Schroeder\*, Susan M Wethington, Catherine H Graham, Donald R Powers**

\*Equal authors

Code by: Anusha Shankar, github/nushiamme; contact: anusha<dot>shankar<at>stonybrook<dot>edu for questions about code

#### Code organisation

The code is organized by data type into 3 scripts, for
1. Plots using respirometry data directly,
2. Temperature plots,
3. Phylogenetic corrections and MCMCglmm models and plots.
Figures not listed here were conceptual figures made in powerpoint.

-   **Torpor\_respirometry\_plots.R** - Needs input files *"Torpor\_individual\_summaries.csv"* and *"Broadbill.csv"*; contains code for
    -   *Figure 3*: VO2 for broad-billed hummingbirds showing inflection point at 14-15 degC
    -   *Supplementary Figure S5*: Duration of torpor per individual per night, as a function of the hour of entry
    -   *Supplementary Figure S6*: Average hourly torpid energy savings relative to normothermy
    -   *Supplementary Figure S7*: Average mass-corrected energy expenditure as a function of minimum hourly chamber temperature
    -   *Supplementary Figure S9*: Total nighttime energy expenditure including rewarming (small black dots), and excluding rewarming (coloured dots), as a function of individual mass
-   **Torpor\_temperature\_plots.R** - Needs input files *"Tc\_AllSites\_summ.csv"* and *"Ta\_AllSites\_summ.csv"*; contains code for
    -   *Figure 2*: Chamber Temp plots by hour, per site
    -   *Supplementary Figure S3*: Ambient temp plots by hour, per site
-   **Torpor\_phylo\_MCMCglmm\_models.R** - Needs input files *"Torpor\_individual\_summaries.csv"* and *"hum294.tre"* hummingbird phylogeny file from McGuire et al. 2014 (contact Anusha Shankar for access to this file). Contains code for all the MCMCglmm models, analyses for the following tables, and plots
    -   *Supplementary Table S2*: Stepwise model DIC values for parameter combinations in the nighttime energy expenditure MCMCglmm models
    -   *Supplementary Table S3*: Comparing MCMCglmm stepwise model results for the rewarming models
    -   *Supplementary Figure S4*: MCMCglmm model for the probability of entry into torpor as a function of individual capture mass
    -   *Supplementary Figure S8*: MCMCglmm model for the nighttime energy expenditure as a function of torpor duration and minimum chamber temperature
    -   *Supplementary Figure S10*: MCMCglmm model for rewarming energy expenditure (kJ) as a function of mass and chamber temperature

#### Packages you will need for respirometry and temperature scripts:

    + ggplot2
    + reshape

#### Packages you will need for phylo\_MCMCglmm script

    + MCMCglmm
    + nlme
    + ape
    + geiger # for treedata() function
    + caper
    + coda
