# Rpath: an open source food web model

Rpath is an implementation of the ecosystem model Ecopath with Ecosim
(EwE; Polovina ([1984](#ref-polovina_model_1984)); Christensen and Pauly
([1992](#ref-christensen_ecopath_1992)); Walters *et al.*
([1997](#ref-walters_structuring_1997)); Christensen and Walters
([2004](#ref-christensen_ecopath_2004)); Walters and Martell
([2004](#ref-walters_fisheries_2004)); Walters *et al.*
([2008](#ref-walters_ecosim_2008))).

Articles for [model
setup](https://noaa-edab.github.io/Rpath/articles/ModelSetup.md),
[running Rpath](https://noaa-edab.github.io/Rpath/articles/RunRpath.md)
and [running
Rsim](https://noaa-edab.github.io/Rpath/articles/RunRsim.md) describe
some of the basic functionality of the package using a fictional
ecosystem, R Ecosystem. Any resemblance to an actual ecosystem is purely
coincidental.

To see the underlying mathematics please refer to Lucey *et al.*
([2020](#ref-lucey_conducting_2020)). For additional features, see
Whitehouse and Aydin ([2020](#ref-whitehouse_assessing_2020)) for
ecosense and Lucey *et al.* ([2021](#ref-lucey_evaluating_2021)) for
management strategy evaluation capability.

## Installation

If necessary, first install the `remotes()` package:
`install.packages("remotes")`

    remotes::install_github("noaa-edab/Rpath",build_vignettes=TRUE)` to install Rpath and the vignettes.

If you experience issues installing the package, please use this
alternative

    pak::pak("noaa-edab/Rpath")

## References

Christensen, V., and Pauly, D. 1992. ECOPATH II — a software for
balancing steady-state ecosystem models and calculating network
characteristics. Ecological Modelling, 61: 169–185.
<https://www.sciencedirect.com/science/article/pii/0304380092900168>
(Accessed 9 April 2025).

Christensen, V., and Walters, C. J. 2004. Ecopath with Ecosim: Methods,
capabilities and limitations. Ecological Modelling, 172: 109–139.
<https://www.sciencedirect.com/science/article/pii/S030438000300365X>
(Accessed 9 April 2025).

Lucey, S. M., Gaichas, S. K., and Aydin, K. Y. 2020. Conducting
reproducible ecosystem modeling using the open source mass balance model
Rpath. Ecological Modelling, 427: 109057.
<https://www.sciencedirect.com/science/article/pii/S0304380020301290>
(Accessed 9 April 2025).

Lucey, S. M., Aydin, K. Y., Gaichas, S. K., Cadrin, S. X., Fay, G.,
Fogarty, M. J., and Punt, A. 2021. Evaluating fishery management
strategies using an ecosystem model as an operating model. Fisheries
Research, 234: 105780.
<https://www.sciencedirect.com/science/article/pii/S0165783620302976>
(Accessed 9 April 2025).

Polovina, J. J. 1984. Model of a coral reef ecosystem. Coral Reefs, 3:
1–11. <https://doi.org/10.1007/BF00306135> (Accessed 9 April 2025).

Walters, C., Christensen, V., and Pauly, D. 1997. Structuring dynamic
models of exploited ecosystems from trophic mass-balance assessments.
Reviews in Fish Biology and Fisheries, 7: 139–172.
<https://doi.org/10.1023/A:1018479526149> (Accessed 9 April 2025).

Walters, C., Martell, S. J. D., Christensen, V., and Mahmoudi, B. 2008.
An Ecosim Model for Exploring Gulf of Mexico Ecosystem Management
Options: Implications of Including Multistanza Life-History Models for
Policy Predictions. Bulletin of Marine Science, 83: 251–271.

Walters, C. J., and Martell, S. J. D. 2004. Fisheries Ecology and
Management Princeton University Press.
<https://press.princeton.edu/books/paperback/9780691115450/fisheries-ecology-and-management>
(Accessed 9 April 2025).

Whitehouse, G. A., and Aydin, K. Y. 2020. Assessing the sensitivity of
three Alaska marine food webs to perturbations: An example of Ecosim
simulations using Rpath. Ecological Modelling, 429: 109074.
<https://www.sciencedirect.com/science/article/pii/S0304380020301460>
(Accessed 9 April 2025).
