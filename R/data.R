

#' AB.params
#'
#' Anchovy Bay (\code{AB.params}) is a sample ecosystem frequently used as an
#' example ecosystem with EwE (see Christensen and Walters, 2024).
#'
#' @format An unbalanced Rpath model object that contains:
#' \describe{
#'   \item{model}{A data.table containing Ecopath unbalanced ecosystem
#'   parameters (base parameters and fisheries).}
#'   \item{diet}{A data.table containing the Ecopath model's diet matrix.}
#'   \item{stanzas}{Parameters for multistanza groups.}
#'   \item{pedigree}{A data.table containing the data quality (pedigree) for the
#'   Ecopath model.}
#'   ...
#' }
#' @source Christensen V, Walters CJ (2024) Ecosystem Modelling with EwE. The
#'   University of British Columbia, Vancouver, BC. doi:10.14288/24d7-ab68
#' @examples
#' # Balance the model
#' AB_bal <- rpath(AB.params)
#'  
"AB.params"

#' REco.params
#'
#' Rpath inputs for a tutorial Rpath ecosystem.
#'
#' @format An unbalanced Rpath model object that contains:
#' \describe{
#'   \item{model}{A data.table containing Ecopath unbalanced ecosystem
#'   parameters (base parameters and fisheries).}
#'   \item{diet}{A data.table containing the Ecopath model's diet matrix.}
#'   \item{stanzas}{Parameters for multistanza groups.}
#'   \item{pedigree}{A data.table containing the data quality (pedigree) for the
#'   Ecopath model.}
#'   ...
#' }
#' @seealso \code{vignette("ModelSetup", package = "Rpath")} for an example
#' of how to build this tutorial ecosystem.
#' @references Lucey SM, Gaichas SK, Aydin KY (2020) Conducting reproducible
#'   ecosystem modeling using the open source mass balance model Rpath. Ecol
#'   Model 427:11. doi:10.1016/j.ecolmodel.2020.109057
#' @examples
#' # Balance the model
#' REco_bal <- rpath(REco.params)
#' 
"REco.params"

#' Eastern Bering Sea 1990s Ecopath model 
#'
#' Rpath inputs for the eastern Bering Sea 1990s Ecopath model (53 biological 
#' groups and 1 fleet).
#'
#' @format An unbalanced Rpath model object that contains:
#' \describe{
#'   \item{model}{A data.table containing Ecopath unbalanced ecosystem
#'   parameters (base parameters and fisheries).}
#'   \item{diet}{A data.table containing the Ecopath model's diet matrix.}
#'   \item{stanzas}{Parameters for multistanza groups.}
#'   \item{pedigree}{A data.table containing the data quality (pedigree) for the
#'   Ecopath model.}
#'   ...
#' }
#' @note This is an aggregated version of the 1990s eastern Bering Sea Ecopath
#'   model of Aydin et al. (2007) and does not include stanzas.
#' @source Whitehouse and Aydin 2020.  Assessing the sensitivity of three Alaska
#'   marine food webs to perturbations: an example of Ecosim simulations using
#'   Rpath. https://doi.org/10.1016/j.ecolmodel.2020.109074
#' @references Aydin KY, Gaichas S, Ortiz I, Kinzey D, Friday N (2007) A
#'   comparison of the Bering Sea, Gulf of Alaska, and Aleutian Islands large
#'   marine ecosystems through food web modeling. U.S. Dep. Commer., NOAA Tech.
#'   Memo. NMFS-AFSC-178.
#' @examples
#' # Balance the model
#' ebs_bal <- rpath(Ecosense.EBS)
#' 
"Ecosense.EBS"

#' Gulf of Alaska (west/central) 1990s Ecopath model 
#'
#' Rpath inputs for the Gulf of Alaska (west/central) 1990s Ecopath model (49 
#' biological groups and 1 fleet).
#'
#' @format An unbalanced Rpath model object that contains:
#' \describe{
#'   \item{model}{A data.table containing Ecopath unbalanced ecosystem
#'   parameters (base parameters and fisheries).}
#'   \item{diet}{A data.table containing the Ecopath model's diet matrix.}
#'   \item{stanzas}{Parameters for multistanza groups.}
#'   \item{pedigree}{A data.table containing the data quality (pedigree) for the
#'   Ecopath model.}
#'   ...
#' }
#' @note This is an aggregated version of the 1990s Gulf of Alaska Ecopath model
#'   of Aydin et al. (2007) and does not include stanzas.
#' @source Whitehouse and Aydin 2020.  Assessing the sensitivity of three Alaska
#'   marine food webs to perturbations: an example of Ecosim simulations using
#'   Rpath. https://doi.org/10.1016/j.ecolmodel.2020.109074
#' @references Aydin KY, Gaichas S, Ortiz I, Kinzey D, Friday N (2007) A
#'   comparison of the Bering Sea, Gulf of Alaska, and Aleutian Islands large
#'   marine ecosystems through food web modeling. U.S. Dep. Commer., NOAA Tech.
#'   Memo. NMFS-AFSC-178.
#' @examples
#' # Balance the model
#' goa_bal <- rpath(Ecosense.GOA)
#' 
"Ecosense.GOA"

#' Eastern Chukchi Sea Ecopath model 
#'
#' Rpath inputs for the eastern Chukchi Sea Ecopath model (52 biological groups 
#' and 1 fleet).
#'
#' @format An unbalanced Rpath model object that contains:
#' \describe{
#'   \item{model}{A data.table containing Ecopath unbalanced ecosystem
#'   parameters (base parameters and fisheries).}
#'   \item{diet}{A data.table containing the Ecopath model's diet matrix.}
#'   \item{stanzas}{Parameters for multistanza groups.}
#'   \item{pedigree}{A data.table containing the data quality (pedigree) for the
#'   Ecopath model.}
#'   ...   
#' }
#' @note This is an aggregated version of the eastern Chukchi Sea Ecopath model
#'   of Whitehouse and Aydin (2016).
#' @source Whitehouse and Aydin 2020.  Assessing the sensitivity of three Alaska
#'   marine food webs to perturbations: an example of Ecosim simulations using
#'   Rpath. https://doi.org/10.1016/j.ecolmodel.2020.109074
#' @references Whitehouse GA, Aydin KY (2016) Trophic structure of the eastern
#'   Chukchi Sea: An updated mass balance food web model. U.S. Dep Commer, NOAA
#'   Tech. Memo. NMFS-AFSC-318. doi:10.7289/V5/TM-AFSC-318
#' @examples
#' # Balance the model
#' ecs_bal <- rpath(Ecosense.ECS)
#' 
"Ecosense.ECS"