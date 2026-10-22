
# ################################################################################
# #' Plot Biomass Fit for a Single Species and Data Source
# #'
# #' Generates a diagnostic plot comparing simulated biomass against observed
# #' survey data for a specific species and a specific data source. It includes
# #' log-normal confidence intervals for the observed data and displays the NLL score.
# #' 
# #' @details
# #' This function extracts the goodness-of-fit metrics from `rsim.fit.obj()` for a
# #' given species and observation source. It calculates the mean observed biomass
# #' adjusted by survey catchability (`survey_q`) and computes the upper and lower
# #' bounds of the 95% confidence interval assuming a log-normal error distribution.
# #' The simulated biomass is plotted as a continuous line, with observed data as
# #' points and error bars.
# #'
# #' The log-scale standard deviation is derived from the arithmetic-scale
# #' coefficient of variation via `sqrt(log(1 + cv^2))`, and the interval is
# #' drawn as `mn * exp(+/-1.96 * sdlog)`. This treats the observed value `mn`
# #' as the median of the log-normal distribution rather than its mean (for
# #' which `E[X] = exp(mu + sigma^2/2)`), a minor and standard approximation
# #' common in fisheries stock-assessment diagnostics.
# #' 
# #' @inheritParams rsim.plot.fit
# #' @param species A single character string representing the species name.
# #' @param datasource A character string identifying the specific survey or data source
# #'   to extract from the objective function data (e.g., `"race_wgoa"`).
# #'
# #' @return Invoked for its side effect of drawing a plot to the active graphics
# #'   device. Returns `NULL`.
# #'
# #'@export
#rsim.plot.fitbio <- function(Rsim.scenario, Rsim.output, species, datasource, small= FALSE) {
#  
#  if (length(species) > 1) {
#    warning("rsim.plot.fitbio expects a single species string. Using only the first element.")
#    species <- species[1]
#  }
#  
#  if(!is.null(Rsim.scenario$fitting$Biomass)){
#    bio.obj <- rsim.fit.obj(Rsim.scenario, Rsim.output)$Biomass
#  } else {
#    bio.obj <- data.frame(Year = character(), Group = character(), Type = character(), 
#                          Stdev = double(), Value = double(), Scale = double(),
#                          Source = character(), ID = character(), obs = double(), 
#                          sd = double(), initial_q = double(), wt = double())
#  }
#  
#  qdat <- bio.obj[bio.obj$Group == species & bio.obj$Source %in% datasource, ]
#  
#  if (nrow(qdat) == 0 && !is.na(datasource)) {
#    warning(sprintf("No biomass observation data found for species '%s' and datasource '%s'. Check your scenario fitting data.", species, datasource))
#  }
#  
#  if(small){
#    pin <- graphics::par("pin") 
#    dyn_cex <- max(0.75, min(0.95, min(pin) *0.2+ 0.6)) #FLAG ####
#    
#    mar_vals <- c(2.8, 1.5, 1, 0.5)
#    axis1_mgp <- c(3, 0.3, 0) #FLAG####
#    axis2_mgp <- c(3, 0.2, 0)
#    
#  } else {
#    dyn_cex <- 0.9
#    mar_vals <- c(4, 3, 2, 1)
#    axis1_mgp <- c(3, 0.6, 0)
#    axis2_mgp <- c(3, 0.6, 0)
#  }
#  
#  oldpar <- graphics::par(mar = mar_vals)
#  on.exit(graphics::par(oldpar))
#  
#  sq <- qdat$survey_q
#  sq[sq == 0 | is.na(sq)] <- 1e-6
#  
#  mn   <- qdat$obs / sq
#  sdlog <- sqrt(log(1.0 + (qdat$sd / sq) * (qdat$sd / sq) / ifelse(mn == 0, 1e-6, mn * mn))) 
#  up <- mn * exp(1.96 * sdlog)
#  dn <- mn / exp(1.96 * sdlog)
#  
#  est <- Rsim.output$annual_Biomass[, species] 
#  tot  <- sum(qdat$fit * qdat$wt)
#  
#  all_vals <- c(up, est)
#  valid_vals <- all_vals[!is.na(all_vals) & !is.infinite(all_vals)]
#  ymax <- if(length(valid_vals) > 0) max(valid_vals) else 1
#  if (ymax <= 0) ymax <- 1
#  
#  graphics::plot(
#    as.numeric(rownames(Rsim.output$annual_Biomass)),
#    est, type = "l", ylim = c(0, ymax), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n"
#  )
#  
#  graphics::axis(1, mgp = axis1_mgp, tck = -0.04, cex.axis = dyn_cex)
#  graphics::axis(2, mgp = axis2_mgp, tck = -0.04, cex.axis = dyn_cex)
#  
#  if(small){
#    graphics::mtext(species, side = 1, cex = dyn_cex*0.8, line = 1.0, adj = 0) #flag ####
#    graphics::mtext(datasource, side = 1, cex = dyn_cex*0.8, line = 1.8, adj = 0) #Flag####
#    
#    if (nrow(qdat) > 0) {
#      graphics::mtext(
#        sprintf("nll:%.3g  q:%.3g (%s)", tot, unique(sq)[1], unique(qdat$Type)[1]),
#        side = 1, line = 2.6, cex = dyn_cex*0.8, adj = 0
#      )
#    }
#  } else {
#    graphics::mtext(
#      paste(datasource, species, "biomass", sprintf("   nll: %.3g", tot)),
#      side = 1, cex = 0.9, line = 1.5, adj = 0
#    )
#    
#    if (nrow(qdat) > 0) {
#      graphics::mtext(
#        sprintf("%s  q: %.3g", unique(qdat$Type)[1], unique(qdat$survey_q)[1]),
#        side = 1, line = 2.5, cex = dyn_cex, adj = 0
#      )
#    }
#  }
#  
#  sp_index <- which(Rsim.scenario$params$spname == species)
#  
#  if (length(sp_index) > 0 && !is.na(Rsim.scenario$params$B_BaseRef[sp_index])) {
#    graphics::abline(h = Rsim.scenario$params$B_BaseRef[sp_index], col = "darkred", lty = 3)
#  }
#  
#  graphics::points(as.numeric(qdat$Year), mn)
#  graphics::segments(as.numeric(qdat$Year), y0 = up, y1 = dn)
#}
#
# #######################################
# #' Plot Catch Fit for a Single Species (Small Panel)
# #' 
# #' Generates a compact diagnostic plot comparing simulated catch against observed
# #' catch data for a specific species. Designed specifically to be nested inside
# #' multi-panel plot grids.
# #'
# #' @details
# #' This function extracts goodness-of-fit metrics from `rsim.fit.obj()` specifically
# #' for the `Catch` component. It calculates log-normal confidence intervals for the
# #' observed catch data and overlays the simulation's continuous catch estimates.
# #' A red dashed reference line is also drawn to represent the base historical catch
# #' derived from the scene's fishing effort and base biomass reference points.
# #'
# #'
# #' @inheritParams rsim.plot.fit
# #' @param species A single character string representing the species name.
# #'
# #' @export
#rsim.plot.fitcatch <- function(Rsim.scenario, Rsim.output, species) {
#  
#  if (length(species) > 1) {
#    warning("rsim.plot.fitcatch expects a single species string. Using only the first element.")
#    species <- species[1]
#  }
#  
#  if(!is.null(Rsim.scenario$fitting$Catch)){
#    catch.obj <- rsim.fit.obj(Rsim.scenario, Rsim.output)$Catch
#  } else {
#    catch.obj <- data.frame(Year = character(), Group = character(), Stdev = double(),
#                            Value = double(), Scale = double(), obs = double(),
#                            sd = double(), wt = double())    
#  }
#  
#  qdat <- catch.obj[catch.obj$Group == species, ]
#  
#  #if (nrow(qdat) == 0) {
#  #  warning(sprintf("No catch observation data found for species '%s'.", species))
#  #}
#  
#  oldpar <- par(mar = c(3, 2, 2, 1))
#  on.exit(par(oldpar))
#  
#  pin <- graphics::par("pin")
#  dyn_cex <- max(0.75, min(0.95, min(pin) * 0.2 + 0.6)) #FLAG ####
#  
#  mn   <- qdat$obs
#  sdlog <- sqrt(log(1.0 + (qdat$sd * qdat$sd) / ifelse(mn == 0, 1e-6, mn * mn)))
#  
#  up <- mn * exp(1.96 * sdlog)
#  dn <- mn / exp(1.96 * sdlog)
#  est  <- Rsim.output$annual_Catch[, species]
#  tot <- sum(qdat$fit * qdat$wt)
#  
#  all_vals <- c(up, est)
#  valid_vals <- all_vals[!is.na(all_vals) & !is.infinite(all_vals)]
#  ymax <- if (length(valid_vals) > 0) max(valid_vals) else 1
#  if (ymax <= 0) ymax <- 1
#  
#  graphics::plot(
#    as.numeric(rownames(Rsim.output$annual_Catch)),
#    est, type = "l", ylim = c(0, ymax), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n"
#  )
#  
#  graphics::axis(1, mgp = c(3, 0.4, 0), tck = -0.04, cex.axis = dyn_cex) #FLAG####
#  graphics::axis(2, mgp = c(3, 0.4, 0), tck = -0.04, cex.axis = dyn_cex) #FLAG####

#  #graphics::mtext(paste(species, sprintf("   nll: %.3g", tot)),
#  #               side = 1, cex = dyn_cex*0.8, line = 1.5, adj = 0)
#  
#  graphics::mtext(species, side = 1, cex = dyn_cex*0.8, line = 1.0, adj = 0)
#  graphics::mtext(sprintf("nll: %.3g", tot), side = 1, cex = dyn_cex*0.8, line = 1.8, adj = 0)
#  
#  graphics::points(as.numeric(qdat$Year), mn)
#  graphics::segments(as.numeric(qdat$Year), y0 = up, y1 = dn)
#  
#  slist <- Rsim.scenario$params$spname[Rsim.scenario$params$FishFrom + 1]
#  tcatch <- sum((Rsim.scenario$params$FishQ * Rsim.scenario$params$B_BaseRef[slist])
#                [slist == species], na.rm = TRUE)
#  if (!is.na(tcatch) && tcatch > 0)
#    abline(h = tcatch, col = "darkred", lty = 3)
#}
#


# ################################################################################
# #' Plot Simulated vs. Observed Catch for a Single Species
# #'
# #' Generates a standard plot comparing the annual simulated catch against historical
# #' catch observations for a specific species, complete with symmetric confidence bounds.
# #'
# #' @details
# #' Unlike the other diagnostic plotting functions that rely on `rsim.fit.obj()`,
# #' this function extracts the historical catch data directly from the scene object
# #' (`scene$fitting$Catch`). Furthermore, it calculates the upper and lower bounds
# #' of the confidence interval using a simple symmetric standard normal calculation
# #' (`mn +/- 1.96 * sd`) rather than assuming a log-normal error distribution.
# #'
# #' @inheritParams rsim.plot.fit 
# #'
# #' @return Invoked for its side effect of drawing a plot to the active graphics
# #'   device. Returns `NULL`.
# #'
# #'@export
# rsim.plot.catch <- function(Rsim.scenario, Rsim.output, species) {
#   qdat <- Rsim.scenario$fitting$Catch[Rsim.scenario$fitting$Catch$Group == species, ]
#   mn   <- qdat$obs
#   up   <- mn + 1.96 * qdat$sd
#   dn   <- mn - 1.96 * qdat$sd
#   est <- Rsim.output$annual_Catch[, species]
#   
#   #tot <- 0 #sum(qdat$fit)
#   ymax <- max(c(up, est), na.rm = TRUE)
#   
#   if (is.infinite(ymax) || is.na(ymax))
#     ymax <- 1
#   
#   graphics::plot(
#     as.numeric(rownames(Rsim.output$annual_Catch)),
#     est ,
#     type = "l",
#     ylim = c(0, ymax),
#     xlab = "Year",
#     ylab = ""
#   )
#   graphics::mtext(
#     side = 2,
#     line = 2.2,
#     paste(species, "catch"),
#     font = 2,
#     cex = 1.0
#   )
#   
#   if (nrow(qdat) > 0) {
#     graphics::points(as.numeric(qdat$Year), mn)
#     graphics::segments(as.numeric(qdat$Year), y0 = up, y1 = dn)
#   }
# }

#  ################################################################################
#  #' Render and Display Rpath Simulation Fit Report
#  #'
#  #' Compiles an R Markdown report (`fit_display.Rmd`) summarizing model fitting
#  #' outputs and opens the resulting HTML report in the RStudio Viewer pane.
#  #'
#  #' @details
#  #' This function acts as a wrapper around `rmarkdown::render()`. It passes the
#  #' balanced model, simulation scene, and fit run objects directly into the
#  #' parameterized R Markdown file `fit_display.Rmd`. The rendered output is stored
#  #' inside an `html/` directory in the current working directory.
#  #' 
#  #' @param bal An `Rpath` mass-balanced model object.
#  #' @param Rsim.scenario An Rsim scenario object containing the fitting data.
#  #' @param Rsim.output An `Rsim` simulation output object resulting from a fitted run.
#  #' @param output A character string defining the base name of the generated HTML file
#  #'   (without the `.html` extension). Defaults to `"test"`.
#  #'
#  #' @return Invoked primarily for its side effect of generating an HTML report file
#  #'   and launching it in the RStudio Viewer pane.
#  #'
#  #'@export
# render.show.fit <- function(bal, Rsim.scenario, Rsim.output, output = "test") {
#   wdir <- file.path(getwd(), "html")
#   
#   if (!dir.exists(wdir))
#     dir.create(wdir, recursive = TRUE)
#   
#   output_clean <- tools::file_path_sans_ext(output)
#   out_name <- paste0(output_clean, ".html")
#   
#   
#   rmd_file <- "fit_display.Rmd"
#   if (!file.exists(rmd_file)) {
#     pkg_file <- system.file("rmd", "fit_display.Rmd", package = "Rpath")
#     if (pkg_file != "")
#       rmd_file <- pkg_file
#   }
#   
#   out_path <- rmarkdown::render(
#     rmd_file,
#     output_file = out_name,
#     output_dir = wdir,
#     params = list(
#       bal = bal,
#       Rsim.scenario = Rsim.scenario,
#       Rsim.output = Rsim.output
#     )
#   )
#   
#   viewer <- getOption("viewer")
#   if (!is.null(viewer)) {
#     viewer(out_path)
#   } else{
#     #FLAG ####
#     #doble check if we have utils as a dependency.
#     utils::browseURL(out_path)
#   }
# }

# #################################################################################
# #' Plot Relative Biomass Over Time
# #'
# #' Generates a line plot showing the relative change in biomass over time for one
# #' or multiple species from an Rpath simulation run. Biomass is scaled relative
# #' to the starting biomass (month 1) of each species.
# #'
# #' @details
# #' This function extracts the monthly biomass data (`out_Biomass`) from an `Rsim`
# #' output object. Depending on the `indplot` flag, it calculates the relative
# #' trajectory for either a single species or a group of species. It dynamically
# #' calculates the legend width to place the legend entirely outside the right
# #' margin of the plotting area.
# #'
# #' @inheritParams rsim.plot.fit
# #' @param spname A character vector of species names to plot. If `indplot = TRUE`,
# #'   this should be a single character string.
# #' @param indplot A logical flag (`TRUE`/`FALSE` or `T`/`F`). If `FALSE` (default),
# #'   `spname` is treated as a vector of column names. If `TRUE`, `spname` is treated
# #'   as a single species name to look up in `params$spname`.
# #' @param ... Additional graphical parameters passed to the `plot()` setup function.
# #'
# #' @return Invoked for its side effect of drawing a plot to the active graphics
# #'   device. Returns `NULL`.
# #'
# #'@export
# rsim.plot.ylim <- function(Rsim.output, spname, indplot = FALSE, ...) {
#   oldpar <- graphics::par(no.readonly = TRUE)
#   on.exit(graphics::par(oldpar))
#   
#   if (indplot == FALSE) {
#     biomass <- Rsim.output$out_Biomass[, spname, drop = FALSE]
#     n <- ncol(biomass)
#     start.bio <- biomass[1, ]
#     start.bio[which(start.bio == 0)] <- 1
#     rel.bio <- matrix(NA, dim(biomass)[1], dim(biomass)[2])
#     for (isp in 1:n)
#       rel.bio[, isp] <- biomass[, isp] / start.bio[isp]
#   }
#   if (indplot == TRUE) {
#     spnum <- which(Rsim.output$params$spname == spname)
#     biomass <- Rsim.output$out_Biomass[, spnum, drop = FALSE]
#     n <- 1
#     st_val <- ifelse(biomass[1] == 0, 1, biomass[1])
#     rel.bio <- biomass / st_val
#   }
#   
#   ymax <- max(rel.bio, na.rm = TRUE) + 0.1 * max(rel.bio, na.rm = TRUE)
#   ymin <- min(rel.bio, na.rm = TRUE) - 0.1 * min(rel.bio, na.rm = TRUE)
#   
#   #xmax <- if(indplot) length(biomass) else nrow(biomass) BD way
#   ifelse(indplot, xmax <- length(biomass), xmax <- nrow(biomass))
#   
#   #Plot relative biomass
#   graphics::par(mar = c(4, 6, 2, 0))
#   
#   line.col <- grDevices::rainbow(n)
#   #Create space for legend
#   graphics::plot.new()
#   l <- graphics::legend(
#     0,
#     0,
#     bty = 'n',
#     spname,
#     plot = FALSE,
#     fill = line.col,
#     cex = 0.6
#   )
#   # calculate right margin width in ndc
#   w <- graphics::grconvertX(l$rect$w, to = 'ndc') - graphics::grconvertX(0, to =
#                                                                            'ndc')
#   
#   graphics::par(omd = c(0, 1 - w, 0, 1))
#   
#   graphics::plot(
#     0,
#     0,
#     xlim = c(0, xmax),
#     ylim = c(ymin, ymax),
#     axes = FALSE,
#     xlab = '',
#     ylab = '',
#     type = 'n',
#     ...
#   )
#   graphics::axis(1)
#   graphics::axis(2, las = TRUE)
#   graphics::box(lwd = 2)
#   graphics::mtext(1,
#                   text = 'Months',
#                   line = 2.5,
#                   cex = 1.8)
#   graphics::mtext(2,
#                   text = 'Relative Biomass',
#                   line = 3,
#                   cex = 1.8)
#   
#   
#   for (i in 1:n) {
#     if (indplot == TRUE)
#       graphics::lines(rel.bio, col = line.col[i], lwd = 3)
#     if (indplot == FALSE)
#       graphics::lines(rel.bio[, i], col = line.col[i], lwd = 3)
#   }
#   
#   graphics::legend(
#     graphics::par('usr')[2],
#     graphics::par('usr')[4],
#     bty = 'n',
#     xpd = NA,
#     spname,
#     fill = line.col,
#     cex = 0.6
#   )
#   
#   graphics::par(oldpar)
# }