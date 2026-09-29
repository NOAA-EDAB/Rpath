################################################################################
#' Render and Display Rpath Simulation Fit Report
#'
#' Compiles an R Markdown report (`fit_display.Rmd`) summarizing model fitting
#' outputs and opens the resulting HTML report in the RStudio Viewer pane.
#'
#' @details
#' This function acts as a wrapper around `rmarkdown::render()`. It passes the
#' balanced model, simulation scene, and fit run objects directly into the
#' parameterized R Markdown file `fit_display.Rmd`. The rendered output is stored
#' inside an `html/` directory in the current working directory.
#'
#' @param bal An `Rpath` mass-balanced model object.
#' @param scene_fit An `Rpath` simulation scene object containing fitting configuration.
#' @param run_fit An `Rsim` simulation output object resulting from a fitted run.
#' @param output A character string defining the base name of the generated HTML file
#'   (without the `.html` extension). Defaults to `"test"`.
#'
#' @return Invoked primarily for its side effect of generating an HTML report file
#'   and launching it in the RStudio Viewer pane.
#'

render.show.fit <- function(bal, scene_fit, run_fit, output = "test") {
  wdir <- file.path(getwd(), "html")
  
  if (!dir.exists(wdir))
    dir.create(wdir, recursive = TRUE)
  
  output_clean <- tools::file_path_sans_ext(output)
  out_name <- paste0(output_clean, ".html")
  
  
  rmd_file <- "fit_display.Rmd"
  if (!file.exists(rmd_file)) {
    pkg_file <- system.file("rmd", "fit_display.Rmd", package = "Rpath")
    if (pkg_file != "")
      rmd_file <- pkg_file
  }
  
  out_path <- rmarkdown::render(
    rmd_file,
    output_file = out_name,
    output_dir = wdir,
    params = list(
      bal = bal,
      scene_fit = scene_fit,
      run_fit = run_fit
    )
  )
  
  viewer <- getOption("viewer")
  if (!is.null(viewer)) {
    viewer(out_path)
  } else{
    #FLAG ####
    #doble check if we have utils as a dependency.
    utils::browseURL(out_path)
  }
}

################################################################################
#' Plot Biomass and Catch Fitting Results
#'
#' Generates multi-panel diagnostic plots comparing simulation run outputs against
#' observed biomass and catch data across a specified list of species.
#'
#' @details
#' This function opens two separate graphic windows. The first displays biomass
#' fits, iterating through all available data sources for each specified species.
#' The second window displays catch fits. Both plots include an overarching title
#' that dynamically calculates and reports the negative log-likelihood (NLL) of
#' the fit.
#'
#' @param scene An `Rpath` simulation scene list object containing fitting parameters
#'   and observation series.
#' @param run An `Rsim` simulation output object (the result of a fitted run).
#' @param species A character vector of species names to include in the plots.
#' @param scene_name An optional character string used to label the main title
#'   of the plots. Defaults to `NULL` (resulting in a blank space).
#'
#' @return Invoked for its side effect of drawing plots to external graphics devices.
#'   Returns `NULL`.
#'

rsim.runplot <- function(scene, run, species, scene_name = NULL, small=TRUE) {
  if (is.null(scene_name)) {
    display_name <- " "
  } else {
    display_name <- scene_name
  }
  
  run.out <- run #rsim.fit.run(values,groups,vartypes,scene,"AB",years,T)
  
  #all_sources <- rsim.fit.list.bio.series(scene_fit)$sources
  all_series <- strsplit(rsim.fit.list.bio.series(scene)$all, ":")
  obs_species <- sapply(all_series, function(x)
    x[1]) # "[[", 1)
  obs_sources <- sapply(all_series, function(x)
    if (length(x) > 1)
      x[2]
    else
      NA) #"[[", 2)
  
  nll_val <- round(rsim.fit.obj(scene, run, FALSE), 2)
  

  total_bio_panels <- sum(sapply(species, function(sp)
    max(1, sum(
      obs_species == sp, na.rm = TRUE
    ))))
  nc_bio <- min(11, ceiling(sqrt(total_bio_panels))) #columns
  nr_bio <- ceiling(total_bio_panels / nc_bio) #row

  dev.new()

  graphics::par(mfrow = c(nr_bio, nc_bio), oma = c(0, 0, 3, 0)) #oma= outside of margin area
  
  for (i in 1:length(species)) {
    sources <- obs_sources[obs_species == species[i] &
                             !is.na(obs_sources)]
    if (length(sources) == 0) {
      rsim.plot.fitbio(scene, run.out, species[i], NA, small=small)
    } else{
      for (j in 1:length(sources)) {
        rsim.plot.fitbio(scene, run.out, species[i], sources[j], small =small)
      }
    }
  }
  graphics::mtext(
    paste("Fitting results for:", display_name, paste(", nll:", nll_val)),
    outer = TRUE,
    side = 3,
    line = 1,
    cex = 1.5,
    font = 1
  )
  
  #catch layout
  nc_catch <- ceiling(sqrt(length(species)))
  nr_catch <- ceiling(length(species) / nc_catch)
  
  
  dev.new()
  graphics::par(mfrow = c(nr_catch, nc_catch),
                oma = c(0, 0, 3, 0))
  for (sp in species) {
    rsim.plot.fitcatch.small(scene, run.out, sp)
  }
  graphics::mtext(
    paste(
      "Catch fitting results for:",
      display_name,
      paste(", nll:", nll_val)
    ),
    outer = TRUE,
    side = 3,
    line = 1,
    cex = 1.5,
    font = 1
  )
  #return(run.out)
  
}

################################################################################
#' Plot Biomass Fit for a Single Species and Data Source
#'
#' Generates a diagnostic plot comparing simulated biomass against observed
#' survey data for a specific species and a specific data source. It includes
#' log-normal confidence intervals for the observed data and displays the NLL score.
#'
#' @details
#' This function extracts the goodness-of-fit metrics from `rsim.fit.obj()` for a
#' given species and observation source. It calculates the mean observed biomass
#' adjusted by survey catchability (`survey_q`) and computes the upper and lower
#' bounds of the 95% confidence interval assuming a log-normal error distribution.
#' The simulated biomass is plotted as a continuous line, with observed data as
#' points and error bars.
#'
#' @param scene An `Rpath` simulation scene object containing the base reference
#'   parameters and data fitting configurations.
#' @param run An `Rsim` simulation output object containing the annual estimated biomass.
#' @param species A character string representing the name of the species/group to plot.
#' @param datasource A character string identifying the specific survey or data source
#'   to extract from the objective function data (e.g., `"race_wgoa"`).
#'
#' @return Invoked for its side effect of drawing a plot to the active graphics
#'   device. Returns `NULL`.
#'

rsim.plot.fitbio <- function(scene, run, species, datasource, small= FALSE) {
  
  if(!is.null(scene$fitting$Biomass)){
    bio.obj <- rsim.fit.obj(scene, run)$Biomass
  } else {
    bio.obj <- data.frame(Year  = character(), 
                          Group = character(), 
                          Type  = character(), 
                          Stdev = double(),
                          Value = double(), 
                          Scale = double(),
                          Source = character(),
                          ID     = character(),
                          obs    = double(),
                          sd     = double(),
                          initial_q = double(),
                          wt        = double())
  }
  
  qdat <- bio.obj[bio.obj$Group == species &
                    bio.obj$Source %in% datasource, ]
  
  if(small){
    mar_vals <- c(2.8, 1.5, 1, 0.5)
    axis1_mgp <- c(3,0,0)
    axis2_mgp <- c(3,0.2,0)
  } else {
    mar_vals <- c(4, 3, 2, 1)
    axis1_mgp <- c(3, 0.6, 0)
    axis2_mgp <- c(3, 0.6, 0)
  }
  
  oldpar <- graphics::par(mar = mar_vals)
  on.exit(graphics::par(oldpar))
  
  # zero/NA catchability fix the division by 0
  sq <- qdat$survey_q
  sq[sq == 0 | is.na(sq)] <- 1e-6
  
  mn   <- qdat$obs / sq
  #FLAG ####
  sdlog <- sqrt(log(1.0 + (qdat$sd / sq) * (qdat$sd / sq) / ifelse(mn ==
                                                                     0, 1e-6, mn * mn))) # review if this is right
  up <- mn * exp(1.96 * sdlog)
  dn <- mn / exp(1.96 * sdlog)
  #up   <- mn + 1.96*qdat$sd / qdat$survey_q #/survey_q
  #dn   <- mn - 1.96*qdat$sd / qdat$survey_q #/survey_q
  
  est <- run$annual_Biomass[, species] # * qdat$survey_q
  tot  <- sum(qdat$fit * qdat$wt)
  
  all_vals <- c(up, est)
  valid_vals <- all_vals[!is.na(all_vals) & !is.infinite(all_vals)]
  ymax <- if(length(valid_vals) > 0) max(valid_vals) else 1
  if (ymax <= 0) ymax <- 1
  
  
  graphics::plot(
    as.numeric(rownames(run$annual_Biomass)),
    est,
    type = "l",
    ylim = c(0, ymax),
    xaxt = "n",
    yaxt = "n",
    xlab = "",
    ylab = "",
    bty = "n"
  )
  graphics::axis(1,
                 mgp = axis1_mgp,
                 tck = -0.04,
                 cex.axis = 0.6)
  graphics::axis(2,
                 mgp = axis2_mgp,
                 tck = -0.04,
                 cex.axis = 0.8)
  
   
  if(small){
    graphics::mtext(species, side = 1, cex = 0.6, line = 0.9, adj = 0)
    graphics::mtext(datasource, side = 1, cex = 0.6, line = 1.7, adj = 0)
    
    if (nrow(qdat) > 0) {
      graphics::mtext(
        sprintf("nl:%.3g  q:%.3g (%s)", tot, unique(sq)[1], unique(qdat$Type)[1]),
        side = 1,
        line = 2.5,
        cex = 0.6,
        adj = 0
      )
    }
  } else {
  graphics::mtext(
    paste(datasource, species, "biomass", sprintf("   nll: %.3g", tot)),
    side = 1,
    cex = 0.9,
    line = 1.5,
    adj = 0
  )
  
  if (nrow(qdat) > 0) {
    graphics::mtext(
      sprintf("%s  q: %.3g", unique(qdat$Type)[1], unique(qdat$survey_q)[1]),
      side = 1,
      line = 2.5,
      cex = 0.9,
      adj = 0
    )
  }
}
  sp_index <- which(scene$params$spname == species)
  
  if (length(sp_index) > 0 &&
      !is.na(scene$params$B_BaseRef[sp_index])) {
    graphics::abline(
      h = scene$params$B_BaseRef[sp_index],
      col = "darkred",
      lty = 3
    )
  }
  
  graphics::points(as.numeric(qdat$Year), mn)
  graphics::segments(as.numeric(qdat$Year), y0 = up, y1 = dn)
}


#######################################
#' Plot Catch Fit for a Single Species (Small Panel)
#'
#' Generates a compact diagnostic plot comparing simulated catch against observed
#' catch data for a specific species. Designed specifically to be nested inside
#' multi-panel plot grids.
#'
#' @details
#' This function extracts goodness-of-fit metrics from `rsim.fit.obj()` specifically
#' for the `Catch` component. It calculates log-normal confidence intervals for the
#' observed catch data and overlays the simulation's continuous catch estimates.
#' A red dashed reference line is also drawn to represent the base historical catch
#' derived from the scene's fishing effort and base biomass reference points.
#'
#' @param scene An `Rpath` simulation scene object containing the base reference
#'   parameters (`FishQ`, `B_BaseRef`, `spname`) and data fitting configurations.
#' @param run An `Rsim` simulation output object containing the annual estimated catch.
#' @param species A character string representing the name of the species/group to plot.
#'
#' @return Invoked for its side effect of drawing a plot to the active graphics
#'   device. Returns `NULL`.
#'

rsim.plot.fitcatch.small <- function(scene, run, species) {
  
  if(!is.null(scene$fitting$Catch)){
     catch.obj <- rsim.fit.obj(scene, run)$Catch
  } else {
    catch.obj <- data.frame(Year  = character(), 
                          Group = character(), 
                          Stdev = double(),
                          Value = double(), 
                          Scale = double(),
                          obs    = double(),
                          sd     = double(),
                          wt        = double())    
  }
  
  qdat <- catch.obj[catch.obj$Group == species, ]
  
  oldpar <- par(mar = c(3, 2, 2, 1))
  on.exit(par(oldpar))
  
  mn   <- qdat$obs
  sdlog <- sqrt(log(1.0 + (qdat$sd * qdat$sd) / ifelse(mn == 0, 1e-6, mn *
                                                         mn)))
  
  up <- mn * exp(1.96 * sdlog)
  dn <- mn / exp(1.96 * sdlog)
  #up   <- mn + 1.96*qdat$sd
  #dn   <- mn - 1.96*qdat$sd
  est  <- run$annual_Catch[, species]
  tot <- sum(qdat$fit * qdat$wt)
  
  all_vals <- c(up, est)
  valid_vals <- all_vals[!is.na(all_vals) & !is.infinite(all_vals)]
  ymax <- if (length(valid_vals) > 0)
    max(valid_vals)
  else
    1
  if (ymax <= 0)
    ymax <- 1
  
  graphics::plot(
    as.numeric(rownames(run$annual_Catch)),
    est,
    type = "l",
    ylim = c(0, ymax),
    xaxt = "n",
    yaxt = "n",
    xlab = "",
    ylab = "",
    bty = "n"
  )
  axis(1,
       mgp = c(3, 0.6, 0),
       tck = -0.04,
       cex.axis = 0.6)
  axis(2,
       mgp = c(3, 0.2, 0),
       tck = -0.04,
       cex.axis = 0.8)
  graphics::mtext(
    paste(species, sprintf("   nll: %.3g", tot)),
    side = 1,
    cex = 0.6,
    line = 1.5,
    adj = 0
  )
  graphics::points(as.numeric(qdat$Year), mn)
  graphics::segments(as.numeric(qdat$Year), y0 = up, y1 = dn)
  
  
  slist <- scene$params$spname[scene$params$FishFrom + 1]
  tcatch <- sum((scene$params$FishQ * scene$params$B_BaseRef[slist])
                [slist == species], na.rm = TRUE)
  if (!is.na(tcatch) &&
      tcatch > 0)
    abline(h = tcatch,
           col = "darkred",
           lty = 3)
}

#################################################################################
#' Plot Full Biomass and Catch Fits for a Single Species
#'
#' Generates a multi-panel figure displaying the simulation fits for a single
#' species across all available biomass observation sources, alongside its catch
#' fitting data.
#'
#' @details
#' This function aggregates the biomass and catch visualization into a single row
#' of plots. It dynamically determines the number of biomass data sources for the
#' specified species and configures the graphics layout (`par(mfrow)`) to fit all
#' biomass plots plus one additional panel for the catch plot. It calculates
#' log-normal confidence intervals for observations and overlays base reference lines.
#'
#' @param scene An `Rpath` simulation scene object containing the base reference
#'   parameters and data fitting configurations.
#' @param run An `Rsim` simulation output object containing the annual estimated
#'   biomass and catch.
#' @param species A character string representing the name of the species/group to plot.
#'
#' @return Invoked for its side effect of drawing a multi-panel plot to the active
#'   graphics device. Returns `NULL`.
#'

rsim.plot.full <- function(scene, run, species) {
  oldpar <- graphics::par(no.readonly = TRUE)
  #no.readonly = logical; if TRUE and there are no other arguments, only parameters are returned which can be set by a subsequent par() call on the same device.
  on.exit(graphics::par(oldpar))
  
  # Biomass plotting
  # Remove TMP here - don't forget.
  fit_obj <- rsim.fit.obj(scene, run)
  bio.obj <- fit_obj$Biomass
  sdat <- bio.obj[bio.obj$Group == species, ]
  #sid <- paste(sdat$Source, sdat$Group, sep=":")
  sidset <- unique(sdat$Source)
  
  graphics::par(
    mfrow = c(1, length(sidset) + 1),
    oma = c(0.0, 0.0, 0.0, 0.0),
    mar = c(3, 2, 2, 1)
  )
  
  for (S in sidset) {
    qdat <- sdat[sdat$Source == S, ]
    sq <- qdat$survey_q
    sq[sq == 0 | is.na(sq)] <- 1e-6
    #survey_q <- 1
    mn   <- qdat$obs / sq
    sdlog <- sqrt(log(1.0 + (qdat$sd / sq) * (qdat$sd / sq) /
                        ifelse(mn == 0, 1e-6, mn * mn)))
    up <- mn * exp(1.96 * sdlog)
    dn <- mn / exp(1.96 * sdlog)
    #up   <- mn + 1.96*qdat$sd / qdat$survey_q #/survey_q
    #dn   <- mn - 1.96*qdat$sd / qdat$survey_q #/survey_q
    
    est <- run$annual_Biomass[, species]
    tot  <- sum(qdat$fit * qdat$wt)
    
    max_up <- max(up[!is.na(up)], na.rm = TRUE)
    max_est <- max(est[!is.na(est)], na.rm = TRUE)
    ymax <- max(c(max_up, max_est), na.rm = TRUE)
    if (is.infinite(ymax) ||
        is.na(ymax))
      ymax <- max(est, 1, na.rm = TRUE)
    
    graphics::plot(
      as.numeric(rownames(run$annual_Biomass)),
      est,
      type = "l",
      ylim = c(0, ymax),
      xaxt = "n",
      yaxt = "n",
      xlab = "",
      ylab = "",
      bty = "n"
    )
    axis(1,
         mgp = c(3, 0.0, 0),
         tck = -0.04,
         cex.axis = 0.6)
    axis(2,
         mgp = c(3, 0.2, 0),
         tck = -0.04,
         cex.axis = 0.8)
    mtext(
      paste(S, species, "B", sprintf("nll: %.3g", tot)),
      side = 1,
      cex.main = 0.9,
      line = 0.8,
      adj = 0
    )
    if (nrow(qdat) > 0) {
      mtext(
        sprintf("%s  q: %.3g", unique(qdat$Type)[1], unique(sq)[1]),
        side = 1,
        line = 1.5,
        cex = 0.9,
        adj = 0
      )
    }
    sp_index <- which(scene$params$spname == species)
    if (length(sp_index) > 0 &&
        !is.na(scene$params$B_BaseRef[sp_index])) {
      abline(
        h = scene$params$B_BaseRef[sp_index],
        col = "darkred",
        lty = 3
      )
    }
    
    points(as.numeric(qdat$Year), mn)
    segments(as.numeric(qdat$Year), y0 = up, y1 = dn)
  }
  
  # Catch plotting
  # REMOVE TMP here
  catch.obj <- fit_obj$Catch
  qdat <- catch.obj[catch.obj$Group == species, ]
  mn   <- qdat$obs
  sdlog <- sqrt(log(1.0 + (qdat$sd * qdat$sd) / ifelse(mn == 0, 1e-6, mn * mn)))
  up <- mn * exp(1.96 * sdlog)
  dn <- mn / exp(1.96 * sdlog)
  #up   <- mn + 1.96*qdat$sd
  #dn   <- mn - 1.96*qdat$sd
  est  <- run$annual_Catch[, species]
  tot <- sum(qdat$fit * qdat$wt)
  
  max_up <- max(up[!is.na(up)], na.rm = TRUE)
  max_est <- max(est[!is.na(est)], na.rm = TRUE)
  ymax <- max(c(max_up, max_est), na.rm = TRUE)
  if (is.infinite(ymax) ||
      is.na(ymax))
    ymax <- max(est, 1, na.rm = TRUE)
  
  
  graphics::plot(
    as.numeric(rownames(run$annual_Catch)),
    est,
    type = "l",
    ylim = c(0, ymax),
    xaxt = "n",
    yaxt = "n",
    xlab = "",
    ylab = "",
    bty = "n"
  )
  axis(1,
       mgp = c(3, 0.2, 0),
       tck = -0.04,
       cex.axis = 0.6)
  axis(2,
       mgp = c(3, 0.2, 0),
       tck = -0.04,
       cex.axis = 0.8)
  graphics::mtext(
    paste(species, "catch", sprintf("   nll: %.3g", tot)),
    side = 1,
    cex.main = 0.9,
    line = 0.8,
    adj = 0
  )
  
  graphics::points(as.numeric(qdat$Year), mn)
  graphics::segments(as.numeric(qdat$Year), y0 = up, y1 = dn)
  
  slist <- scene$params$spname[scene$params$FishFrom + 1]
  tcatch <- sum((scene$params$FishQ * scene$params$B_BaseRef)[slist == species], na.rm =
                  TRUE)
  if (!is.na(tcatch) &&
      tcatch > 0)
    abline(h = tcatch,
           col = "darkred",
           lty = 3)
}


#################################################################################
#' Plot Relative Biomass Over Time
#'
#' Generates a line plot showing the relative change in biomass over time for one
#' or multiple species from an Rpath simulation run. Biomass is scaled relative
#' to the starting biomass (month 1) of each species.
#'
#' @details
#' This function extracts the monthly biomass data (`out_Biomass`) from an `Rsim`
#' output object. Depending on the `indplot` flag, it calculates the relative
#' trajectory for either a single species or a group of species. It dynamically
#' calculates the legend width to place the legend entirely outside the right
#' margin of the plotting area.
#'
#' @param Rsim.output An `Rsim` simulation output object containing `out_Biomass`
#'   and base parameters (`params$spname`).
#' @param spname A character vector of species names to plot. If `indplot = TRUE`,
#'   this should be a single character string.
#' @param indplot A logical flag (`TRUE`/`FALSE` or `T`/`F`). If `FALSE` (default),
#'   `spname` is treated as a vector of column names. If `TRUE`, `spname` is treated
#'   as a single species name to look up in `params$spname`.
#' @param ... Additional graphical parameters passed to the `plot()` setup function.
#'
#' @return Invoked for its side effect of drawing a plot to the active graphics
#'   device. Returns `NULL`.
#'

rsim.plot.ylim <- function(Rsim.output, spname, indplot = FALSE, ...) {
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  
  if (indplot == FALSE) {
    # KYA April 2020 this seems incorrect? REVIEWD BD
    #biomass <- Rsim.output$out_Biomass[, 2:ncol(Rsim.output$out_Biomass)]
    biomass <- Rsim.output$out_Biomass[, spname, drop = FALSE]
    n <- ncol(biomass)
    start.bio <- biomass[1, ]
    start.bio[which(start.bio == 0)] <- 1
    rel.bio <- matrix(NA, dim(biomass)[1], dim(biomass)[2])
    for (isp in 1:n)
      rel.bio[, isp] <- biomass[, isp] / start.bio[isp]
  }
  if (indplot == TRUE) {
    spnum <- which(Rsim.output$params$spname == spname)
    biomass <- Rsim.output$out_Biomass[, spnum, drop = FALSE]
    n <- 1
    st_val <- ifelse(biomass[1] == 0, 1, biomass[1])
    rel.bio <- biomass / st_val
  }
  
  ymax <- max(rel.bio, na.rm = TRUE) + 0.1 * max(rel.bio, na.rm = TRUE)
  ymin <- min(rel.bio, na.rm = TRUE) - 0.1 * min(rel.bio, na.rm = TRUE)
  
  #xmax <- if(indplot) length(biomass) else nrow(biomass) BD way
  ifelse(indplot, xmax <- length(biomass), xmax <- nrow(biomass))
  
  #Plot relative biomass
  graphics::par(mar = c(4, 6, 2, 0))
  
  line.col <- grDevices::rainbow(n)
  #Create space for legend
  graphics::plot.new()
  l <- graphics::legend(
    0,
    0,
    bty = 'n',
    spname,
    plot = FALSE,
    fill = line.col,
    cex = 0.6
  )
  # calculate right margin width in ndc
  w <- graphics::grconvertX(l$rect$w, to = 'ndc') - graphics::grconvertX(0, to =
                                                                           'ndc')
  
  graphics::par(omd = c(0, 1 - w, 0, 1))
  
  graphics::plot(
    0,
    0,
    xlim = c(0, xmax),
    ylim = c(ymin, ymax),
    axes = F,
    xlab = '',
    ylab = '',
    type = 'n',
    ...
  )
  graphics::axis(1)
  graphics::axis(2, las = T)
  graphics::box(lwd = 2)
  graphics::mtext(1,
                  text = 'Months',
                  line = 2.5,
                  cex = 1.8)
  graphics::mtext(2,
                  text = 'Relative Biomass',
                  line = 3,
                  cex = 1.8)
  
  
  for (i in 1:n) {
    if (indplot == TRUE)
      graphics::lines(rel.bio, col = line.col[i], lwd = 3)
    if (indplot == FALSE)
      graphics::lines(rel.bio[, i], col = line.col[i], lwd = 3)
  }
  
  graphics::legend(
    graphics::par('usr')[2],
    graphics::par('usr')[4],
    bty = 'n',
    xpd = NA,
    spname,
    fill = line.col,
    cex = 0.6
  )
  
  graphics::par(oldpar)
}

################################################################################
#' Plot Simulated vs. Observed Catch for a Single Species
#'
#' Generates a standard plot comparing the annual simulated catch against historical
#' catch observations for a specific species, complete with symmetric confidence bounds.
#'
#' @details
#' Unlike the other diagnostic plotting functions that rely on `rsim.fit.obj()`,
#' this function extracts the historical catch data directly from the scene object
#' (`scene$fitting$Catch`). Furthermore, it calculates the upper and lower bounds
#' of the confidence interval using a simple symmetric standard normal calculation
#' (`mn +/- 1.96 * sd`) rather than assuming a log-normal error distribution.
#'
#' @param scene An `Rpath` simulation scene object containing the historical catch
#'   fitting data under `scene$fitting$Catch`.
#' @param run An `Rsim` simulation output object containing the matrix of annual
#'   catch estimates (`annual_Catch`).
#' @param species A character string representing the name of the species/group to plot.
#'
#' @return Invoked for its side effect of drawing a plot to the active graphics
#'   device. Returns `NULL`.
#'
#'@export
rsim.plot.catch <- function(scene, run, species) {
  qdat <- scene$fitting$Catch[scene$fitting$Catch$Group == species, ]
  mn   <- qdat$obs
  up   <- mn + 1.96 * qdat$sd
  dn   <- mn - 1.96 * qdat$sd
  est <- run$annual_Catch[, species]
  
  #tot <- 0 #sum(qdat$fit)
  ymax <- max(c(up, est), na.rm = TRUE)
  
  if (is.infinite(ymax) || is.na(ymax))
    ymax <- 1
  
  graphics::plot(
    as.numeric(rownames(run$annual_Catch)),
    est ,
    type = "l",
    ylim = c(0, ymax),
    xlab = "Year",
    ylab = ""
  )
  graphics::mtext(
    side = 2,
    line = 2.2,
    paste(species, "catch"),
    font = 2,
    cex = 1.0
  )
  
  if (nrow(qdat) > 0) {
    graphics::points(as.numeric(qdat$Year), mn)
    graphics::segments(as.numeric(qdat$Year), y0 = up, y1 = dn)
  }
}

################################################################################
#' Plot Simulated vs. Observed Biomass for a Single Species
#'
#' Generates a standard diagnostic plot comparing the annual simulated biomass
#' against historical survey observations for a specific species, complete with
#' symmetric confidence bounds.
#'
#' @details
#' This function extracts goodness-of-fit metrics from `rsim.fit.obj()` for the
#' `Biomass` component. Unlike `rsim.plot.fitbio` (which calculates log-normal
#' errors), this function calculates symmetric confidence intervals directly using
#' standard error scaled by catchability (`mn +/- 1.96 * sd / survey_q`). It also
#' draws a red dashed reference line representing the species' initial biomass
#' (`B_Initial`).
#'
#' @param scene An `Rpath` simulation scene object containing the base reference
#'   parameters (specifically `B_Initial`) and data fitting configurations.
#' @param run An `Rsim` simulation output object containing the matrix of annual
#'   biomass estimates (`annual_Biomass`).
#' @param species A character string representing the name of the species/group to plot.
#'
#' @return Invoked for its side effect of drawing a plot to the active graphics
#'   device. Returns `NULL`.
#'
#'@export
rsim.plot.biomass <- function(scene, run, species) {
  bio.obj <- rsim.fit.obj(scene, run)$Biomass
  qdat <- bio.obj[bio.obj$Group == species, ]
  sq <- qdat$survey_q
  sq[sq == 0 | is.na(sq)] <- 1e-6
  
  #survey_q <- 1
  mn   <- qdat$obs / sq
  up   <- mn + 1.96 * qdat$sd / sq
  dn   <- mn - 1.96 * qdat$sd / sq
  tot  <- sum(qdat$fit)
  
  est <- run$annual_Biomass[, species]
  
  ymax <- max(c(up, est), na.rm = TRUE)
  if (is.infinite(ymax) || is.na(ymax))
    ymax <- 1
  
  graphics::plot(
    as.numeric(rownames(run$annual_Biomass)),
    est,
    type = "l",
    ylim = c(0, ymax),
    xlab = "",
    ylab = ""
  )
  graphics::mtext(
    side = 2,
    line = 2.2,
    paste(species, "biomass"),
    font = 2,
    cex = 0.8
  )
  #cat(species, tot, "\n");
  if (nrow(qdat) > 0) {
    q_str <- paste(unique(round(sq, 3)), collapse = ", ")
    graphics::mtext(
      sprintf("NLL: %.3g", tot),
      side = 1,
      line = 2,
      cex = 0.8
    )
    graphics::mtext(paste("  q:", q_str),
                    side = 1,
                    line = 3,
                    cex = 0.8)
    
    x_years <- as.numeric(as.character(qdat$Year))
    
    graphics::points(x_years, mn)
    graphics::segments(
      x0 = x_years,
      y0 = up,
      x1 = x_years,
      y1 = dn
    )
    
  }
  
  sp_index <- which(scene$params$spname == species)
  
  if (length(sp_index) > 0) {
    # Extract the value using the first match (in case of duplicates)
    b_ref <- scene$params$B_BaseRef[sp_index[1]] #changed it from B_Initial to B_BaseRef
    
    # Strictly check that the value exists, is a single number, and is not NA
    if (!is.null(b_ref) && length(b_ref) == 1 && !is.na(b_ref)) {
      graphics::abline(h = b_ref,
                       col = "darkred",
                       lty = 3)
    }
  }
}
