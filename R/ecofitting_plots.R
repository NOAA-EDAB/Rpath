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
#' @param Rsim.scenario An `Rsim` scenario object containing the fitting data.
#' @param Rsim.output An `Rsim` simulation output object resulting from a fitted run.
#' @param species.list A character vector of species names to include in the plots.
#' @param scene.name An optional character string used to label the main title
#'   of the plots. Defaults to `NULL` (resulting in a blank space).
#'
#' @return Invoked for its side effect of drawing plots to external graphics devices.
#'   Returns `NULL`.
#'
#'@export
rsim.plot.fit <- function(Rsim.scenario, Rsim.output, species.list, scene.name = NULL) {
  if(!is.character(species.list)){
    stop("Argument 'species.list' must be a vector.")
  }
  
  run.out <- Rsim.output
  
  missing_sp <- setdiff(species.list, colnames(run.out$annual_Biomass))
  if(length(missing_sp) > 0){
    warning(sprintf("The following species were not found in the simulation output and will be skipped: %s. Please check spelling.",
                    paste(missing_sp, collapse = ", ")))
    species.list <- intersect(species.list, colnames(run.out$annual_Biomass))
    if(length(species.list) == 0) stop("No valid species remain to plot.")
  }
  
  if(is.null(scene.name)) {
    display_name <- " "
  } else {
    display_name <- scene.name
  }
  
  all_series <- strsplit(rsim.fit.list.bio.series(Rsim.scenario)$all, ":")
  obs_species <- sapply(all_series, function(x) x[1]) 
  obs_sources <- sapply(all_series, function(x) if (length(x) > 1) x[2] else NA) 
  
  nll_val <- round(rsim.fit.obj(Rsim.scenario, Rsim.output, FALSE), 2)
  
  total_bio_panels <- sum(sapply(species.list, function(sp) max(1, sum(obs_species == sp, na.rm = TRUE))))
  nc_bio <- min(11, ceiling(sqrt(total_bio_panels))) #columns
  nr_bio <- ceiling(total_bio_panels / nc_bio) #row
  
  #fetch biomass obj data, this is part of the old rsim.plot.fitbio
  
  if(!is.null(Rsim.scenario$fitting$Biomass)){
    bio.obj <- rsim.fit.obj(Rsim.scenario, Rsim.output)$Biomass
  } else {
    bio.obj <- data.frame(Year = character(), Group = character(), Type = character(), 
                          Stdev = double(), Value = double(), Scale = double(),
                          Source = character(), ID = character(), obs = double(), 
                          sd = double(), initial_q = double(), wt = double())
  }
  
  
  #fetch catch obj data, same thing as above, but now for catch
  
  if(!is.null(Rsim.scenario$fitting$Catch)){
    catch.obj <- rsim.fit.obj(Rsim.scenario, Rsim.output)$Catch
  } else {
    catch.obj <- data.frame(Year = character(), Group = character(), Stdev = double(),
                            Value = double(), Scale = double(), obs = double(),
                            sd = double(), wt = double())   
  }
  
  #Helper for dynamic window plotting (same as the old functions)
  get_par_settings <- function() {
    pin <- graphics::par("pin") 
    dyn_cex <- max(0.75, min(0.95, min(pin) * 0.2 + 0.6))
      list(dyn_cex = dyn_cex, mar_vals = c(2.8, 1.5, 1, 0.5),
           axis1_mgp = c(3, 0.3, 0), axis2_mgp = c(3, 0.2, 0))
  }
  
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  
  #BIOMASS PLOTS
  
  grDevices::dev.new(noRStudioGD=TRUE)
  graphics::par(mfrow = c(nr_bio, nc_bio), 
                oma = c(0, 0, 3, 0)) #oma= outside of margin area

  
  for(sp in species.list) {
    sources <- obs_sources[obs_species == sp & !is.na(obs_sources)]
    plot_sources <- if (length(sources) == 0) NA else sources
    
    for(datasource in plot_sources) {
      qdat <- bio.obj[bio.obj$Group == sp & bio.obj$Source %in% datasource, ]
      
      #FLAG ### 
      # decide if we want this warnings, they are very annoying!! BD
      #if (nrow(qdat) == 0 && !is.na(datasource)) {
      #  warning(sprintf("No biomass observation data found for species '%s' and datasource '%s'.", sp, datasource))
      #}
      
      p_set <- get_par_settings()
      graphics::par(mar = p_set$mar_vals)
      
      sq <- qdat$survey_q
      sq[sq == 0 | is.na(sq)] <- 1e-6
      
      mn    <- qdat$obs / sq
      sdlog <- sqrt(log(1.0 + (qdat$sd / sq)^2 / ifelse(mn == 0, 1e-6, mn^2))) 
      up <- mn * exp(1.96 * sdlog)
      dn <- mn / exp(1.96 * sdlog)
      
      est <- run.out$annual_Biomass[, sp] 
      tot <- sum(qdat$fit * qdat$wt)
      
      all_vals <- c(up, est)
      valid_vals <- all_vals[!is.na(all_vals) & !is.infinite(all_vals)]
      ymax <- if(length(valid_vals) > 0) max(valid_vals) else 1
      if(ymax <= 0) ymax <- 1
      
      graphics::plot(as.numeric(rownames(run.out$annual_Biomass)), est, 
                     type = "l", ylim = c(0, ymax), xaxt = "n", yaxt = "n", 
                     xlab = "", ylab = "", bty = "n")
      
      graphics::axis(1, mgp = p_set$axis1_mgp, tck = -0.04, cex.axis = p_set$dyn_cex)
      graphics::axis(2, mgp = p_set$axis2_mgp, tck = -0.04, cex.axis = p_set$dyn_cex)
      
     
        graphics::mtext(sp, side = 1, cex = p_set$dyn_cex*0.8, line = 1.0, adj = 0) 
        graphics::mtext(datasource, side = 1, cex = p_set$dyn_cex*0.8, line = 1.8, adj = 0) 
        if (nrow(qdat) > 0) {
          graphics::mtext(sprintf("nll:%.3g  q:%.3g (%s)", tot, unique(sq)[1], unique(qdat$Type)[1]),
                          side = 1, line = 2.6, cex = p_set$dyn_cex*0.8, adj = 0)
        }
      
      sp_index <- which(Rsim.scenario$params$spname == sp)
      if (length(sp_index) > 0 && !is.na(Rsim.scenario$params$B_BaseRef[sp_index])) {
        graphics::abline(h = Rsim.scenario$params$B_BaseRef[sp_index], col = "darkred", lty = 3)
      }
      
      graphics::points(as.numeric(qdat$Year), mn)
      graphics::segments(as.numeric(qdat$Year), y0 = up, y1 = dn)
    }
  }
  
  graphics::mtext(paste("Fitting results for:", display_name, paste(", nll:", nll_val)),
                  outer = TRUE, side = 3, line = 1, cex = 1.5, font = 1)
  
  #catch plots
  nc_catch <- ceiling(sqrt(length(species.list)))
  nr_catch <- ceiling(length(species.list) / nc_catch)
  
  grDevices::dev.new(noRStudioGD=TRUE)
  graphics::par(mfrow = c(nr_catch, nc_catch), oma = c(0, 0, 3, 0))
  for(sp in species.list) {
    qdat <- catch.obj[catch.obj$Group == sp, ]
    graphics::par(mar= c(3,2,2,1))
    
    pin <- graphics::par("pin")
    dyn_cex <- max(0.75, min(0.95, min(pin) *0.2+0.6))
    
    mn   <- qdat$obs
    sdlog <- sqrt(log(1.0 + (qdat$sd * qdat$sd) / ifelse(mn == 0, 1e-6, mn * mn)))
    
    up <- mn * exp(1.96 * sdlog)
    dn <- mn / exp(1.96 * sdlog)
    est  <- run.out$annual_Catch[, sp]
    tot <- sum(qdat$fit * qdat$wt)
    
    all_vals <- c(up, est)
    valid_vals <- all_vals[!is.na(all_vals) & !is.infinite(all_vals)]
    ymax <- if(length(valid_vals) > 0) max(valid_vals) else 1
    if(ymax <= 0) ymax <- 1
    
    graphics::plot(
      as.numeric(rownames(run.out$annual_Catch)),
      est, type = "l", ylim = c(0, ymax), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n"
    )
    
    graphics::axis(1, mgp = c(3, 0.4, 0), tck = -0.04, cex.axis = dyn_cex) #FLAG####
    graphics::axis(2, mgp = c(3, 0.4, 0), tck = -0.04, cex.axis = dyn_cex)
    
    graphics::mtext(sp, side = 1, cex = dyn_cex*0.8, line = 1.0, adj = 0)
    graphics::mtext(sprintf("nll: %.3g", tot), side = 1, cex = dyn_cex*0.8, line = 1.8, adj = 0)
    
    graphics::points(as.numeric(qdat$Year), mn)
    graphics::segments(as.numeric(qdat$Year), y0 = up, y1 = dn)
    
    slist <- Rsim.scenario$params$spname[Rsim.scenario$params$FishFrom + 1]
    tcatch <- sum((Rsim.scenario$params$FishQ * Rsim.scenario$params$B_BaseRef[slist])[slist == sp], na.rm = TRUE)
    if (!is.na(tcatch) && tcatch > 0) {
      graphics::abline(h = tcatch, col = "darkred", lty = 3)
    }
  }
  
  graphics::mtext(paste("Catch fitting results for:", display_name, paste(", nll:", nll_val)),
                  outer = TRUE, side = 3, line = 1, cex = 1.5, font = 1)
  }

#################################################################################
#' Plot Full Biomass and Catch Fits for a Single Species
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
#' @inheritParams rsim.plot.fit 
#' @param species A single character string representing the species name.
#'
#' @return Invoked for its side effect of drawing a multi-panel plot to the active
#' graphics device. Returns `NULL`.
#'
#' @export
rsim.plot.fitspecies <- function(Rsim.scenario, Rsim.output, species) {
  
  if (length(species) > 1) {
    warning("rsim.plot.full expects a single species string. Using only the first element.")
    species <- species[1]
  }
  
  grDevices::dev.new(noRStudioGD = TRUE)
  
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(try(graphics::par(oldpar), silent=TRUE))
  
  # Biomass plotting
  fit_obj <- rsim.fit.obj(Rsim.scenario, Rsim.output)
  bio.obj <- fit_obj$Biomass
  sdat <- bio.obj[bio.obj$Group == species, ]
  sidset <- unique(sdat$Source)
  
  graphics::par(
    mfrow = c(1, max(1, length(sidset)) + 1),
    oma = c(0.0, 0.0, 0.0, 0.0),
    mar = c(4.5, 2.5, 2, 1)
  )
  
  pin <- graphics::par("pin")
  dyn_cex <- max(0.8, min(1.0, min(pin) * 0.2 + 0.6)) #FLAG####
  
  for (S in sidset) {
    qdat <- sdat[sdat$Source == S, ]
    sq <- qdat$survey_q
    sq[sq == 0 | is.na(sq)] <- 1e-6
    mn   <- qdat$obs / sq
    sdlog <- sqrt(log(1.0 + (qdat$sd / sq) * (qdat$sd / sq) / ifelse(mn == 0, 1e-6, mn * mn)))
    up <- mn * exp(1.96 * sdlog)
    dn <- mn / exp(1.96 * sdlog)
    
    est <- Rsim.output$annual_Biomass[, species]
    tot  <- sum(qdat$fit * qdat$wt)
    
    all_vals <- c(up, est)
    valid_vals <- all_vals[!is.na(all_vals) & !is.infinite(all_vals)]
    ymax <- if (length(valid_vals) > 0) max(valid_vals) else 1
    if (ymax <= 0) ymax <- 1
    
    graphics::plot(
      as.numeric(rownames(Rsim.output$annual_Biomass)),
      est, type = "l", ylim = c(0, ymax), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n"
    )
    
    graphics::axis(1, mgp = c(3, 1, 0), tck = -0.04, cex.axis = dyn_cex) #FLAG ####
    graphics::axis(2, mgp = c(3, 1, 0), tck = -0.04, cex.axis = dyn_cex) #FLAG ####
    
    graphics::mtext(paste(S, species, " B", sprintf(" nll: %.3g", tot)),
                    side = 1, cex= dyn_cex*0.85, line = 2.2, adj = 0)
    
    if (nrow(qdat) > 0) {
      graphics::mtext(sprintf("%s  q: %.3g", unique(qdat$Type)[1], unique(sq)[1]),
                      side = 1, line = 3.2, cex = dyn_cex*0.85, adj = 0) #FLAG ####
    }
    
    sp_index <- which(Rsim.scenario$params$spname == species)
    if (length(sp_index) > 0 && !is.na(Rsim.scenario$params$B_BaseRef[sp_index])) {
      graphics::abline(h = Rsim.scenario$params$B_BaseRef[sp_index], col = "darkred", lty = 3)
    }
    
    graphics::points(as.numeric(qdat$Year), mn)
    graphics::segments(as.numeric(qdat$Year), y0 = up, y1 = dn)
  }
  
  # Catch plotting
  catch.obj <- fit_obj$Catch
  qdat <- catch.obj[catch.obj$Group == species, ]
  
  if (nrow(qdat) == 0) {
    graphics::plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
    graphics::mtext(paste(species, "catch (No Data)"), side = 1, cex = dyn_cex*0.85, line = 2.2)
  } else {
    mn   <- qdat$obs
    sdlog <- sqrt(log(1.0 + (qdat$sd * qdat$sd) / ifelse(mn == 0, 1e-6, mn * mn)))
    up <- mn * exp(1.96 * sdlog)
    dn <- mn / exp(1.96 * sdlog)
    
    est  <- Rsim.output$annual_Catch[, species]
    tot <- sum(qdat$fit * qdat$wt)
    
    all_vals <- c(up, est)
    valid_vals <- all_vals[!is.na(all_vals) & !is.infinite(all_vals)]
    ymax <- if (length(valid_vals) > 0) max(valid_vals) else 1
    if (ymax <= 0) ymax <- 1
    
    graphics::plot(
      as.numeric(rownames(Rsim.output$annual_Catch)),
      est, type = "l", ylim = c(0, ymax), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n"
    )
    
    graphics::axis(1, mgp = c(3, 1, 0), tck = -0.04, cex.axis = dyn_cex)
    graphics::axis(2, mgp = c(3, 1, 0), tck = -0.04, cex.axis = dyn_cex)
    
    graphics::mtext(paste(species, "catch", sprintf("   nll: %.3g", tot)),
                    side = 1, cex= dyn_cex*0.85, line = 2.2, adj = 0)
    
    graphics::points(as.numeric(qdat$Year), mn)
    graphics::segments(as.numeric(qdat$Year), y0 = up, y1 = dn)
    
    slist <- Rsim.scenario$params$spname[Rsim.scenario$params$FishFrom + 1]
    tcatch <- sum((Rsim.scenario$params$FishQ * Rsim.scenario$params$B_BaseRef[slist])[slist == species], na.rm = TRUE)
    if (!is.na(tcatch) && tcatch > 0) {
      graphics::abline(h = tcatch, col = "darkred", lty = 3)
    }
  }
}

