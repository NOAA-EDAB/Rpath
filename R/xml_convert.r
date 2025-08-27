#' Creates an Rpath object from an EwE exported model (EIIXML format)
#'
#' @description
#' Reads in the eiixml model format (XML) created from exporting a model from Ecobase GUI.
#' Translates the XML to an rpath object which can then be balanced
#'
#' @param eiifile Full path to exported EwE XML file
#' @param verbose Logical. Use for debugging. If TRUE, prints out useful content
#' @param export Logical. Use for debugging. If TRUE, exports the data frames created in the function
#'
#' @return An Rpath object (list) with the following components:
#' \item{stanzas}{}
#' \item{pedigree}{}
#' \item{diet}{}
#' \item{model}{}
#'
#'
#' @export
create.rpath.from.eiixml <- function(
  eiifile,
  verbose = FALSE,
  export = FALSE
) {
  # Import the xml file and parse it into a list of data frames
  parsed_object <- import.eiixml(eiifile, verbose, export)

  # Extract and order group, gear and stanza names to create Rpath object---------

  # Order of groups on spreadsheets/in Rpath is listed in the
  # Sequence column.  This is different than the GroupID order
  # so a lookup is needed.  GroupID is the key used for From and To
  # flows so a lookup table is needed
  seq_look <- parsed_object$ewe_EcopathGroup$Sequence
  ord <- order(parsed_object$ewe_EcopathGroup$Sequence)
  ford <- order(parsed_object$ewe_EcopathFleet$Sequence)

  # Sort data tables by sequence order
  ordgroups <- parsed_object$ewe_EcopathGroup[ord, ]
  ordfleets <- parsed_object$ewe_EcopathFleet[ford, ]

  # Get and clean names and Types for living groups and gear
  bio_names <- janitor::make_clean_names(ordgroups$GroupName)
  gear_names <- janitor::make_clean_names(ordfleets$FleetName)

  if (sum(gear_names %in% bio_names) > 0) {
    warning(
      "There are Fleets with the same name as a bio group.
      These fleets will be renamed to avoid confusion.
      The string `_fleet` will be appended.\n
      Fleet names that match bio group names: ",
      paste0(gear_names[gear_names %in% bio_names], collapse = ",")
    )
  }
  gear_names <- ifelse(
    gear_names %in% bio_names,
    paste(gear_names, "fleet", sep = '_'),
    gear_names
  )
  det_names <- janitor::make_clean_names(ordgroups$GroupName[
    ordgroups$Type == 2
  ])
  live_names <- janitor::make_clean_names(ordgroups$GroupName[
    ordgroups$Type != 2
  ])

  g_names <- c(bio_names, gear_names)
  g_types <- c(ordgroups$Type, rep(3, length(gear_names)))
  row.names(ordgroups) <- bio_names

  # Also make a couple of named vectors associating names with groupID
  # (lookups for species and gear, using ID# as a character index)
  pnames <- bio_names
  names(pnames) <- ordgroups$GroupID
  gnames <- gear_names
  names(gnames) <- ordfleets$FleetID

  # Create the stanza table only if the model has stanzas defined in the xml

  if (!(identical(parsed_object$ewe_Stanza, NA))) {
    stanza_info <- make_stanza_table(
      parsed_object$ewe_Stanza,
      parsed_object$ewe_StanzaLifeStage,
      pnames,
      gnames
    )
    stanza_name_only <- stanza_info$stanza_name_only
    ordstanzas <- stanza_info$ordstanzas
    ordstages <- stanza_info$ordstages
  } else {
    # create a default stanza object
    stanza_name_only <- NULL
  }

  #-------------------------------------------------------------------------------
  # Create the Unbalanced Rpath object here
  unbal <- create.rpath.params(
    group = g_names,
    type = g_types,
    stgroup = stanza_name_only
  )

  # return(list(
  #   g_names = g_names,
  #   g_types = g_types,
  #   stanza_names_only = stanza_name_only
  # ))
  # ------------------------------------------------------------------------------
  # Populate unbal Rpath object with values

  # Fill group vectors
  gear_na <- rep(NA, length(gear_names))
  unbal$model$Biomass <- as.numeric(c(vec_na(ordgroups$Biomass), gear_na))
  unbal$model$PB <- as.numeric(c(vec_na(ordgroups$ProdBiom), gear_na))
  unbal$model$QB <- as.numeric(c(vec_na(ordgroups$ConsBiom), gear_na))
  unbal$model$EE <- as.numeric(c(vec_na(ordgroups$EcoEfficiency), gear_na))
  unbal$model$ProdCons <- as.numeric(c(vec_na(ordgroups$ProdCons), gear_na))
  unbal$model$BioAcc <- as.numeric(c(vec_na(ordgroups$BiomAcc), gear_na))
  unbal$model$Unassim <- as.numeric(c(vec_na(ordgroups$Unassim), gear_na))
  unbal$model$DetInput <- as.numeric(c(vec_na(ordgroups$DtImports), gear_na))

  #EwE output seems to have a few 0s or other numbers saved that should
  #be NAs, specifically with detritus.  Ensuring those don't sneak through.
  unbal$model$DetInput[ordgroups$Type != 2] <- NA
  unbal$model$EE[ordgroups$Type == 2] <- NA

  ewe_ordgroups <<- ordgroups

  # DIET TABLE---------------------------------------
  # TODO: where are diet imports in EwE XML data?
  diet_table <- parsed_object$ewe_EcopathDietComp
  diet_table$pred_name <- pnames[as.character(
    parsed_object$ewe_EcopathDietComp$PredID
  )]
  diet_table$prey_name <- pnames[as.character(
    parsed_object$ewe_EcopathDietComp$PreyID
  )]

  ppmat <- data.frame(unbal$diet[, -1])
  row.names(ppmat) <- unbal$diet$Group
  for (i in 1:length(diet_table$Diet)) {
    #cat(i,diet_table$prey_name[i], diet_table$pred_name[i],"\n")
    if (
      diet_table$prey_name[i] %in%
        rownames(ppmat) &
        diet_table$pred_name[i] %in% colnames(ppmat)
    ) {
      ppmat[
        diet_table$prey_name[i],
        diet_table$pred_name[i]
      ] <- diet_table$Diet[i]
    }
  }
  # need to convert matrices to data frames or data.table is unhappy
  unbal$diet[, 2:ncol(unbal$diet)] <- ppmat

  # add diet import
  impdiet <- ordgroups$ImpVar
  names(impdiet) <- rownames(ordgroups)
  predlist <- names(unbal$diet)[2:ncol(unbal$diet)]
  for (p in predlist) {
    unbal$diet[Group == "Import", p] <- impdiet[p]
  }

  # DETRITUS FATE FOR NON-GEAR -------------------------------------------------
  detframe <- data.frame(unbal$model)[, det_names]
  if (is.null(dim(detframe))) {
    detframe <- as.data.frame(detframe)
    colnames(detframe) <- det_names
  }
  row.names(detframe) <- g_names
  for (i in 1:length(diet_table$Diet)) {
    if (diet_table$prey_name[i] %in% det_names) {
      detframe[
        diet_table$pred_name[i],
        diet_table$prey_name[i]
      ] <- diet_table$DetritusFate[i]
    }
  }

  # eiixmls have some detritus fates NA/blank, replace with 0s
  detframe[is.na(detframe)] <- 0

  unbal$model[, det_names] <- detframe

  # CATCH AND DISCARDS --------------------------------------------------------
  if (any(is.na(parsed_object$ewe_EcopathCatch))) {
    catch_table <- data.frame(
      DiscardMortality = 1,
      Discards = 0.0,
      FleetID = 1,
      GroupID = seq(1, (length(pnames))),
      Landing = 0.0,
      Price = 0.0
    )
  } else {
    catch_table <- parsed_object$ewe_EcopathCatch
  }
  catch_table$gear_name <- gnames[as.character(catch_table$FleetID)]
  catch_table$group_name <- pnames[as.character(catch_table$GroupID)]

  discard_names <- paste(gear_names, "disc", sep = '.')

  fishframe <- data.frame(unbal$model)[, gear_names]
  if (is.null(dim(fishframe))) {
    fishframe <- as.data.frame(fishframe)
    colnames(fishframe) <- gear_names
  }

  row.names(fishframe) <- unbal$model$Group
  discardframe <- data.frame(unbal$model)[, discard_names]
  if (is.null(dim(discardframe))) {
    discardframe <- as.data.frame(discardframe)
    colnames(discardframe) <- discard_names
  }
  row.names(discardframe) <- unbal$model$Group

  for (i in 1:length(catch_table$Landing)) {
    if (
      catch_table$gear_name[i] %in%
        colnames(fishframe) &
        catch_table$group_name[i] %in% rownames(fishframe)
    ) {
      fishframe[
        catch_table$group_name[i],
        catch_table$gear_name[i]
      ] <- catch_table$Landing[i]
      discardframe[
        catch_table$group_name[i],
        paste(catch_table$gear_name[i], "disc", sep = '.')
      ] <- catch_table$Discards[i]
    }
  }

  unbal$model[, gear_names] <- fishframe
  unbal$model[, discard_names] <- discardframe

  # DETRITUS FATE FOR GEAR -----------------------------------------------------
  if (any(is.na(parsed_object$ewe_EcopathDiscardFate))) {
    # If no fleets/catch, route fleet1 detritus to last detritus group
    last_detritus <- as.numeric(names(which(
      pnames == det_names[length(det_names)]
    )))
    fate_table <- data.frame(
      DiscardFate = 1.0,
      FleetID = 1,
      GroupID = last_detritus
    )
  } else {
    fate_table <- parsed_object$ewe_EcopathDiscardFate
  }
  fate_table$gear_name <- gnames[as.character(fate_table$FleetID)]
  fate_table$group_name <- pnames[as.character(fate_table$GroupID)]

  fateframe <- data.frame(unbal$model)[unbal$model$Group %in% gnames, det_names]
  if (is.null(dim(fateframe))) {
    fateframe <- as.data.frame(fateframe)
    colnames(fateframe) <- det_names
  }

  row.names(fateframe) <- gnames
  for (i in 1:length(fate_table$DiscardFate)) {
    if (
      fate_table$gear_name[i] %in%
        rownames(fateframe) &
        fate_table$group_name[i] %in% colnames(fateframe)
    ) {
      fateframe[
        fate_table$gear_name[i],
        fate_table$group_name[i]
      ] <- fate_table$DiscardFate[i]
    }
  }

  unbal$model[
    (length(bio_names) + 1):(length(bio_names) + length(gear_names)),
    det_names
  ] <- fateframe

  # Configure the stanzas object
  if (!(identical(parsed_object$ewe_Stanza, NA))) {
    # STANZAS ----------------------------------------------------------------------

    unbal$stanzas$stgroups$Wmat <- ordstanzas[
      unbal$stanzas$stgroups$StanzaGroup,
      "WmatWinf"
    ]
    unbal$stanzas$stgroups$BAB <- ordstanzas[
      unbal$stanzas$stgroups$StanzaGroup,
      "BABsplit"
    ]
    unbal$stanzas$stgroups$RecPower <- ordstanzas[
      unbal$stanzas$stgroups$StanzaGroup,
      "RecPower"
    ]

    row.names(ordstages) <- ordstages$gname
    unbal$stanzas$stindiv$StanzaNum <- ordstages[
      unbal$stanzas$stindiv$Group,
      "Sequence"
    ]
    unbal$stanzas$stindiv$Leading <- ifelse(
      ordstages[unbal$stanzas$stindiv$Group, "LeadingLifeStage"] ==
        ordstages[unbal$stanzas$stindiv$Group, "Sequence"],
      TRUE,
      FALSE
    )
    unbal$stanzas$stindiv$First <- ordstages[
      unbal$stanzas$stindiv$Group,
      "AgeStart"
    ]

    # Not sure how Mort is stored outside of PB, or if it's different.
    pb <- unbal$model$PB
    names(pb) <- unbal$model$Group
    unbal$stanzas$stindiv$Z <- pb[unbal$stanzas$stindiv$Group]

    # Looping through in this weird way doesn't change the data.tables data type
    # (default logical) for the Last column.  (more foolishness)
    unbal$stanzas$stindiv[, Last := as.numeric(Last)]
    unbal$stanzas$stgroups[, VBGF_Ksp := as.numeric(VBGF_Ksp)]

    for (gg in 1:max(unbal$stanzas$stindiv$StGroupNum)) {
      # Finding vonBK is a mess involving all three tables (stored on group table)
      # Look up Leading Group and Stanza Group Number to get name on main table.
      unbal$stanzas$stgroups[StGroupNum == gg, "VBGF_Ksp"] <-
        ordgroups[
          as.character(unbal$stanzas$stindiv[
            StGroupNum == gg & Leading == T,
            "Group"
          ]),
        ]$vbK

      # Now loop through to add Last Month to first month
      szs <- unbal$stanzas$stindiv[StGroupNum == gg, ]
      for (ss in 1:(max(szs$StanzaNum) - 1)) {
        unbal$stanzas$stindiv[
          StGroupNum == gg & StanzaNum == ss,
          "Last"
        ] <- szs[szs$StanzaNum == ss + 1, "First"] - 1
      }
      ss <- ss + 1
      unbal$stanzas$stindiv[
        StGroupNum == gg & StanzaNum == ss,
        "Last"
      ] <- 999
    }

    # TODO RPATH - why does order matter here?
    unbal$stanzas$stindiv <- unbal$stanzas$stindiv[order(StGroupNum, StanzaNum)]
    #Not sure why Ksp is stored on the indiv table not the main group table
    #assuming the value for the leading stanzas is correct.
  }
  if (0) {
    cat("Created unbal (unbalanced ecopath) from ", eiifile, "\n")
  }
  return(unbal)
}

############################################################################
############################################################################
############################################################################
############################################################################
#' Reads in EwE exported XML file and parses into data frames
#'
#' @description
#' Parses the exported eiixml file (XML format, exported using Ecobase GUI) file into a list of data frames, one for each table in the XML.
#' To be useful in Rpath this needs to be further processed into an Rpath object.
#'
#' @inheritParams create.rpath.from.eiixml
#'
#' @return A list of data frames, one for each table in the XML file.
#'
#' @export
import.eiixml <- function(eiifile, verbose = F, export = F) {
  # Warn as you go, not at the end
  options(warn = 1)

  # create list object used to store output
  eweobject <- list()

  # Read and parse xml file
  dat <- xml2::read_xml(eiifile)

  # Create a list of Tables (main eiixml output data structure is "Table")
  tables <- xml2::xml_find_all(dat, ".//Table")

  # Loop through the nodes (Tables) and read into data frames
  for (node in tables) {
    # how to select a named node from the list of tables:
    node_name <- xml2::xml_attr(node, "Name")
    #cat(node_name,""); #flush.console()
    # Parse the Columns attribute of each Table which is var name/var type.
    # This messy command splits column names by , and by : then discarding
    # the : which is variable type.  I'm not really sure how the `[[` works but
    # it's choosing the 1st element from each list.
    cols <- unlist(lapply(
      strsplit(unlist(strsplit(xml2::xml_attr(node, "Columns"), ",")), ":"),
      `[[`,
      1
    ))
    # Get rows of data and split into a data frame with that data
    rawdat <- xml2::xml_text(xml2::xml_find_all(node, "Row"))

    # special reformatting of some nodes node to avoid comma issues
    # We can't split by commas, like we can for other nodes, because
    # there are some fields with commas in the text (within quotes),
    # so we need to grep these lines,
    # and replace the commas with whitespace for all occurrences

    if (
      node_name %in%
        c(
          "EcopathModel",
          "Pedigree",
          "EcopathGroup",
          "EcopathFleet",
          "Auxillary",
          "UpdateLog"
        ) &&
        (length(rawdat) > 0)
    ) {
      input_string <- rawdat
      for (irow in 1:length(input_string)) {
        row_string <- input_string[irow]
        match_list <- gregexpr("\"[^\"]*\"", row_string)
        quoted_content <- regmatches(row_string, match_list)
        modified_quoted_content <- lapply(quoted_content, function(x) {
          gsub(",", " ", x, fixed = TRUE)
        })
        modified_string <- row_string
        regmatches(modified_string, match_list) <- modified_quoted_content
        row_string <- paste0(modified_string, " ")
        rawdat[irow] <- row_string
      }
    }

    if (length(rawdat) > 0) {
      if (verbose) {
        cat(node_name, "")
      }
      #dtable <- data.frame(matrix(as.numeric(unlist(strsplit(rawdat,","))), ncol=length(cols), byrow=T))
      unlisted_dat <- (unlist(strsplit(rawdat, ",")))
      if (length(unlisted_dat) %% length(cols) != 0) {
        warning("Matrix alignment (comma issue?) in table ", node_name, ":")
      }
      dtable <- type.convert(
        data.frame(matrix(unlisted_dat, ncol = length(cols), byrow = T)),
        as.is = T
      )
      if (verbose) {
        cat(length(rawdat), "rows\n")
      }
      names(dtable) <- cols
      eweobject[[paste0("ewe_", node_name)]] <- dtable
      #assign(paste0("eweobject$ewe_", node_name), dtable)
      if (export) {
        do.call("<<-", list(paste0("ewe_", node_name), dtable))
      }
    } else {
      dtable <- data.frame(matrix(ncol = length(cols), byrow = T))[-1, ]
      # Assign default value of NA to all variables that dont exist for model
      eweobject[[paste0("ewe_", node_name)]] <- NA
      #assign(paste0("eweobject$ewe_", node_name), NULL)
      #cat("0\n")
    }
    # Use these here if you want to make empty tables where there's no data
    #names(dtable) <- cols
    #assign(paste("ewe_",node_name,sep="") , dtable)
  } # end node loop

  ## Create dummy values for missing components
  if (identical(eweobject$ewe_EcopathFleet, NA)) {
    warning("No Fleets are present. A dummy fleet will be created")
    eweobject$ewe_EcopathFleet <- data.frame(
      FixedCost = 0,
      FleetID = 1,
      FleetName = "Fleet1",
      NominalEffort = NA,
      PoolColor = 0,
      SailingCost = 1,
      Sequence = 1,
      VariableCost = 0
    )
  }

  return(eweobject)

  ## ------ ALL RAW DATA from XML should be in R data frames at this point -------
}

############################################################################
############################################################################
############################################################################
############################################################################
#' Function to replace negative values with NAs
#'
#' Utility function
#' @noRd
vec_na <- function(vec) {
  return(ifelse(vec < (-0.1), NA, vec))
}
############################################################################
############################################################################
############################################################################
############################################################################
#' Create the Stanza table
#'
#' @description
#' Creates a table of stanzas and their life stages from the EwE Ecobase eiixml
#'
#' @param ewe_Stanza Stanza data from XML node
#' @param ewe_StanzaLifeStage Stanza life stage data from XML node
#' @param pnames Names of biological groups in the model
#' @param gnames Names of gear groups in the model
#'
#' @return list
#' \item{ordstanzas}{}
#' \item{stanza_name_only}{}
#' \item{ordstages}{}
#'
#' @noRd
make_stanza_table <- function(ewe_Stanza, ewe_StanzaLifeStage, pnames, gnames) {
  # Get and clean stanza names, order stanza table
  sord <- order(ewe_Stanza$StanzaID)
  ordstanzas <- ewe_Stanza[sord, ]
  ordstanzas$stanza_names <- janitor::make_clean_names(ordstanzas$StanzaName)
  row.names(ordstanzas) <- as.character(ordstanzas$StanzaID)

  # Combine the life stages and stanza info into one table (making some duplication
  # as stanza parameters will have multiple copies).
  ordstages <- cbind(
    ewe_StanzaLifeStage,
    ordstanzas[as.character(ewe_StanzaLifeStage$StanzaID), ]
  )

  ordstages$gname <- pnames[as.character(ordstages$GroupID)]

  stanza_column <- rep(NA, length(pnames))
  names(stanza_column) <- names(pnames)
  stanza_column[as.character(ordstages$GroupID)] <- ordstages$stanza_names
  stanza_name_only <- c(as.character(stanza_column), rep(NA, length(gnames)))

  # Rename ordered tables to use cleaned names, not GroupIDs
  row.names(ordstanzas) <- ordstanzas$stanza_names

  return(list(
    ordstanzas = ordstanzas,
    stanza_name_only = stanza_name_only,
    ordstages = ordstages
  ))
}
