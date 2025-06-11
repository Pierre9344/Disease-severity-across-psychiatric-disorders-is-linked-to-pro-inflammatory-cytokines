GetAndFormatClinicalData <- function(NonImputedFile, ImputedFile) {
  load(NonImputedFile)
  load(ImputedFile)
  # add the visits (the row are in order)
  pat_nonimp_proposal %<>%
    group_by(v1_id) %>%
    mutate(visit = row_number()) %>%
    ungroup()

  pat_imp_proposal$v1_scid_dsm_dx_cat <- NA
  pat_imp_proposal$gsa_id <- NA
  pat_imp_proposal$smRNAome_id <- NA
  pat_imp_proposal$sex <- NA
  pat_nonimp_proposal$v1_scid_dsm_dx_cat <- as.character(pat_nonimp_proposal$v1_scid_dsm_dx_cat)
  for (i in 1:nrow(pat_imp_proposal)) {
    for (j in 1:nrow(pat_nonimp_proposal)) {
      if (pat_imp_proposal$v1_id[i] == pat_nonimp_proposal$v1_id[j]) {
        pat_imp_proposal$v1_scid_dsm_dx_cat[i] <- pat_nonimp_proposal$v1_scid_dsm_dx_cat[j]
        pat_imp_proposal$gsa_id[i] <- pat_nonimp_proposal$gsa_id[j]
        pat_imp_proposal$smRNAome_id[i] <- pat_nonimp_proposal$smRNAome_id[j]
        pat_imp_proposal$sex[i] <- as.character(pat_nonimp_proposal$v1_sex[j])
        break
      }
    }
  }
  # Correct the age that were not correctly recorded
  # The loop next was run outside the pipeline to detect the ages that were not recorded correctly

  # for (i in unique(pat_imp_proposal$v1_id)) {
  #   v1 <- pat_imp_proposal %>% filter(v1_id == i, visit == 1)
  #   v2 <- pat_imp_proposal %>% filter(v1_id == i, visit == 2)
  #   v3 <- pat_imp_proposal %>% filter(v1_id == i, visit == 3)
  #   v4 <- pat_imp_proposal %>% filter(v1_id == i, visit == 4)
  #   if (v1$age > v2$age) {
  #     print(i)
  #   } else if (v2$age > v3$age) {
  #     print(i)
  #   } else if (v3$age > v4$age) {
  #     print(i)
  #   }
  # }

  # The loop next was used to correct the age
  # informations was computed by the IPPG using the visit date and birthday dates of the participants

  for (i in 1:nrow(pat_imp_proposal)) {
    if (pat_imp_proposal$v1_id[i] == "digt726" && pat_imp_proposal$visit[i] == 3) {
      pat_imp_proposal$age[i] <- 77
    } else if (pat_imp_proposal$v1_id[i] == "gzzb222" && pat_imp_proposal$visit[i] == 2) {
      pat_imp_proposal$age[i] <- 64
    } else if (pat_imp_proposal$v1_id[i] == "iyxn298" && pat_imp_proposal$visit[i] == 3) {
      pat_imp_proposal$age[i] <- 36
    } else if (pat_imp_proposal$v1_id[i] == "pjug754" && pat_imp_proposal$visit[i] == 2) {
      pat_imp_proposal$age[i] <- 37
    } else if (pat_imp_proposal$v1_id[i] == "qmdk332" && pat_imp_proposal$visit[i] == 3) {
      pat_imp_proposal$age[i] <- 74
    } else if (pat_imp_proposal$v1_id[i] == "uaao152" && pat_imp_proposal$visit[i] == 4) {
      pat_imp_proposal$age[i] <- 62
    } else if (pat_imp_proposal$v1_id[i] == "xpfu752" && pat_imp_proposal$visit[i] == 3) {
      pat_imp_proposal$age[i] <- 45
    } else if (pat_imp_proposal$v1_id[i] == "oegm588" && pat_imp_proposal$visit[i] == 2) {
      pat_imp_proposal$age[i] <- 50
    } else if (pat_imp_proposal$v1_id[i] == "oegm588" && pat_imp_proposal$visit[i] == 3) {
      pat_imp_proposal$age[i] <- 50
    }
  }

  # An individual was recorded with an aberrant imputed gaf value inferior to 0 while it should be between 0 and 100
  # We removed this individual from the analysis as the gaf will be used by other function
  pat_imp_proposal %<>% filter(gaf >= 0)

  # we keep only the individual with 4 visits and patients with BD, SCZ, TSCZA or MDD
  # (other classifications are too different clinicaly and represent a really small proportion)
  toKeep <- as.data.frame(table(pat_imp_proposal$v1_id)) %>%
    filter(Freq == 4)
  pat_imp_proposal %<>%
    filter(v1_id %in% toKeep$Var1) %>%
    filter(v1_scid_dsm_dx_cat %in% c(
      "Bipolar-I Disorder", "Bipolar-II Disorder",
      "Depression", "ICD-10 Schizophrenia",
      "Schizoaffective Disorder", "Schizophrenia"
    )) %>%
    droplevels()
  # Individual annotated as ICD-10 Schizophrenia were diagnosed according to ICD-10 criteria.
  # We change it to "Schizophrenia" to simplify the analysis
  for (i in 1:nrow(pat_imp_proposal)) {
    if (pat_imp_proposal$v1_scid_dsm_dx_cat[i] == "ICD-10 Schizophrenia") {
      pat_imp_proposal$v1_scid_dsm_dx_cat[i] <- "Schizophrenia"
    }
  }

  # Datatype conversion
  # We convert the data to get the same format as in the PsyCourse codebook

  pat_imp_proposal$curr_paid_empl <- factor(pat_imp_proposal$curr_paid_empl, levels = c("N", "Y"))
  pat_imp_proposal$partner <- factor(pat_imp_proposal$partner, levels = c("N", "Y"))
  pat_imp_proposal$med_pst_wk <- factor(pat_imp_proposal$med_pst_wk, levels = c(1, 2, 3, 4, 5, 6))
  pat_imp_proposal$sex <- factor(pat_imp_proposal$sex, levels = c("F", "M"))
  pat_imp_proposal$idsc_itm1 <- factor(pat_imp_proposal$idsc_itm1, levels = c(0, 1, 2, 3))
  pat_imp_proposal$idsc_itm10 <- factor(pat_imp_proposal$idsc_itm10, levels = c(0, 1, 2, 3))
  pat_imp_proposal$idsc_itm15 <- factor(pat_imp_proposal$idsc_itm15, levels = c(0, 1, 2, 3))
  pat_imp_proposal$idsc_itm16 <- factor(pat_imp_proposal$idsc_itm16, levels = c(0, 1, 2, 3))
  pat_imp_proposal$idsc_itm17 <- factor(pat_imp_proposal$idsc_itm17, levels = c(0, 1, 2, 3))
  pat_imp_proposal$idsc_itm18 <- factor(pat_imp_proposal$idsc_itm18, levels = c(0, 1, 2, 3))
  pat_imp_proposal$idsc_itm19 <- factor(pat_imp_proposal$idsc_itm19, levels = c(0, 1, 2, 3))
  pat_imp_proposal$idsc_itm2 <- factor(pat_imp_proposal$idsc_itm2, levels = c(0, 1, 2, 3))
  pat_imp_proposal$idsc_itm20 <- factor(pat_imp_proposal$idsc_itm20, levels = c(0, 1, 2, 3))
  pat_imp_proposal$idsc_itm21 <- factor(pat_imp_proposal$idsc_itm21, levels = c(0, 1, 2, 3))
  pat_imp_proposal$idsc_itm22 <- factor(pat_imp_proposal$idsc_itm22, levels = c(0, 1, 2, 3))
  pat_imp_proposal$idsc_itm23 <- factor(pat_imp_proposal$idsc_itm23, levels = c(0, 1, 2, 3))
  pat_imp_proposal$idsc_itm24 <- factor(pat_imp_proposal$idsc_itm24, levels = c(0, 1, 2, 3))
  pat_imp_proposal$idsc_itm25 <- factor(pat_imp_proposal$idsc_itm25, levels = c(0, 1, 2, 3))
  pat_imp_proposal$idsc_itm27 <- factor(pat_imp_proposal$idsc_itm27, levels = c(0, 1, 2, 3))
  pat_imp_proposal$idsc_itm28 <- factor(pat_imp_proposal$idsc_itm28, levels = c(0, 1, 2, 3))
  pat_imp_proposal$idsc_itm29 <- factor(pat_imp_proposal$idsc_itm29, levels = c(0, 1, 2, 3))
  pat_imp_proposal$idsc_itm30 <- factor(pat_imp_proposal$idsc_itm30, levels = c(0, 1, 2, 3))
  pat_imp_proposal$idsc_itm3 <- factor(pat_imp_proposal$idsc_itm3, levels = c(0, 1, 2, 3))
  pat_imp_proposal$idsc_itm4 <- factor(pat_imp_proposal$idsc_itm4, levels = c(0, 1, 2, 3))
  pat_imp_proposal$idsc_itm5 <- factor(pat_imp_proposal$idsc_itm5, levels = c(0, 1, 2, 3))
  pat_imp_proposal$idsc_itm6 <- factor(pat_imp_proposal$idsc_itm6, levels = c(0, 1, 2, 3))
  pat_imp_proposal$idsc_itm7 <- factor(pat_imp_proposal$idsc_itm7, levels = c(0, 1, 2, 3))
  pat_imp_proposal$idsc_itm8 <- factor(pat_imp_proposal$idsc_itm8, levels = c(0, 1, 2, 3))
  pat_imp_proposal$idsc_itm9 <- factor(pat_imp_proposal$idsc_itm9, levels = c(0, 1, 2, 3))
  pat_imp_proposal$panss_p1 <- factor(pat_imp_proposal$panss_p1, levels = c(1, 2, 3, 4, 5, 6, 7))
  pat_imp_proposal$panss_p2 <- factor(pat_imp_proposal$panss_p2, levels = c(1, 2, 3, 4, 5, 6, 7))
  pat_imp_proposal$panss_p3 <- factor(pat_imp_proposal$panss_p3, levels = c(1, 2, 3, 4, 5, 6, 7))
  pat_imp_proposal$panss_p4 <- factor(pat_imp_proposal$panss_p4, levels = c(1, 2, 3, 4, 5, 6, 7))
  pat_imp_proposal$panss_p5 <- factor(pat_imp_proposal$panss_p5, levels = c(1, 2, 3, 4, 5, 6, 7))
  pat_imp_proposal$panss_p6 <- factor(pat_imp_proposal$panss_p6, levels = c(1, 2, 3, 4, 5, 6, 7))
  pat_imp_proposal$panss_p7 <- factor(pat_imp_proposal$panss_p7, levels = c(1, 2, 3, 4, 5, 6, 7))
  pat_imp_proposal$panss_n1 <- factor(pat_imp_proposal$panss_n1, levels = c(1, 2, 3, 4, 5, 6, 7))
  pat_imp_proposal$panss_n2 <- factor(pat_imp_proposal$panss_n2, levels = c(1, 2, 3, 4, 5, 6, 7))
  pat_imp_proposal$panss_n3 <- factor(pat_imp_proposal$panss_n3, levels = c(1, 2, 3, 4, 5, 6, 7))
  pat_imp_proposal$panss_n4 <- factor(pat_imp_proposal$panss_n4, levels = c(1, 2, 3, 4, 5, 6, 7))
  pat_imp_proposal$panss_n5 <- factor(pat_imp_proposal$panss_n5, levels = c(1, 2, 3, 4, 5, 6, 7))
  pat_imp_proposal$panss_n6 <- factor(pat_imp_proposal$panss_n6, levels = c(1, 2, 3, 4, 5, 6, 7))
  pat_imp_proposal$panss_n7 <- factor(pat_imp_proposal$panss_n7, levels = c(1, 2, 3, 4, 5, 6, 7))
  pat_imp_proposal$panss_g1 <- factor(pat_imp_proposal$panss_g1, levels = c(1, 2, 3, 4, 5, 6, 7))
  pat_imp_proposal$panss_g2 <- factor(pat_imp_proposal$panss_g2, levels = c(1, 2, 3, 4, 5, 6, 7))
  pat_imp_proposal$panss_g3 <- factor(pat_imp_proposal$panss_g3, levels = c(1, 2, 3, 4, 5, 6, 7))
  pat_imp_proposal$panss_g4 <- factor(pat_imp_proposal$panss_g4, levels = c(1, 2, 3, 4, 5, 6, 7))
  pat_imp_proposal$panss_g5 <- factor(pat_imp_proposal$panss_g5, levels = c(1, 2, 3, 4, 5, 6, 7))
  pat_imp_proposal$panss_g6 <- factor(pat_imp_proposal$panss_g6, levels = c(1, 2, 3, 4, 5, 6, 7))
  pat_imp_proposal$panss_g7 <- factor(pat_imp_proposal$panss_g7, levels = c(1, 2, 3, 4, 5, 6, 7))
  pat_imp_proposal$panss_g8 <- factor(pat_imp_proposal$panss_g9, levels = c(1, 2, 3, 4, 5, 6, 7))
  pat_imp_proposal$panss_g9 <- factor(pat_imp_proposal$panss_g9, levels = c(1, 2, 3, 4, 5, 6, 7))
  pat_imp_proposal$panss_g10 <- factor(pat_imp_proposal$panss_g10, levels = c(1, 2, 3, 4, 5, 6, 7))
  pat_imp_proposal$panss_g11 <- factor(pat_imp_proposal$panss_g11, levels = c(1, 2, 3, 4, 5, 6, 7))
  pat_imp_proposal$panss_g12 <- factor(pat_imp_proposal$panss_g12, levels = c(1, 2, 3, 4, 5, 6, 7))
  pat_imp_proposal$panss_g13 <- factor(pat_imp_proposal$panss_g13, levels = c(1, 2, 3, 4, 5, 6, 7))
  pat_imp_proposal$panss_g14 <- factor(pat_imp_proposal$panss_g14, levels = c(1, 2, 3, 4, 5, 6, 7))
  pat_imp_proposal$panss_g15 <- factor(pat_imp_proposal$panss_g15, levels = c(1, 2, 3, 4, 5, 6, 7))
  pat_imp_proposal$panss_g16 <- factor(pat_imp_proposal$panss_g16, levels = c(1, 2, 3, 4, 5, 6, 7))
  pat_imp_proposal$ymrs_itm1 <- factor(pat_imp_proposal$ymrs_itm1, levels = c(0, 1, 2, 3, 4))
  pat_imp_proposal$ymrs_itm2 <- factor(pat_imp_proposal$ymrs_itm2, levels = c(0, 1, 2, 3, 4))
  pat_imp_proposal$ymrs_itm3 <- factor(pat_imp_proposal$ymrs_itm3, levels = c(0, 1, 2, 3, 4))
  pat_imp_proposal$ymrs_itm4 <- factor(pat_imp_proposal$ymrs_itm4, levels = c(0, 1, 2, 3, 4))
  pat_imp_proposal$ymrs_itm5 <- factor(pat_imp_proposal$ymrs_itm5, levels = c(0, 2, 4, 6, 8))
  pat_imp_proposal$ymrs_itm6 <- factor(pat_imp_proposal$ymrs_itm6, levels = c(0, 2, 4, 6, 8))
  pat_imp_proposal$ymrs_itm7 <- factor(pat_imp_proposal$ymrs_itm7, levels = c(0, 1, 2, 3, 4))
  pat_imp_proposal$ymrs_itm8 <- factor(pat_imp_proposal$ymrs_itm8, levels = c(0, 2, 4, 6, 8))
  pat_imp_proposal$ymrs_itm9 <- factor(pat_imp_proposal$ymrs_itm9, levels = c(0, 2, 4, 6, 8))
  pat_imp_proposal$ymrs_itm10 <- factor(pat_imp_proposal$ymrs_itm10, levels = c(0, 1, 2, 3, 4))
  pat_imp_proposal$ymrs_itm11 <- factor(pat_imp_proposal$ymrs_itm11, levels = c(0, 1, 2, 3, 4))
  pat_imp_proposal$v1_scid_dsm_dx_cat <- factor(pat_imp_proposal$v1_scid_dsm_dx_cat)

  # The lines bellow correspond to covariates with missing values that we can't use for our analysis
  pat_imp_proposal$idsc_itm26 <- NULL
  pat_imp_proposal$idsc_itm11 <- NULL
  pat_imp_proposal$idsc_itm12 <- NULL
  pat_imp_proposal$idsc_itm13 <- NULL
  pat_imp_proposal$idsc_itm14 <- NULL
  pat_imp_proposal$nrpsy_lng <- NULL
  pat_imp_proposal$nrpsy_mtv <- NULL
  pat_imp_proposal$nrpsy_tmt_A_err <- NULL
  pat_imp_proposal$nrpsy_tmt_B_err <- NULL

  return(pat_imp_proposal)
}

GetClassificationDimensions <- function(ClinicalData) {
  # create a list with, for each dimensions, the name or the start of the names of the items recorded during the PsyCourse Study.
  ls <- list()
  ls[["DIM"]][["PANSS"]][["Names"]] <- c("panss_p", "panss_n", "panss_g")
  ls[["DIM"]][["YMRS"]][["Names"]] <- c("ymrs_itm")
  ls[["DIM"]][["IDSC"]][["Names"]] <- c("idsc_itm")
  ls[["DIM"]][["COGNITION"]][["Names"]] <- c("nrpsy_")
  ls[["DIM"]][["FUNCTIONNING"]][["Names"]] <- c(
    "partner", "Antidepressants", "Antipsychotics",
    "Mood_stabilizers", "Tranquilizers", "gaf"
  )

  # compute the different possible combinations of those variables we can use with longmixr
  for (i in 1:5) {
    ls[["VarComb"]] %<>% append(utils::combn(names(ls[["DIM"]]), i, simplify = F) %>%
      sapply(function(x) paste(x, collapse = "_"))) %>% unique()
  }

  # Reduce the dimensions
  for (i in names(ls[["DIM"]])) {
    if (i == "COGNITION") {
      ls[["DIM"]][[i]][["dim"]] <- ClinicalData %>%
        select(tidyr::starts_with(ls[["DIM"]][[i]][["Names"]])) %>%
        stats::prcomp(scale = T)
      ls[["DIM"]][[i]][["comp"]] <- as.data.frame(ls[["DIM"]][[i]][["dim"]]$x[, 1])
    } else {
      ls[["DIM"]][[i]][["dim"]] <- ClinicalData %>%
        select(tidyr::starts_with(ls[["DIM"]][[i]][["Names"]])) %>%
        FactoMineR::FAMD(ncp = 10, graph = F)
      ls[["DIM"]][[i]][["comp"]] <- as.data.frame(ls[["DIM"]][[i]][["dim"]]$ind$coord[, 1])
    }
    colnames(ls[["DIM"]][[i]][["comp"]]) <- i
  }

  # create a dataframe with the dimensions to give to longmixr
  ls[["DIM_CLUSTER"]] <- bind_cols(
    data.frame(
      ID = ClinicalData$v1_id,
      visit = ClinicalData$visit,
      age = ClinicalData$age,
      sex = as.numeric(ClinicalData$sex)
    ),
    PANSS = ls[["DIM"]][["PANSS"]][["comp"]],
    IDSC = ls[["DIM"]][["IDSC"]][["comp"]],
    YMRS = ls[["DIM"]][["YMRS"]][["comp"]],
    FUNCTIONNING = ls[["DIM"]][["FUNCTIONNING"]][["comp"]],
    COGNITION = ls[["DIM"]][["COGNITION"]][["comp"]]
  )
  return(ls)
}


GenerateResiduals <- function(x, age, sex, ID, correctAge = T, correctSEX = F) {
  data <- data.frame(x = x, age = age, ID = ID)
  resid <- NA
  if (correctAge & correctSEX) {
    model <- lme4::lmer(x ~ age * sex + (1 | ID), data = data)
    resid <- BiocGenerics::residuals(model, type = "response")
    names(resid) <- NULL
  } else if (correctAge) {
    model <- lme4::lmer(x ~ age + (1 | ID), data = data)
    resid <- BiocGenerics::residuals(model, type = "response")
    names(resid) <- NULL
  } else if (correctSEX) {
    model <- lme4::lmer(x ~ sex + (1 | ID), data = data)
    resid <- BiocGenerics::residuals(model, type = "response")
    names(resid) <- NULL
  } else {
    resid <- x
  }
  resid
}


DimensionsCorrection <- function(dim_data) {
  dim_data <- dim_data[["DIM_CLUSTER"]]
  dim_data$sex <- factor(ifelse(dim_data$sex == 1, "F", "M"), levels = c("F", "M"))

  dim_data_resid <- dim_data %>%
    mutate(across(all_of("COGNITION"),
      ~ GenerateResiduals(
        x = .x, age = age, sex = sex,
        ID = ID, correctAge = F, correctSEX = T
      ),
      .names = "{.col}_resid"
    )) %>%
    mutate(across(all_of("PANSS"),
      ~ GenerateResiduals(
        x = .x, age = age, sex = sex,
        ID = ID, correctAge = F, correctSEX = F
      ),
      .names = "{.col}_resid"
    )) %>%
    mutate(across(all_of("IDSC"),
      ~ GenerateResiduals(
        x = .x, age = age, sex = sex,
        ID = ID, correctAge = F, correctSEX = T
      ),
      .names = "{.col}_resid"
    )) %>%
    mutate(across(all_of("YMRS"),
      ~ GenerateResiduals(
        x = .x, age = age, sex = sex,
        ID = ID, correctAge = F, correctSEX = F
      ),
      .names = "{.col}_resid"
    )) %>%
    mutate(across(all_of("FUNCTIONNING"),
      ~ GenerateResiduals(
        x = .x, age = age, sex = sex,
        ID = ID, correctAge = T, correctSEX = F
      ),
      .names = "{.col}_resid"
    )) %>%
    select(ID, visit, ends_with("resid"), age, sex) %>%
    mutate(sex = as.numeric(sex))
  return(dim_data_resid)
}


LongmixrClustering <- function(dim_corrected, dim_to_use, seed, permut, max_cluster) {
  cat("\n\nStarting longmixr: ", dim_to_use, "\n\n", sep = "")
  # convert the sex variable as a factor for the longmixr algorithm
  dim_corrected$sex <- factor(ifelse(dim_corrected$sex == 1, "F", "M"), levels = c("F", "M"))

  list_models <- lapply(dim_to_use %>% strsplit("_") %>% unlist() %>% paste0("_resid"), function(x) {
    flexmix::FLXMRmgcv(as.formula(paste0(x, " ~ .")))
  })

  ls <- list()
  ls[["cluster_model"]] <- longitudinal_consensus_cluster(
    data = dim_corrected,
    id_column = "ID", max_k = max_cluster,
    reps = permut, p_item = 0.8,
    seed = seed,
    model_list = list_models,
    flexmix_formula = as.formula("~s(visit, k = 4) | ID"),
    final_linkage = "ward.D2", verbose = T
  )
  ls[["dimensions"]] <- dim_to_use
  ls[["permut"]] <- permut
  ls[["seed"]] <- seed
  list(result = ls)
}


GetClusterAssignments <- function(clusters, VarClustering, ImputedVar, NonImputedVar,
                                  Centers, PRS_MD, PRS_BD, PRS_SCZ, PRS_TREATMENT,
                                  PRS_PFACTOR, BMI) {
  assignments <- NA
  FirstLoop <- TRUE
  for (i in names(clusters)) {
    if (FirstLoop) {
      assignments <- longmixr::get_clusters(clusters[[i]][["cluster_model"]]) %>%
        mutate(across(
          starts_with("assignment_num_clus"),
          ~ recode(.x, `1` = "A", `2` = "B", `3` = "C", `4` = "D", `5` = "E")
        ))
      colnames(assignments) %<>% stringr::str_replace_all(
        "assignment_num_clus",
        clusters[[i]][["dimensions"]]
      )
      FirstLoop <- FALSE
    } else {
      assignments2 <- longmixr::get_clusters(clusters[[i]][["cluster_model"]]) %>%
        mutate(across(
          starts_with("assignment_num_clus"),
          ~ recode(.x, `1` = "A", `2` = "B", `3` = "C", `4` = "D", `5` = "E")
        ))
      colnames(assignments2) %<>% stringr::str_replace_all(
        "assignment_num_clus",
        clusters[[i]][["dimensions"]]
      )
      assignments %<>% base::merge(assignments2, by = "ID")
    }
  }


  df <- merge(assignments, VarClustering, by.x = "ID", by.y = "v1_id") %>%
    select("ID", visit, all_of(colnames(assignments)), sex)
  load(NonImputedVar)
  load(ImputedVar)
  load(Centers)
  pat_imp_proposal %<>% filter(v1_id %in% df$ID)
  pat_nonimp_proposal %<>% filter(v1_id %in% df$ID)
  df <- merge(df, pat_imp_proposal,
    by.x = c("ID", "visit"),
    by.y = c("v1_id", "visit"),
    all.y = F
  )

  # add the non imputed value that were not used for the clustering
  pat_nonimp_proposal2 <- pat_nonimp_proposal %>%
    select(v1_id | !any_of(colnames(df))) %>%
    select(-v1_evr_ill_drg.1)
  for (i in colnames(pat_nonimp_proposal2)) {
    if (!(i %in% colnames(df))) {
      # we only want the variables that are not already in our dataframe
      print(i)
      df[i] <- NA
      for (j in 1:nrow(df)) {
        for (k in 1:nrow(pat_nonimp_proposal2)) {
          # Check if the IDs are the same. No need to check for the visit are those value are only measured once
          if (pat_nonimp_proposal2$v1_id[k] == df$ID[j]) {
            df[j, i] <- ifelse(is.factor(pat_nonimp_proposal2[, i]),
              as.character(pat_nonimp_proposal2[k, i]),
              pat_nonimp_proposal2[k, i]
            )
            break
          }
        }
      }
    }
  }
  # remove duplicated column
  df$v1_id.1 <- NULL
  cent$v1_stat <- NULL
  cent$v1_sex <- NULL
  # add the clinical centers informations
  for (i in colnames(cent)) {
    if (!(i %in% c("v1_id", "v1_sex"))) {
      # we already have the IDs and sex of the individuals
      if (!(i %in% colnames(df))) {
        print(i)
        df[i] <- NA
        for (j in 1:nrow(df)) {
          for (k in 1:nrow(cent)) {
            if (cent$v1_id[k] == df$ID[j]) {
              df[j, i] <- cent[k, i]
              break
            }
          }
        }
      }
    }
  }

  # Add polygenic risk scores and reformat some columns
  PRS_mdd <- read.delim(PRS_MD)
  PRS_mdd <- subset(PRS_mdd, PRS_mdd$gsa_id %in% pat_nonimp_proposal$gsa_id)
  PRS_bd <- read.delim(PRS_BD)
  PRS_bd <- subset(PRS_bd, PRS_bd$gsa_id %in% pat_nonimp_proposal$gsa_id)
  PRS_scz <- read.delim(PRS_SCZ)
  PRS_scz <- subset(PRS_scz, PRS_scz$gsa_id %in% pat_nonimp_proposal$gsa_id)
  # used phi = 1e.02 to get a broad view of the associated variants as we want to compare the PRS across the clusters that will be defined by longmix
  PRS <- data.frame(
    PRS_bd = PRS_bd$PsyCourse_MimicSS_bip3_noPsyCourse_phi_1e.02_score,
    PRS_SCZ = PRS_scz$PsyCourse_MimicSS_scz3_phi_1e.02_score,
    PRS_MD = PRS_mdd$PsyCourse_MimicSS_mdd2018_phi_1e.02_score
  )
  PRS$PRS_sum <- rowSums(PRS)
  PRS$gsa_id <- PRS_bd$gsa_id

  df <- df %>%
    left_join(PRS, by = "gsa_id") %>%
    mutate(
      nrpsy_lng = as.factor(nrpsy_lng),
      nrpsy_mtv = as.factor(nrpsy_mtv),
      v4_opcrit = factor(v4_opcrit,
        levels = c(-999, 1, 2, 3, 4, 5),
        labels = c(
          "Not estimable",
          "single episode with good remission",
          "multiple episodes with good remission between episodes",
          "multiple episodes with partial remission between episodes",
          "ongoing chronic disease",
          "ongoing chronic disease with deterioration"
        )
      ),
      scid_suic_ide = factor(scid_suic_ide,
        levels = c(1, 2, 3, 4),
        label = c(
          "only fleeting thoughts",
          "serious thoughts (details were elaborated)",
          "persistent thoughts (serious and less serious)",
          "serious and persistent thoughts"
        )
      ),
      scid_suic_thght_mth = factor(scid_suic_thght_mth,
        levels = c(1, 2, 3),
        label = c(
          "no",
          "yes, but no details",
          "yes, including details"
        )
      ),
      v1_lftm_alc_dep = as.factor(v1_lftm_alc_dep),
      curr_paid_empl = factor(curr_paid_empl, level = c("N", "Y")),
      partner = factor(partner, level = c("N", "Y")),
      v3_cts_els_dic = as.factor(v3_cts_els_dic),
      v3_cts_1_dic = as.factor(v3_cts_1_dic),
      v3_cts_2_dic = as.factor(v3_cts_2_dic),
      v3_cts_3_dic = as.factor(v3_cts_3_dic),
      v3_cts_4_dic = as.factor(v3_cts_4_dic),
      v3_cts_5_dic = as.factor(v3_cts_5_dic),
      v1_ever_smkd = as.factor(v1_ever_smkd),
      v1_evr_ill_drg = as.factor(v1_evr_ill_drg),
      cur_psy_trm = factor(cur_psy_trm,
        levels = c(1, 2, 3, 4),
        labels = c("no", "yes, outpatient", "yes, day patient", "yes, inpatient")
      ),
      v1_scid_dsm_dx_cat = as.factor(v1_scid_dsm_dx_cat),
      v1_alc_pst12_mths = factor(v1_alc_pst12_mths,
        levels = as.character(c(1:7)),
        labels = c(
          "never", "only on special occasions",
          "once per month or less",
          "two to four times per month",
          "two to three times per week",
          "four times per week or several times but not daily",
          "daily"
        )
      ),
      v1_1st_ep = as.factor(v1_1st_ep),
      v1_sex = sex,
      v1_age_smk = as.numeric(v1_age_smk)
    ) %>%
    select(-v1_stat, -v1_scid_dsm_dx, -sex)

  PRS_treatment <- read.table(PRS_TREATMENT, header = T) %>%
    filter(FID %in% df$gsa_id)
  PRS_pfactor <- read.table(PRS_PFACTOR, header = T) %>%
    filter(FID %in% df$gsa_id)
  for (i in 1:nrow(df)) {
    for (j in 1:nrow(PRS_treatment)) {
      if (!is.na(df$gsa_id[i]) && df$gsa_id[i] == PRS_treatment$FID[j]) {
        df$PRS_TREATMENT_RES[i] <- PRS_treatment$SCORESUM[j]
        break
      }
    }
    for (j in 1:nrow(PRS_pfactor)) {
      if (!is.na(df$gsa_id[i]) && df$gsa_id[i] == PRS_pfactor$FID[j]) {
        df$PRS_PFACTOR_RES[i] <- PRS_pfactor$SCORESUM[j]
        break
      }
    }
  }

  # Add BMI
  load(BMI)
  df$BMI <- NA
  for (i in 1:nrow(df)) {
    for (j in 1:nrow(cluster_assignments_bmi)) {
      if (df$ID[i] == cluster_assignments_bmi$ID[j]) {
        if (df$visit[i] == 1) {
          df$BMI[i] <- cluster_assignments_bmi$v1_bmi[j]
          break
        } else if (df$visit[i] == 2) {
          df$BMI[i] <- cluster_assignments_bmi$v2_bmi[j]
          break
        } else if (df$visit[i] == 3) {
          df$BMI[i] <- cluster_assignments_bmi$v3_bmi[j]
          break
        } else if (df$visit[i] == 4) {
          df$BMI[i] <- cluster_assignments_bmi$v4_bmi[j]
          break
        }
      }
    }
  }
  return(df)
}


ForcePublicationClusters <- function(ComputedCluster, PublicationCluster) {
  load(PublicationCluster)
  if (identical(ComputedCluster$ID, CLUSTER$ID)) {
    ComputedCluster$PANSS_YMRS_IDSC_COGNITION_FUNCTIONNING <- CLUSTER$PANSS_YMRS_IDSC_COGNITION_FUNCTIONNING
  } else {
    for (i in 1:nrow(ComputedCluster)) {
      for (j in 1:nrow(CLUSTER)) {
        if (ComputedCluster$ID[i] == CLUSTER$ID[j]) {
          ComputedCluster$PANSS_YMRS_IDSC_COGNITION_FUNCTIONNING[i] <- CLUSTER$PANSS_YMRS_IDSC_COGNITION_FUNCTIONNING[j]
          break
        }
      }
    }
  }
  return(ComputedCluster)
}

# PROTEOME

ProtReadNPX <- function(NPX, CLUSTER, PSYCOURSE_NON_IMP) {
  load(PSYCOURSE_NON_IMP)
  pat_nonimp_proposal %<>%
    group_by(v1_id) %>%
    mutate(visit = row_number()) %>%
    ungroup() %>%
    filter(visit == 1)

  data <- OlinkAnalyze::read_NPX(NPX) %>%
    base::merge(
      CLUSTER %>%
        filter(visit == 1) %>%
        select(ID, PANSS_YMRS_IDSC_COGNITION_FUNCTIONNING) %>%
        rename(cluster = PANSS_YMRS_IDSC_COGNITION_FUNCTIONNING),
      by.x = "SampleID", by.y = "ID"
    ) %>%
    base::merge(
      pat_nonimp_proposal %>%
        select(v1_id, v1_scid_dsm_dx_cat, v1_sex, age, v4_opcrit),
      by.x = "SampleID", by.y = "v1_id"
    ) %>%
    droplevels() %>%
    mutate(
      cluster = factor(cluster, levels = c("A", "B"), labels = c("Low Severity", "High Severity")),
      diagnosis = as.factor(v1_scid_dsm_dx_cat),
      sex = as.factor(v1_sex),
      v1_sex = NULL,
      v1_scid_dsm_dx_cat = NULL
    )

  # non_variable_proteins <- data %>%
  #  group_by(Assay) %>%
  #  summarize(variance = var(NPX, na.rm = TRUE)) %>%
  #  filter(!is.na(variance)) %>%
  #  filter(variance >= quantile(variance)[2])
  # data %<>% filter(Assay %in% non_variable_proteins$Assay)
  return(data)
}







ComputeOlinkCor <- function(NPX, NPX_LIMMA, CLUSTER) {
  # for debugging
  # NPX <- tar_read(olink_npx)
  # NPX_LIMMA <- tar_read(olink_limma)
  # CLUSTER <- tar_read(cluster_assignments2)
  column_to_keep <- sapply(CLUSTER, is.numeric) | sapply(CLUSTER, is.integer)
  CLUSTER %<>%
    filter(visit == 1, ID %in% NPX$SampleID) %>%
    column_to_rownames("ID") %>%
    select(-visit) %>%
    select(any_of(names(column_to_keep)[column_to_keep])) %>%
    select(starts_with("nrpsy"), age)

  SIGNIF <- NPX_LIMMA %>% filter(qvalue <= 0.05)
  NPX %<>%
    select(SampleID, Assay, NPX) %>%
    filter(!is.na(NPX)) %>%
    reshape(idvar = "SampleID", timevar = "Assay", direction = "wide") %>%
    remove_rownames() %>%
    column_to_rownames("SampleID")
  colnames(NPX) %<>% str_remove("NPX\\.")


  # summary(rownames(NPX) == rownames(CLUSTER))
  # Mode    TRUE
  # logical     176

  PROTEINS <- c()
  VARIABLE <- c()
  CORRELATION <- c()
  PARTIAL_COR_C <- c()
  PARTIAL_COR_P <- c()
  PVALUE <- c()
  PROTEINS2 <- c()
  VARIABLE2 <- c()
  NPX_VALUE <- c()
  ASSESSMENT <- c()

  for (PROT in colnames(NPX)) {
    for (VAR in colnames(CLUSTER)) {
      CLUSTER2 <- merge(CLUSTER, NPX, by = 0) %>%
        filter(!is.na(!!sym(VAR)))

      PROTEINS %<>% append(PROT)
      VARIABLE %<>% append(VAR)
      cor <- cor.test(CLUSTER %>% pull(VAR),
        NPX %>% pull(PROT),
        method = "spearman"
      )
      partial_cor <- ppcor::pcor.test(CLUSTER2 %>% pull(VAR),
        CLUSTER2 %>% pull(PROT),
        CLUSTER2$age,
        method = "spearman"
      )
      CORRELATION %<>% append(cor$estimate)
      PARTIAL_COR_C %<>% append(partial_cor$estimate)
      PARTIAL_COR_P %<>% append(partial_cor$p.value)
      PVALUE %<>% append(cor$p.value)
    }
  }
  # Create the dataframe with the correlation and their pvalue.
  COR_DATA <- data.frame(
    PROTEIN = PROTEINS,
    VARIABLE = VARIABLE,
    COR = CORRELATION,
    PVALUE = PVALUE,
    PARTIAL_COR = PARTIAL_COR_C,
    PARTIAL_P = PARTIAL_COR_P
  ) %>%
    mutate(
      PADJ = p.adjust(PVALUE, method = "BH"),
      PARTIAL_PADJ = p.adjust(PARTIAL_P, method = "BH"),
      LIMMA = case_when(PROTEIN %in% row.names(SIGNIF) ~ "Significant", .default = "Not significant")
    )
  return(COR_DATA)
}

ComputeOlinkCorTreatment <- function(NPX, NPX_LIMMA, CLUSTER) {
  # for debugging
  # NPX <- tar_read(olink_npx)
  # NPX_LIMMA <- tar_read(olink_limma)
  # CLUSTER <- tar_read(cluster_assignments2)
  column_to_keep <- sapply(CLUSTER, is.numeric) | sapply(CLUSTER, is.integer)
  CLUSTER %<>%
    filter(visit == 1, ID %in% NPX$SampleID) %>%
    column_to_rownames("ID") %>%
    select(-visit) %>%
    select(any_of(c("Antidepressants", "Antipsychotics", "Mood_stabilizers", "Tranquilizers")))

  SIGNIF <- NPX_LIMMA %>% filter(qvalue <= 0.05)
  NPX %<>%
    select(SampleID, Assay, NPX) %>%
    filter(!is.na(NPX)) %>%
    reshape(idvar = "SampleID", timevar = "Assay", direction = "wide") %>%
    remove_rownames() %>%
    column_to_rownames("SampleID")
  colnames(NPX) %<>% str_remove("NPX\\.")


  # summary(rownames(NPX) == rownames(CLUSTER))
  # Mode    TRUE
  # logical     176

  PROTEINS <- c()
  VARIABLE <- c()
  CORRELATION <- c()
  PVALUE <- c()


  for (PROT in colnames(NPX)) {
    for (VAR in colnames(CLUSTER)) {
      CLUSTER2 <- merge(CLUSTER, NPX, by = 0) %>%
        filter(!is.na(!!sym(VAR)))
      PROTEINS %<>% append(PROT)
      VARIABLE %<>% append(VAR)
      cor <- cor.test(CLUSTER %>% pull(VAR),
        NPX %>% pull(PROT),
        method = "spearman"
      )
      CORRELATION %<>% append(cor$estimate)
      PVALUE %<>% append(cor$p.value)
    }
  }
  # Create the dataframe with the correlation and their pvalue.
  COR_DATA <- data.frame(
    PROTEIN = PROTEINS,
    VARIABLE = VARIABLE,
    COR = CORRELATION,
    PVALUE = PVALUE
  ) %>%
    mutate(
      PADJ = p.adjust(PVALUE, method = "BH"),
      LIMMA = case_when(PROTEIN %in% row.names(SIGNIF) ~ "Significant", .default = "Not significant")
    )
  return(COR_DATA)
}


OlinkLimma <- function(NPX_LONG) {
  # NPX_LONG <- tar_read(olink_npx)
  df_wide <- NPX_LONG %>%
    filter(!is.na(NPX)) %>%
    select(SampleID, Assay, NPX) %>%
    pivot_wider(names_from = SampleID, values_from = NPX)
  # Store protein names as rownames
  olink_matrix <- as.data.frame(df_wide)
  rownames(olink_matrix) <- olink_matrix$Assay
  olink_matrix %<>%
    mutate(Assay = NULL)
  metadata <- NPX_LONG %>% filter(!duplicated(SampleID))
  # Ensure your sample order in metadata matches the columns of olink_matrix
  olink_matrix <- olink_matrix[, metadata$SampleID]
  ## Create the design matrix with cluster and age
  metadata$cluster <- factor(metadata$cluster) # ensure it's a factor
  # limma analysis
  require(limma)
  require(variancePartition)
  require(edgeR)
  require(ashr)
  metadata <- metadata %>%
    dplyr::filter(SampleID %in% colnames(olink_matrix)) %>%
    dplyr::arrange(match(SampleID, colnames(olink_matrix)))
  # Ensure sample IDs match
  stopifnot(all(metadata$sample_id == colnames(olink_matrix)))
  # Convert to ExpressionSet
  expr_set <- Biobase::ExpressionSet(assayData = as.matrix(olink_matrix))
  # Define formula (example)
  # varPart <- fitExtractVarPartModel(expr_set, form, metadata)
  # Include covariates that matter
  design <- model.matrix(~ age + cluster, data = metadata)
  colnames(design)[3] <- "cluster"
  # Fit model and contrast
  fit2 <- lmFit(olink_matrix, design)
  fit2 <- limma::eBayes(fit2)
  top <- limma::topTable(fit2, coef = "cluster", number = Inf)
  top$betahat <- top$logFC
  top$sebehat <- sqrt(fit2$s2.post) * (fit2$stdev.unscaled %>%
    as.data.frame() %>%
    pull("cluster"))

  ash_fit <- ashr::ash(
    betahat = top$logFC,
    sebetahat = sqrt(fit2$s2.post) * (fit2$stdev.unscaled %>%
      as.data.frame() %>%
      pull("cluster"))
  )
  ash_res <- ash_fit$result
  # we apply an additional filter to avoid low psterior uncertainty (PosteriorSD < 0.02) to eliminate findings that results only from overconfidend shrinkage (https://pmc.ncbi.nlm.nih.gov/articles/PMC5379932/).
  ash_res_2 <- merge(ash_res, top %>% select(logFC, t, P.Value, adj.P.Val), by = 0) %>%
    filter(PosteriorSD > 0.02) %>%
    column_to_rownames("Row.names")
  return(ash_res_2)
  # ash_res_2 %>%
  #  mutate(SIGNIF = qvalue <= 0.05) %>%
  #  ggplot(aes(logFC, -log10(qvalue), colour = SIGNIF, label = ifelse(SIGNIF, Row.names, ""))) +
  #  geom_point() +
  #  geom_label_repel()
  #  ggtitle("volcano of the qvalue")
}

ComputeOlinkGO <- function(RES, THRESHOLD = 0.05) {
  genes <- NA
  if (is.data.frame(RES)) {
    genes <- RES %>%
      filter(qvalue <= THRESHOLD) %>%
      row.names()
  } else if (is.vector(RES)) {
    genes <- RES
  } else {
    stop(paste0("The RES argument must be a dataframe with a 'qvalue' column and genes names as row names, or a vector of genes names!"))
  }
  require(org.Hs.eg.db)
  ls <- list()
  ls[["QVALUE_THRESHOLD"]] <- THRESHOLD
  ls[["enrichGO"]] <- clusterProfiler::enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  # ls[["GO_barplot"]] <- barplot(ls[["enrichGO"]], showCategory = 10, title = "Top 10 GO Biological Processes")
  # ls[["GO_dotplot"]] <- dotplot(ls[["enrichGO"]], showCategory = 10, title = "Top 10 GO Biological Processes")
  # ls[["GO_cnplot"]] <- enrichplot::cnetplot(ls[["enrichGO"]], showCategory = 10)
  # ls[["GO_emapplot"]] <- enrichplot::emapplot(ls[["enrichGO"]], showCategory = 10)
  list(result = ls)
}

CorrectGOList <- function(GO_LIST) {
  if (!is.list(GO_LIST)) {
    error("The input of the CorrectGOList function must be a list object!")
  }
  ls <- list()
  for (NAMES in names(GO_LIST)) {
    ls[[paste0("q_", GO_LIST[[NAMES]][["QVALUE_THRESHOLD"]])]] <- GO_LIST[[NAMES]][["enrichGO"]]
  }
  return(ls)
}
