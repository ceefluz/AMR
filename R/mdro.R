

# Function to detect MDRO according to the Dutch guidelines (https://www.rivm.nl/dsresource?objectid=b6b99580-44e2-4b9c-8183-52871e61764f&type=org&disposition=inline)


mdro <- function(x) {
  
  # count "R" occurence
  
  count_resistant <- function(group_selected) {
    sum(stri_count(as.character(group_selected), fixed = "R"), na.rm = TRUE)
  }
  
  # filter for agents so only available columns are actually evaluated (based on ablist)
  
  fluoro <- ablist$umcg[ablist$atc_group2 == "Fluorochinolonen"][!is.na(ablist$umcg[ablist$atc_group2 == "Fluorochinolonen"])]
  
  fluoro_columns <- names(x)[names(x) %in% fluoro]
  
  amino <- ablist$umcg[ablist$atc_group1 == "Aminoglycosiden"][!is.na(ablist$umcg[ablist$atc_group1 == "Aminoglycosiden"])]
  
  amino_columns <- names(x)[names(x) %in% amino]
  
  sulfo_columns <- ablist$umcg[ablist$atc == "J01EE01"]
  
  pip_columns <- ablist$umcg[ablist$atc %in% c("J01CR05", "J01CA12")]
  
  cefta_columns <- ablist$umcg[ablist$atc == "J01DD02"]
  
  penicillines <- ablist$umcg[ablist$atc_group1 == "Betalactam-antibiotica, penicillines"][!is.na(ablist$umcg[ablist$atc_group1 == "Betalactam-antibiotica, penicillines"])]
  
  penicillines_columns <- names(x)[names(x) %in% penicillines]
  
  vanco_columns <- ablist$umcg[ablist$atc == "J01XA01"]
  
  
  # ESBL or CARB
  
  esbl_count <- if (!is.na(x$esbl)) {
    if (x$esbl == TRUE) {
      1
    }
    else {
      0
    }
  }
  else {
    0
  }
  
  carb_count <- if (!is.na(x$carb)) {
    if (x$carb == TRUE) {
      1
    } else {
      0
    }
  }
  else {
    0
  }
  
  group_a_count <- esbl_count + carb_count
  
  # Fluorochinolones
  
  fluoro_select <- select(x, .columns = fluoro_columns)
  
  fluoro_count <- count_resistant(fluoro_select)
  
  # Aminoglycosides
  
  amino_select <- select(x, .columns = amino_columns)
  
  amino_count <- count_resistant(amino_select)
  
  # Sulfonamides and trimethoprim
  
  sulfo_select <- select(x, .columns = sulfo_columns)
  
  sulfo_count <- count_resistant(sulfo_select)
  
  # Piperacilline
  
  pip_select <- select(x, .columns = pip_columns)
  
  pip_count <- count_resistant(pip_select)
  
  # Ceftazidime
  
  cefta_select <- select(x, .columns = cefta_columns)
  
  cefta_count <- count_resistant(cefta_select)
  
  
  ## Enterobacteriaceae
  
  if (x$family == "Enterobacteriaceae") {
    if (group_a_count >= 1) {
      return(TRUE)
    } else {
      if (fluoro_count >= 1 & amino_count >= 1) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    }
  }
  
  ## acinetobacter species
  
  if (x$fullname == "Acinetobacter species") {
    if (x$carb == TRUE) {
      TRUE
    } else {
      if (fluoro_count >= 2 & amino_count >= 1) { # >= 2 due to intrinsic resistancy against norfloxacine
        return(TRUE)
      } else {
        return(FALSE)
      }
    }
  }
  
  ## Stenotrophomonas maltophilia
  
  if (x$fullname == "Stenotrophomonas maltophilia") {
    if (sulfo_count == 1) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  
  ## Pseudomonas aeruginosa
  
  pseudo_fluoro <- if (fluoro_count >= 1) {
    1
  } else {
    0
  }
  pseudo_amino <- if (amino_count >= 1) {
    1
  } else {
    0
  }
  pseudo_carb <- if (x$carb == TRUE) {
    1
  } else {
    0
  }
  
  pseudo_sum <- sum(
    pseudo_fluoro,
    pseudo_amino,
    pseudo_carb,
    cefta_count,
    pip_count
  )
  
  if (x$fullname == "Pseudomonas aeruginosa") {
    if (pseudo_sum >= 3) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  
  
  ## Streptococcus pneumoniae
  
  if (x$fullname == "Streptococcus pneumoniae") {
    if (penicillines_columns >= 1 | vanco_columns >= 1) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  
  ## Enterococcus faecium
  
  if (x$fullname == "Enterococcus faecium") {
    if (penicillines_columns >= 1 & vanco_columns >= 1) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  
  else {
    return(FALSE)
  }
  
}
