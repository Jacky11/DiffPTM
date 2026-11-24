#################################################
# Filters on the number of quantified values
#
# For each condition k (n_k samples), obs_k = number of quantified values of a
# protein (value > 0 with enough unique peptides, see step 1 of the pipeline).
#
#   required_k = max(min_observed, ceiling(min_percent_observed / 100 * n_k))
#   sufficiently quantified in k : obs_k >= required_k
#   absent from k                : obs_k <= absent_max
#   present in k                 : obs_k >= present_min_k (default: required_k)
#
# Categories:
#   present1 / present2 : absent from one condition and present in the other
#                         (not tested, reported as more abundant where present)
#   tested              : not absent from any condition, and sufficiently quantified
#                         in at least one condition (quantified_in = "one") or in
#                         each condition (quantified_in = "each")
#   removed             : all other proteins
#################################################

filter_thresholds <- function(p, n1, n2) {
  req <- c(max(p$min_observed, ceiling(p$min_percent_observed / 100 * n1 - 1e-9)),
           max(p$min_observed, ceiling(p$min_percent_observed / 100 * n2 - 1e-9)))
  pm <- if (is.null(p$present_min) || is.na(p$present_min)) req else rep(p$present_min, 2)
  am <- p$absent_max
  if (am >= min(req))
    stop("absent_max (", am, ") must be lower than the number of quantified values required in a condition (",
         min(req), ").", call. = FALSE)
  if (any(pm <= am))
    stop("present_min must be greater than absent_max.", call. = FALSE)
  list(required = req, present_min = pm, absent_max = am, n = c(n1, n2), mode = p$quantified_in)
}

classify_proteins <- function(obs1, obs2, th, c1, c2) {
  req <- th$required; pm <- th$present_min; am <- th$absent_max
  Q1 <- obs1 >= req[1]; Q2 <- obs2 >= req[2]
  A1 <- obs1 <= am;     A2 <- obs2 <= am
  present1 <- A2 & !A1 & obs1 >= pm[1]
  present2 <- A1 & !A2 & obs2 >= pm[2]
  tested <- !A1 & !A2 & (if (th$mode == "each") Q1 & Q2 else Q1 | Q2)
  category <- ifelse(tested, "tested", ifelse(present1, "present1", ifelse(present2, "present2", "removed")))

  absent_txt <- function(cond) if (am == 0) paste0("no quantified value in ", cond) else
    paste0("at most ", n_values(am), " in ", cond)
  reason <- rep(NA_character_, length(obs1))
  for (i in which(category == "removed")) {
    reason[i] <- if (A1[i] && A2[i]) {
      if (am == 0) paste0("no quantified value in ", c1, " nor in ", c2) else
        paste0("at most ", n_values(am), " in ", c1, " and in ", c2)
    } else if (A2[i]) {
      paste0(absent_txt(c2), " and only ", n_values(obs1[i]), " in ", c1, " (", pm[1], " required to be considered present in ", c1, ")")
    } else if (A1[i]) {
      paste0(absent_txt(c1), " and only ", n_values(obs2[i]), " in ", c2, " (", pm[2], " required to be considered present in ", c2, ")")
    } else if (th$mode == "each") {
      paste(c(if (!Q1[i]) paste0("only ", n_values(obs1[i]), " in ", c1, " (", req[1], " required)"),
              if (!Q2[i]) paste0("only ", n_values(obs2[i]), " in ", c2, " (", req[2], " required)")),
            collapse = "; ")
    } else {
      paste0("fewer than ", req[1], " quantified values in ", c1, " (", obs1[i], ") and fewer than ",
             req[2], " in ", c2, " (", obs2[i], ")")
    }
  }
  list(category = category, reason = reason)
}

# "1 value" / "2 values"
n_values <- function(k, what = "quantified value") paste0(k, " ", what, if (k > 1) "s" else "")

# Methods text describing the filters with the values used in the analysis
filter_methods_text <- function(p, th, c1, c2) {
  req <- th$required; n <- th$n; pm <- th$present_min; am <- th$absent_max
  rule <- if (p$min_percent_observed > 0)
    paste0(" (i.e. at least ", p$min_observed, " values and at least ", p$min_percent_observed, "% of the samples of the condition)") else ""
  suff <- if (req[1] == req[2] && n[1] == n[2]) {
    paste0("A protein is considered sufficiently quantified in a condition if it has at least ",
           n_values(req[1]), " out of the ", n[1], " samples of this condition", rule, ".")
  } else {
    paste0("A protein is considered sufficiently quantified in ", c1, " if it has at least ", n_values(req[1]),
           " out of ", n[1], " samples, and in ", c2, " if it has at least ", n_values(req[2]), " out of ", n[2], " samples", rule, ".")
  }
  absent <- if (am == 0) "no quantified value" else paste0("at most ", n_values(am))
  pm_txt <- if (pm[1] == pm[2]) paste0("at least ", n_values(pm[1])) else
    paste0("at least ", n_values(pm[1]), " (", c1, ") or ", n_values(pm[2]), " (", c2, ")")
  list(
    definition = paste0("A quantified value is a non-zero intensity measured with at least ",
                        p$min_unique_peptides, " unique peptide", if (p$min_unique_peptides > 1) "s" else "", "."),
    tested = paste0(suff, " A protein is tested if it is sufficiently quantified in ",
                    if (th$mode == "each") "each of the two conditions" else "at least one of the two conditions",
                    " and has at least ", n_values(am + 1), " in each condition."),
    present = paste0("A protein with ", absent, " in one condition and ", pm_txt, " in the other condition is considered ",
                     "present in the latter and absent from the former. These proteins are not tested; they are reported ",
                     "as more abundant in the condition where they are present."),
    removed = "All other proteins are removed (sheet without_enough_values, with the reason for each protein)."
  )
}
