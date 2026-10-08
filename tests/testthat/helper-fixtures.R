feature_fixture <- function() {
  tibble::tibble(
    WellId = paste0("W", 1:6),
    FeatureA = c(1, 1, 1, 4, 5, 6),
    FeatureB = c(-3, -2, -1, 1, 2, 3)
  )
}

metadata_fixture <- function() {
  tibble::tibble(
    WellId = paste0("W", 1:6),
    Plate = c("P1", "P1", "P1", "P2", "P2", "P2"),
    Compound = c("DrugA", "DrugA", NA, "DrugB", "DrugB", NA),
    Concentration = c(0.1, 1, NA, 0.1, 1, NA)
  )
}
