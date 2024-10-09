#' Predict Treatment Outcomes with a Trained Decision Tree
#'
#' This function utilizes a trained decision tree model (ODT) to predict treatment
#' outcomes for test data based on patient sensitivity data and features, such as
#' mutations or gene expression profiles.
#'
#' @param tree A trained decision tree object created by the `trainTree` function.
#' @param PatientSensitivity A numeric vector representing drug response values.
#' Higher values indicate greater resistance and, consequently, lower sensitivity
#' to treatment. Depending on the interpretation of these values, users may need to
#' adjust the sign of this data.
#' @param PatientData A matrix of patient features. This may contain:
#' \itemize{
#'   \item Binary mutation data (1 indicates the presence of a mutation).
#'   \item Continuous data from gene expression profiles.
#' }
#'
#' @return A factor representing the assigned treatment for each node in the
#' decision tree based on the provided patient data and sensitivity.
#'
#' @examples
#' \dontrun{
#'   # Example 1: Prediction using mutation data
#'   data(DataODT.rda)
#'   ODTmut <- trainTree(PatientData = mutations_w12, 
#'                       PatientSensitivity = drug_response_w12, 
#'                       minbucket = 10)
#'   ODT_mutpred <- predictTree(tree = ODTmut, 
#'                               PatientSensitivity = drug_response_w34, 
#'                               PatientData = mutations_w34)
#'
#'   # Example 2: Prediction using gene expression data
#'   data(DataODT.rda)
#'   ODTExp <- trainTree(PatientData = expression_w34, 
#'                        PatientSensitivity = drug_response_w34, 
#'                        minbucket = 20)
#'   ODT_EXPpred <- predictTree(tree = ODTExp, 
#'                               PatientSensitivity = drug_response_w12, 
#'                               PatientData = expression_w12)
#' }
#'
#' @import partykit
#' @export
predictTree <- function(tree, PatientSensitivity, PatientData) {
  # Check if tree is of the correct class
  if (!inherits(tree, "party")) {
    stop("The 'tree' parameter must be a trained decision tree object of class 'party'.")
  }
  
  # Adjust PatientData based on its unique values
  if (length(unique(c(unlist(PatientData)))) == 2) {
    PatientData <- PatientData - min(PatientData) + 1L
    mode(PatientData) <- "integer"
  } else {
    PatientData <- t(PatientData)
  }
  
  # Predict treatments based on the decision tree
  treatments <- unlist(nodeapply(tree,
                                 predict.party(tree, as.data.frame(PatientData)), 
                                 info_node))
  
  # Match treatments with sensitivity data
  TratamientoTree <- match(treatments, colnames(PatientSensitivity))
  TratamientoTree <- factor(TratamientoTree, levels = 1:ncol(PatientSensitivity))
  
  return(TratamientoTree)
}
