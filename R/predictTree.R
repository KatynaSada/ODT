#' Predict Treatment Outcomes with a Trained Decision Tree
#'
#' This function utilizes a trained decision tree model (ODT) to predict treatment
#' outcomes for test data based on patient sensitivity data and features, such as
#' mutations or gene expression profiles.
#'
#' @param tree A trained decision tree object created by the `trainTree` function.
#' @param PatientData A matrix representing patient features, where rows correspond to patients/samples
#'                    and columns correspond to genes/features. This matrix can contain:
#'                    \itemize{
#'                      \item Binary mutation data (e.g., presence/absence of mutations).
#'                      \item Continuous data from gene expression profiles (e.g., expression levels).
#'                    }
#' @param PatientSensitivityTrain A matrix containing the drug response values of the **training dataset**. 
#'                                 In this matrix, rows correspond to patients, and columns correspond to drugs. 
#'                                 This matrix is used solely for extracting treatment names and is not 
#'                                 used in the prediction process itself.
#' @return A factor representing the assigned treatment for each node in the
#' decision tree based on the provided patient data and sensitivity.
#'
#' @examples
#' \dontrun{
#'   # Example 1: Prediction using mutation data
#'   data(DataODT.rda)
#'   ODTmut <- trainTree(PatientData = mutations_w12, 
#'                       PatientSensitivityTrain = drug_response_w12, 
#'                       minbucket = 10)
#'   ODT_mutpred <- predictTree(tree = ODTmut, 
#'                               PatientSensitivityTrain = drug_response_w34, 
#'                               PatientData = mutations_w34)
#'
#'   # Example 2: Prediction using gene expression data
#'   data(DataODT.rda)
#'   ODTExp <- trainTree(PatientData = expression_w34, 
#'                        PatientSensitivityTrain = drug_response_w34, 
#'                        minbucket = 20)
#'   ODT_EXPpred <- predictTree(tree = ODTExp, 
#'                               PatientSensitivityTrain = drug_response_w12, 
#'                               PatientData = expression_w12)
#' }
#'
#' @import partykit
#' @export
predictTree <- function(tree, PatientData, PatientSensitivityTrain) {
  # Check if tree is of the correct class
  if (!inherits(tree, "party")) {
    stop("The 'tree' parameter must be a trained decision tree object of class 'party'.")
  }
  
  # Adjust PatientData based on its unique values
  if (length(unique(c(unlist(PatientData)))) == 2) {
    PatientData <- PatientData - min(PatientData) + 1L
    mode(PatientData) <- "integer"
  } 
  
  # Predict treatments based on the decision tree
  treatments <- unlist(nodeapply(tree,
                                 predict.party(tree, as.data.frame(PatientData)), 
                                 info_node))
  
  # Match treatments with sensitivity data
  TratamientoTree <- match(treatments, colnames(PatientSensitivityTrain))
  TratamientoTree <- factor(TratamientoTree, levels = 1:ncol(PatientSensitivityTrain))
  
  return(TratamientoTree)
}
