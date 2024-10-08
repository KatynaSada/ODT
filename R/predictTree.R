#' predictTree function
#' Predict Treatment Outcomes with a Trained Decision Tree
#'
#' This function predicts treatment outcomes for test data using a trained
#' decision tree model (ODT). The predictions are based on the provided
#' patient sensitivity data and patient features (mutations or gene expression).
#'
#' @param tree A trained decision tree object.
#' @param PatientSensitivity A numeric vector representing drug response values.
#' Higher values indicate greater resistance and thus lower sensitivity to the treatment.
#' Depending on the interpretation of these response values, users may need to 
#' adjust the sign of the data accordingly.
#' @param PatientData A matrix containing patient features. This can represent
#' either binary mutation data (where 1 indicates the presence of a mutation) 
#' or continuous data from gene expression profiles.
#'
#' @return A party object representing the predicted tree, with treatments assigned 
#' to each node based on the provided patient data and sensitivity.
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
#' @export
predictTree <- function(tree, PatientSensitivity, PatientData) {
  if (length(unique(c(unlist(PatientData)))) == 2) {
    PatientData <- PatientData - min(PatientData) + 1L
    PatientData <- apply(PatientData, 2, as.integer)
  } else{
    PatientData <- t(PatientData)
  }
  
  treatments <- unlist(nodeapply(tree,
                                 predict.party(tree, as.data.frame(PatientData)), 
                                 info_node))
  TratamientoTree <- match(treatments, colnames(PatientSensitivity))
  TratamientoTree <- factor(TratamientoTree, levels = 1:ncol(PatientSensitivity))
  return(TratamientoTree)
}