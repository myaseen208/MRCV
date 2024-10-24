#' @name    farmer1
#' @aliases farmer1
#' @docType data
#' @title Data for One SRCV and One MRCV from the Kansas Farmer Survey
#' @description
#' Responses for one SRCV and one MRCV from a survey of Kansas farmers. This
#' data was first examined by Loughin and Scherer (1998) and subsequently
#' examined by papers such as Agresti and Liu (1999) and Bilder, Loughin, and
#' Nettleton (2000).
#'
#' @format The data frame contains the following 6 columns:
#'
#' Column 1, labeled \code{Ed}, corresponds to the respondent's highest
#' attained level of education.  A total of r = 5 levels of education are
#' represented in the data. \itemize{ \item\code{1}: High school \item\code{2}:
#' Vocational school \item\code{3}: 2-Year college \item\code{4}: 4-Year
#' college \item\code{5}: Other } Columns 2-6 correspond to the response
#' categories for the question "What are your primary sources of veterinary
#' information?"  Binary responses (1 = Yes, 0 = No) are provided for each
#' category. \itemize{ \item\code{Y1}: Professional consultant \item\code{Y2}:
#' Veterinarian \item\code{Y3}: State or local extension service
#' \item\code{Y4}: Magazines \item\code{Y5}: Feed companies and representatives
#' }
#' @references Agresti, A. and Liu, I.-M. (1999) Modeling a categorical
#' variable allowing arbitrarily many category choices.  \emph{Biometrics},
#' \bold{55}, 936--943.
#'
#' Bilder, C., Loughin, T., and Nettleton, D. (2000) Multiple marginal
#' independence testing for pick any/c variables.  \emph{Communications in
#' Statistics--Theory and Methods}, \bold{29}, 1285--1316.
#'
#' Loughin, T. M. and Scherer, P. N. (1998) Testing for association in
#' contingency tables with multiple column responses.  \emph{Biometrics},
#' \bold{54}, 630--637.
#' @source Loughin, T. M. and Scherer, P. N. (1998) Testing for association in
#' contingency tables with multiple column responses.  \emph{Biometrics},
#' \bold{54}, 630--637.
#' @examples
#' data(farmer1)
NULL
