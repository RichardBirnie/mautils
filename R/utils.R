# A collection of utility functions

#' Take a data frame and convert all factor variables to character variables
#'
#' @param df A data frame with at least one character column
#' @return The same data frame but with factors converted to character
#' @examples
#' data(iris)
#' class(iris$Species)
#' [1] "factor"
#' iris = factorToCharacter(iris)
#' class(iris$Species)
#' [1] "character"
factorToCharacter = function(df) {
  df[sapply(df, is.factor)] <- lapply(df[sapply(df, is.factor)],
                                      as.character)
  df
}

#' Capitalise the first letter of each word in a string
#'
#' @param s A vector of one or more character strings
#' @return A vector of one or more character strings with the first letter of each word capitalised
#' @examples
#' # A single character string with two words
#' capwords('one two')
#' # Three character strings. The first has two words the last two
#'  have one word each
#' capwords(c('one two', 'three', 'four'))
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

#' A simple wrapper around \code{write.xlsx} from the \code{xlsx} package
#'
#' @param x A \code{data.frame} to be saved as a single sheet in an excel file
#' @param file the path to the output file. This should include the '.xlsx' file
#'   extension
#' @param sheetName a character string with the sheet name.
#' @param col.names a logical value indicating if the column names of \code{x}
#'   are to be written along with \code{x} to the file.
#' @param row.names a logical value indicating whether the row names of \code{x}
#'   are to be written along with \code{x} to the file.
#' @param append a logical value indicating if \code{x} should be appended to an
#'   existing file. If \code{TRUE} the file is read from disk.
#' @param showNA a logical value. If set to FALSE, NA values will be left as
#'   empty cells
#' @param overwrite a logical value indicating if \code{x} should overwrite the
#'   contents of an existing sheet in the workbook
#'
#' @details This function is a basic wrapper around \code{write.xlsx} from the
#'   \code{xlsx} package designed to modify some of the default behaviour. The
#'   main difference is the addition of the \code{overwrite} argument. This is
#'   intended to allow the user to overwrite sheets that already exist in the
#'   workbook. If the Excel file specified by \code{file} already exists and
#'   includes a worksheet with same name as \code{sheetName} then this sheet
#'   will be deleted and replaced if \code{overwrite=TRUE}. If
#'   \code{overwrite=FALSE} then the sheet is appended to the end of the file.
#'
#'   The default for \code{append} in the underlying \code{write.xlsx} function
#'   is \code{FALSE}. If the file you want to write to already exists you should
#'   set this to \code{append=TRUE}. If the file exists and \code{append=FALSE}
#'   the existing file will be deleted and replaced. In most cases
#'   \code{append=TRUE} is probably what you want.
#'
#' @examples
#' #basic usage for a new file
#' #creates a temporary file in the current directory
#' data(iris)
#' f = paste0('./temp','.xlsx')
#' saveXLSX(x = iris, file = f, sheetName = 'Test')
#'
#' #using the same file as above
#' #use a different data set but keep the same name
#' #if the file already exists set overwrite = TRUE and append = TRUE
#' data(CO2)
#' saveXLSX(x = CO2, file = f, sheetName = 'Test', append = TRUE,
#' overwrite=TRUE)
#'
#' #if overwrite = FALSE then the new sheet is appended at
#' #the end of the file
#' #use another different data set.
#' data(Indometh)
#' saveXLSX(x = Indometh, file = f, sheetName = 'Test', append = TRUE,
#' overwrite=FALSE)
#'
#' @seealso \code{\link[xlsx]{write.xlsx}}
saveXLSX = function(x, file, sheetName = "Sheet1", col.names = TRUE, row.names = FALSE,
                    append = FALSE, showNA = FALSE, overwrite = TRUE) {

  #if the file exists load the workbook and check if it has a sheet with
  #the same name as 'sheetName'
  if (file.exists(file)) {
    wb = xlsx::loadWorkbook(file)
    sheets = xlsx::getSheets(wb)
    sheetExists = sheetName %in% names(sheets)
    if (sheetExists) {
      #if the worksheet already exists decide if we want to overwrite it
      if(overwrite) {
        #remove the existing sheet and replace it with the new one
        xlsx::removeSheet(wb, sheetName = sheetName)
        xlsx::saveWorkbook(wb, file)
        xlsx::write.xlsx(
          as.data.frame(x), file = file, sheetName = sheetName, col.names = col.names, row.names =
            row.names, showNA = showNA, append = append
        )
      } else {
        #drop the sheetname and append this at the end of the file
        sheet = paste0('Sheet', length(sheets) + 1)
        xlsx::write.xlsx(
          as.data.frame(x), file = file, sheetName=sheet, col.names = col.names, row.names =
            row.names, showNA = showNA, append = TRUE
        )
      }
    } else {
      xlsx::write.xlsx(
        as.data.frame(x), file = file, sheetName = sheetName, col.names = col.names, row.names =
          row.names, showNA = showNA, append = append
      )
    }
  } else {
    xlsx::write.xlsx(
      as.data.frame(x), file = file, sheetName = sheetName, col.names = col.names, row.names =
        row.names, showNA = showNA, append = append
    )
  }
}
