#' @name loadBWfile
#' @title loadBWfile
#' @description loadBWfile function used to read bigwig files
#' @param file Your bigwig file name.
#' @param sample Your bigwig file sample name.
#' @param group1 Your bigwig file group name, like 'Input','IP'.
#' @param group2 Your bigwig file treat name, like 'control','treat1','treat2'.
#'
#' @return Reture a dataframe all bigwig contents.
#' @export
#'
#' @examples
#' \dontrun{
#' file <- c(
#'   "control.input.bw", "control.ip.bw",
#'   "t1.input.bw", "t1.ip.bw",
#'   "t2.input.bw", "t2.ip.bw"
#' )
#' samp <- c(
#'   "control.input", "control.ip",
#'   "t1.input", "t1.ip",
#'   "t2.input", "t2.ip"
#' )
#' group1 <- c(rep(c("Input", "IP"), 3))
#' group2 <- c(rep("Control", 2), rep("T1", 2), rep("T2", 2))
#'
#' # 2.load bw files
#' allBw <- loadBWfile(file = file, sample = samp, group1 = group1, group2 = group2)
#' }
#'
# define viriables
globalVariables(c("strand", "width"))


# define function2
loadBWfile <- function(file = NULL,
                       sample = NULL,
                       group1 = NULL,
                       group2 = NULL) {
  # loop for read all files
  purrr::map_df(1:length(file), function(x) {
    # import single bigwig file
    bw <- rtracklayer::import.bw(file[x]) %>%
      data.frame() %>%
      dplyr::select(-strand, -width)

    # add name
    bw$sample <- sample[x]

    # add meta info

    # add group
    bw$group1 <- group1[x]

    # add sample name
    bw$group2 <- group2[x]

    return(bw)
  }) -> all_bw
}
