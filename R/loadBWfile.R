#' @title loadBWfile
#' @description loadBWfile function used to read bigwig files.
#' @param file Your bigwig file name.
#' @param name Your bigwig file sample name.
#' @param samp Your bigwig file group name, like 'Input','IP'.
#' @param group Your bigwig file treat name, like 'control','treat1','treat2'.
#'
#' @return
#' @export
#'
#' @examples
#'\dontrun{
#'file <- c("control.input.bw" ,"control.ip.bw",
#'          "t1.input.bw","t1.ip.bw",
#'          "t2.input.bw","t2.ip.bw")
#'sname <- c("control.input" ,"control.ip",
#'           "t1.input","t1.ip",
#'           "t2.input","t2.ip")
#'samp <- c(rep(c('Input','IP'),3))
#'grou <- c(rep('Control',2),rep('T1',2),rep('T2',2))
#'
#'# 2.load bw files
#'allBw <- loadBWfile(file = file,name = sname,samp = samp,group = grou)
#'}

# define function2
loadBWfile <- function(file = NULL,
                       name = NULL,
                       samp = NULL,
                       group = NULL) {
  # loop for read all files
  map_df(1:length(file), function(x) {
    # import single bigwig file
    bw <- import.bw(file[x]) %>%
      data.frame() %>%
      select(-strand, -width)

    # add name
    bw$name <- name[x]

    # add meta info

    # add group
    bw$group <- group[x]

    # add sample name
    bw$sample <- samp[x]

    return(bw)
  }) -> all_bw
}
