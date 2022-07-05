#' @name plotTrack
#' @title plotTrack
#' @description plotTrack is used to visulaize gene track from bigwig files
#' @param gtfFile Your annotation file with GTF format, eg 'Mus_musculus.GRCm38.102.gtf'.
#' @param gene Which gene used to plot, eg 'Actb'.
#' @param bigwigFile loadBWfile function output results.
#' @param multiple Multiple gene isoforms will be ploted when it is TRUE, default is 'FALSE'.
#' @param myTransId Which transcript will be plot, can be multiple vectors. Parameter 'multiple' should be TRUE when you want to plot 2 or more transcripts.
#' @param strucCol Gene colors, defult is '#336699'.
#' @param arrowCol Arrow colors on gene structure, defult is '#336699'.
#' @param strucHeight Gene structure height, defult is 0.25.
#' @param uped Extend from start site, defult is 1000 bp.
#' @param downed Extend from end site, defult is 1000 bp.
#' @param addfacetCol Whether add more facets, defult is 'FALSE'.
#' @param facetFill Facets fill colors.
#' @param borderCol Facets border colors.
#' @param sampleAes Which variable should be used to aes track plot.
#' @param trackCol Supply your own colors vectors to change track colors.
#' @param facetVars Facet by this variable names. Can be multiple.
#' @param relHeight Relative height of track plot to the the gene structure, default is 10.
#' @param annoTextSize Facet text size, default is 10.
#' @param base_size ggplot base size, default is 12.
#' @param ylab Y axis label, default is 'Read coverage'.
#'
#' @return Return a track plot.
#' @export
#'
#' @examples
#' \dontrun{
#' # examples
#' plotTrack(
#'   gtfFile = gtf,
#'   gene = "Actb",
#'   arrowCol = "black",
#'   bigwigFile = allBw,
#'   sampleAes = "group1",
#'   facetVars = "sample"
#' )
#'
#' # change aes variable
#' plotTrack(
#'   gtfFile = gtf,
#'   gene = "Actb",
#'   arrowCol = "black",
#'   bigwigFile = allBw,
#'   sampleAes = "sample",
#'   facetVars = "sample"
#' )
#'
#' # change track colors
#' plotTrack(
#'   gtfFile = gtf,
#'   gene = "Actb",
#'   arrowCol = "black",
#'   bigwigFile = allBw,
#'   sampleAes = "sample",
#'   facetVars = "sample",
#'   trackCol = pal_lancet()(6)
#' )
#'
#' # draw multiple genes
#' plotTrack(
#'   gtfFile = gtf,
#'   gene = "Tnf",
#'   arrowCol = "black",
#'   bigwigFile = allBw,
#'   sampleAes = "group1",
#'   facetVars = "sample",
#'   multiple = TRUE,
#'   myTransId = c("ENSMUST00000025263", "ENSMUST00000167924")
#' )
#'
#' # change facet fill colors
#' plotTrack(
#'   gtfFile = gtf,
#'   gene = "Actb",
#'   arrowCol = "black",
#'   bigwigFile = allBw,
#'   sampleAes = "group1",
#'   addfacetCol = TRUE,
#'   facetVars = "sample",
#'   facetFill = pal_d3()(6),
#'   borderCol = rep("white", 6)
#' )
#'
#' # add one facet
#' plotTrack(
#'   gtfFile = gtf,
#'   gene = "Actb",
#'   arrowCol = "black",
#'   bigwigFile = allBw,
#'   sampleAes = "group1",
#'   facetVars = c("group2", "group1"),
#'   addfacetCol = TRUE,
#'   facetFill = c(pal_d3()(3), pal_npg()(6)),
#'   borderCol = rep("white", 9)
#' )
#'
#' # add more one facet
#' plotTrack(
#'   gtfFile = gtf,
#'   gene = "Actb",
#'   arrowCol = "black",
#'   bigwigFile = allBw,
#'   sampleAes = "group1",
#'   facetVars = c("group2", "group1", "sample"),
#'   addfacetCol = TRUE,
#'   facetFill = c(pal_d3()(3), pal_npg()(6), pal_locuszoom()(6)),
#'   borderCol = rep("white", 15)
#' )
#' }
#'
# define viriables
globalVariables(c(
  "cdsLen", "desc", "end", "exonLen", "gene_name", "score", "seqnames",
  "start", "strand", "transcript_id", "transcript_name", "type", "width"
))

# define function
plotTrack <- function(gtfFile = NULL,
                      gene = NULL,
                      bigwigFile = NULL,
                      multiple = FALSE,
                      myTransId = NULL,
                      strucCol = "#336699",
                      arrowCol = "#336699",
                      strucHeight = 0.25,
                      uped = 1000,
                      downed = 1000,
                      addfacetCol = FALSE,
                      facetFill = NULL,
                      borderCol = NULL,
                      sampleAes = NULL,
                      trackCol = NULL,
                      facetVars = NULL,
                      relHeight = 10,
                      annoTextSize = 10,
                      base_size = 12,
                      ylab = "Read coverage") {
  ############################################################
  # filter gene
  ginfo <- gtfFile %>% dplyr::filter(gene_name == gene)

  # transcript_id
  tid <- unique(ginfo$transcript_id) %>%
    stats::na.omit() %>%
    as.character()

  # exon and cds length
  purrr::map_df(tid, function(x) {
    tmp <- ginfo %>% dplyr::filter(transcript_id == x)

    # 1.transcript region
    rg <- tmp %>%
      dplyr::filter(type == "transcript") %>%
      dplyr::select(start, end)

    # 2.get exon length
    exonLen <- tmp %>%
      dplyr::filter(type == "exon") %>%
      dplyr::select(width) %>%
      sum()

    # 3.test gene whether is protein coding
    # and order by cds and exon length
    eType <- unique(tmp$type)
    if ("CDS" %in% eType) {
      cdsLen <- tmp %>%
        dplyr::filter(type == "CDS") %>%
        dplyr::select(width) %>%
        sum()
    } else {
      cdsLen <- 0
    }

    # 4.output
    tmpRes <- data.frame(
      gene = unique(tmp$gene_name),
      tid = x,
      exonLen = exonLen,
      cdsLen = cdsLen,
      chr = unique(tmp$seqnames),
      rg
    )
    return(tmpRes)
  }) %>%
    dplyr::arrange(desc(cdsLen), desc(exonLen)) -> lengthInfo
  print(lengthInfo)

  ############################################################
  # prepare gene strcture

  # 1.select which transcript to plot
  if (multiple == FALSE) {
    # filter gene
    target_gene <- gtfFile %>%
      dplyr::filter(transcript_id %in% lengthInfo$tid[1])

    # transcript name info
    print(paste0(
      "Choosed transcript ",
      lengthInfo$tid[1],
      " to draw gene structure!"
    ))
  } else if (multiple == TRUE) {
    # filter gene
    target_gene <- gtfFile %>%
      dplyr::filter(transcript_id %in% myTransId)
  }

  # 2.get exon info
  target_exon <- target_gene %>% dplyr::filter(type == "exon")

  # 3.get cds info
  target_cds <- target_gene %>% dplyr::filter(type == "CDS")

  # 4.plot
  gene_strcture <- ggplot2::ggplot(
    target_exon,
    ggplot2::aes(xstart = start, xend = end, y = transcript_name)
  ) +
    ggtranscript::geom_range(
      fill = strucCol,
      height = strucHeight,
      color = NA
    ) +
    ggtranscript::geom_range(
      data = target_cds,
      fill = strucCol,
      height = 2 * strucHeight,
      color = NA
    ) +
    ggtranscript::geom_intron(
      data = ggtranscript::to_intron(target_exon, "transcript_name"),
      ggplot2::aes(strand = strand),
      color = arrowCol,
      size = 0.3,
      arrow.min.intron.length = 0
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 10, 0, "mm"))

  ############################################################
  geneRion <- ginfo %>% dplyr::filter(type == 'gene')

  # filter gene region
  regionBW <- bigwigFile %>%
    dplyr::filter(seqnames %in% lengthInfo$chr) %>%
    dplyr::filter(start >= (geneRion$start - uped) &
      end <= (geneRion$end + downed))

  ############################################################

  # plot track
  p1 <-
    ggplot2::ggplot(regionBW) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = start,
        xmax = end,
        ymin = 0,
        ymax = score,
        fill = get(sampleAes),
        color = get(sampleAes)
      ),
      show.legend = F
    ) +
    ggplot2::scale_y_continuous(position = "right") +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::ggtitle(unique(lengthInfo$gene)) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.placement = "outside",
      legend.position = "top",
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, size = base_size + 2),
      panel.spacing.y = ggplot2::unit(0, "mm"),
      strip.text = ggplot2::element_text(size = annoTextSize)
    ) +
    ggplot2::ylab(ylab) +
    ggplot2::xlab("") +
    ggplot2::coord_cartesian(expand = 0)

  # whether change trcak colors
  if (is.null(trackCol)) {
    p2 <- p1
  } else {
    p2 <- p1 +
      ggplot2::scale_fill_manual(
        name = "",
        values = trackCol
      ) +
      ggplot2::scale_color_manual(
        name = "",
        values = trackCol
      )
  }

  # whether add facet colors
  if (addfacetCol == FALSE) {
    p3 <- p2 +
      ggh4x::facet_nested_wrap(
        facets = facetVars,
        dir = "v",
        ncol = 1,
        strip.position = "left"
      )
  } else if (addfacetCol == TRUE) {
    strip <- ggh4x::strip_nested(
      background_y =
        ggh4x::elem_list_rect(
          fill = facetFill,
          color = borderCol
        )
    )
    # plot
    p3 <- p2 +
      ggh4x::facet_nested_wrap(
        facets = facetVars,
        dir = "v",
        ncol = 1,
        strip.position = "left",
        strip = strip
      )
  }

  ##########################################
  # combine
  # pcombined <-
    # p3 %>% aplot::insert_bottom(gene_strcture, height = relHeight)
  pcombined <- gene_strcture %>% aplot::insert_top(p3, height = relHeight)

  return(pcombined)
}
