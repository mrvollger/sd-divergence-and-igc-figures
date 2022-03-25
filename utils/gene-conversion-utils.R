#!/usr/bin/env Rscript
if (!require("karyoploteR")) BiocManager::install("karyoploteR")
if (!require("GenomicRanges")) BiocManager::install("GenomicRanges")
if (!require("argparse")) BiocManager::install("argparse")
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("ggnewscale")) install.packages("ggnewscale")
if (!require("ggrepel")) install.packages("ggrepel")
if (!require("data.table")) install.packages("data.table")
if (!require("glue")) install.packages("glue")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
if (!require("scales")) install.packages("scales")
if (!require("cowplot")) install.packages("cowplot")
if (!require("argparse")) install.packages("argparse")
library(ggpubr)
library(ggforce)
library(IRanges)
library(openxlsx)
library(data.table)
library(circlize)
library(glue)
library(ggExtra)
library(HelloRanges)
odir <<- "figures"

my_ggsave <- function(filename, plot = last_plot(), ...) {
    filename <- glue(filename)
    ggsave(filename, plot = plot, ...)
    basename <- tools::file_path_sans_ext(filename)
    table_out <- paste0(basename, ".datatable.tsv")
    data <- apply(plot$data, 2, as.character)
    print(head(data))
    write.table(
        data,
        file = table_out,
        sep = "\t",
        row.names = F,
        quote = F
    )
}

acceptor_columns <- c("#reference_name.liftover", "reference_start.liftover", "reference_end.liftover")
donor_columns <- c("reference_name", "reference_start", "reference_end")

convert_to_sd_coords <- function(in_df, cols, converter) {
    # cols = c(2,5,6); in_df=merged_df
    conv <- converter$converter %>%
        arrange(seqnames, start, end) %>%
        data.table()

    start_col <- cols[2]
    end_col <- cols[3]
    col_names <- colnames(in_df)
    o <- findOverlaps(
        toGRanges(in_df[, ..cols]),
        toGRanges(conv)
    )
    print(length(unique(queryHits(o))))
    print(length(queryHits(o)))

    colnames(conv)[1:4] <- c("conv-nm", "conv-st", "conv-en", "conv-width")

    o_df <- cbind(
        in_df[queryHits(o)],
        conv[subjectHits(o)]
    )
    max_en <- o_df$new_start + o_df$`conv-width`
    min_st <- o_df$new_start

    new_start <- o_df[, ..start_col][[1]] - o_df$subtract
    new_start <- pmax(new_start, min_st)
    new_start <- pmin(new_start, max_en)
    o_df[[eval(start_col)]] = new_start

    new_end <- (round(o_df[, ..end_col]) - round(o_df$subtract))[[1]]
    new_end <- pmax(new_end, min_st)
    new_end <- pmin(new_end, max_en)
    o_df[[eval(end_col)]] = new_end

    as.data.table(o_df[, ..col_names])
}


make_sd_coord_system <- function(sds,
                                 spacer = 100e3,
                                 min_size = 50e3,
                                 second_filter = 1e5,
                                 cens = CYTOFULL[CYTOFULL$gieStain == "acen"]) {
    ranges <- GenomicRanges::reduce(toGRanges(sds))

    for (i in seq(1)) {
        r_gaps <- GenomicRanges::gaps(ranges)
        n_gaps <- length(r_gaps)
        gap_is_smaller_than_joining_sds <- width(r_gaps) / 2 < width(ranges[1:n_gaps]) | width(r_gaps) / 2 < width(ranges[2:(n_gaps + 1)])
        ranges <- GenomicRanges::reduce(c(ranges, r_gaps[gap_is_smaller_than_joining_sds]))
        ranges <- ranges[width(ranges) >= min_size]
    }
    ranges <- ranges[width(ranges) >= second_filter]
    ranges <- GenomicRanges::reduce(ranges + spacer / 2)
    start(ranges[start(ranges) < 0]) <- 0
    sum(width(ranges)) / 1e6
    sum(width(GenomicRanges::reduce(ranges))) / 1e6
    ranges

    tmp <- setdiff(cens, ranges) %>%
        as.data.table() %>%
        group_by(seqnames) %>%
        filter(width == max(width) & width > 1e6) %>%
        mutate(end = start + 1e6) %>%
        dplyr::select(-width, -strand) %>%
        data.table()
    dcens <- setdiff(toGRanges(tmp), ranges)
    dcens$gieStain <- "acen"
    width(dcens)
    ranges <- sort(c(ranges, dcens))
    sum(width(ranges))
    sum(width(GenomicRanges::reduce(ranges)))
    is.unsorted(ranges)
    # intersect(cens,ranges)

    start <- 0
    chrm <- ""
    starts <- NULL
    for (idx in 1:length(ranges)) {
        range <- ranges[idx]
        cur_chrm <- as.character(seqnames(range))
        if (cur_chrm != chrm) {
            chrm <- cur_chrm
            start <- 0
        }
        starts <- c(starts, start)
        start <- start + width(range)
    }
    length(starts)
    length(ranges)
    odf <- as.data.table(ranges)
    odf$new_start <- starts
    odf$subtract <- odf$start - odf$new_start
    cyto <- odf %>%
        mutate(
            name = paste(round(start / 1e6, 1), round(end / 1e6, 1), sep = "-"),
            start = start - subtract,
            end = end - subtract,
            n = 1:nrow(odf),
            gieStain = case_when(gieStain == "acen" ~ "acen", n %% 2 == 0 ~ "gneg", TRUE ~ "gpos50")
        ) %>%
        dplyr::select(-width, -strand, -new_start, -subtract, -n) %>%
        data.frame()
    genome <- cyto %>%
        group_by(seqnames) %>%
        summarise(start = min(start), end = max(end)) %>%
        data.frame()
    print(sum(genome$end - genome$start) / 1e6)
    print(nrow(cyto))
    list(genome = genome, cyto = cyto, converter = odf)
}


make_gc_ideogram <- function(a, d, chromosomes, outfile, height = 36, width = 45) {
    samples <- unique(a$label)
    n_samples <- length(samples)
    samples
    s <- 1.5
    pdf(outfile, height = height, width = width)
    pp <- getDefaultPlotParams(plot.type = 3)
    pp$leftmargin <- pp$leftmargin * 2
    pp$rightmargin <- pp$rightmargin / 5
    pp$data1height <- 200
    pp$data2height <- 200
    pp$topmargin <- pp$topmargin / n_samples * 20
    pp$bottommargin <- 0
    pp$ideogramlateralmargin <- pp$ideogramlateralmargin * 3
    pp$ideogramheight <- pp$ideogramheight
    pp$data1inmargin <- pp$ideogramheight
    pp$data2inmargin <- 0 # pp$ideogramheight / 4
    offset <- -(pp$ideogramheight + pp$data1inmargin + pp$data2inmargin) / pp$data1height
    kp <- plotKaryotype(
        cytobands = converter$cyto,
        genome = converter$genome, plot.type = 3,
        plot.params = pp, chromosomes = chromosomes
    )
    # kpAddBaseNumbers(kp,minor.tick.dist = 1e6, minor.tick.len = 20, minor.tick.col = "black",tick.dist = 5e6, tick.len = 40, tick.col="red", cex=0.5)
    i <- 1
    for (sample in rev(sort(samples))) {
        print(paste(sample, i * 100 / n_samples))
        h <- 1.0 / n_samples

        if (i %% 2 == 0 & F) {
            kpDataBackground(kp, r0 = h * i - h, r1 = h * i, color = "#ededed")
        } else if (F) {
            kpDataBackground(kp, r0 = h * i - h, r1 = h * i, color = "#e1e1e1")
        }
        ## add label
        kpText(kp,
            chr = chromosomes[[1]],
            x = 0, y = 0.5, labels = sample,
            r0 = h * i - h, r1 = h * i, pos = 2
        )
        kpLines(kp,
            chr = converter$genome$seqnames,
            x = c(converter$genome$start, converter$genome$end),
            y = c(0.5, 0.5), r0 = h * i - h, r1 = h * i, col = "gray"
        )

        is_sample <- a$label == sample
        gra <- toGRanges(a[acceptor_left & is_sample])
        grd <- toGRanges(d[acceptor_left & is_sample])
        gra2 <- toGRanges(a[(!acceptor_left) & is_sample])
        grd2 <- toGRanges(d[(!acceptor_left) & is_sample])

        if (length(gra) > 0 && length(grd) > 0) {
            kpPlotLinks(kp,
                data = gra, data2 = grd, col = acceptor_color,
                r0 = h * i - h / 2, r1 = h * i,
                border = transparent(acceptor_color, 0.75)
            )
        }

        if (length(gra2) > 0 && length(grd2) > 0) {
            kpPlotLinks(kp,
                data = gra2, data2 = grd2, col = donor_color,
                r0 = offset - h * i + h / 2, r1 = offset - h * i,
                data.panel = 2,
                border = transparent(donor_color, 0.75)
            )
        }
        i <- i + 1
    }

    label_rows <- 1
    label.df <- converter$cyto %>%
        group_by(seqnames) %>%
        filter(seqnames %in% chromosomes) %>%
        arrange(-(end - start)) %>%
        mutate(y = (1:n() %% 2) * 0.05) %>%
        head(width * 2) %>%
        data.table()
    for (i in seq(label_rows)) {
        r0 <- (i - 1) / label_rows
        r1 <- (i - 1) / label_rows
        if (nrow(label.df) > 0) {
            kpPlotMarkers(kp,
                data = toGRanges(label.df),
                labels = label.df$name,
                y = 0.05 + label.df$y,
                data.panel = 2,
                max.iter = 1000,
            )
        }
    }

    kpPlotRegions(kp, data = toGRanges(GENES_V1.1.sd), r0 = -pp$data1inmargin / pp$data1height * 0.5, r1 = 0, avoid.overlapping = FALSE, border = NA)
    kpPlotRegions(kp,
        data = toGRanges(morb_sd),
        r0 = -pp$data1inmargin / pp$data1height, r1 = -pp$data1inmargin / pp$data1height * 0.5,
        col = transparent(morb_sd$color, 0.25), avoid.overlapping = TRUE, border = NA
    )

    dev.off()
    rm(kp)
    gc()
}


add_gene_list_col <- function(
    df, 
    cols=c(1,2,3), 
    merge=TRUE, 
    genelist=GENES,
    gene_cols=c(1,2,3), 
    outcol="genes",
    gene_name_col="gene"
) {
    dim(df)
    o = GenomicRanges::findOverlaps(toGRanges(df[,..cols]), toGRanges(genelist[,..gene_cols]))
    group_cols = colnames(df)
    odf = cbind(df[queryHits(o)], genelist[subjectHits(o), ..gene_name_col])
    missing = setdiff(seq(nrow(df)), queryHits(o))

    odf %>% 
        group_by_at(vars(group_cols)) %>%
        summarise(crazy=paste(unique(!!as.symbol(gene_name_col)), collapse = ";")) %>%
        bind_rows(., df[missing,]) %>%
        rename(!!outcol := crazy) %>%
        data.table()
}
