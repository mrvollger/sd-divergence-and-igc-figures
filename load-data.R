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
library(gt)
# library(multidplyr)
library("tidylog", warn.conflicts = FALSE)
# cluster <- new_cluster(16)
# cluster_library(cluster, "dplyr")
odir <<- "figures"
source("utils/snv-setup.R")

if (F) {
    acceptor_color <- "darkblue"
    donor_color <- "darkorange"
    igc_colors <- list(Acceptor = acceptor_color, Donor = donor_color)
    acceptor_columns <- c(
        "#reference_name.liftover",
        "reference_start.liftover",
        "reference_end.liftover"
    )
    donor_columns <- c(
        "reference_name",
        "reference_start",
        "reference_end"
    )
    load("Rdata/plotutils.data")
    source("utils/gene-conversion-utils.R")
    rm(ALL_ALN)
    rm(ALL_GENES)
    rm(ALLPAV)
    rm(RM)
    rm(GENES)
    rm(SEDEF)
    rm(SEDEF_38)
    rm(SEDEF_CELARA)
    rm(tmpIndelDEL)
    rm(tmpIndelINS)
    rm(tmpINS)
    rm(tmpINV)
    rm(tmpSNVs)
    gc()

    #
    # read in the genes
    #
    pli <- fread("data/gnomad.v2.1.1.lof_metrics.by_gene.txt")
    pli_cols <- c(
        "gene_id", "pLI", "exac_pLI",
        "cds_length", "num_coding_exons"
    )
    gene_bed <- fread("data/CHM13.combined.v4.bb.bed") %>%
        mutate(gene = gsub("-\\d+$", "", name)) %>%
        rename(
            chr = `#chrom`,
            start = chromStart,
            end = chromEnd
        ) %>%
        mutate(
            gene_id = gsub(".[0-9]+$", "", sourceGene),
            gene_name = gsub("CHM13_", "", gene)
        ) %>%
        merge(
            pli[, ..pli_cols],
            by = "gene_id", all.x = TRUE
        ) %>%
        data.table()

    gene_bed6 <- gene_bed %>%
        separate_rows(blockSizes,
            chromStarts, exonFrames,
            convert = TRUE, sep = ","
        ) %>%
        group_by(sourceTranscript, name) %>%
        mutate(
            end = start + chromStarts + blockSizes,
            start = start + chromStarts,
            exonNum = seq(n()),
            gene = paste(gene, exonNum, sep = "__")
        ) %>%
        data.table()
    dim(gene_bed)
    dim(gene_bed6)
    GENE_BED <- gene_bed
    GENE_BED6 <- gene_bed6


    # these two files are the same but the second has a "group" column
    # f <- "./data/gene-conversion/acceptor_gene_conversion_windows.bed.gz"
    f <- "data/gene-conversion/merged_acceptor.bed.gz"
    df <- fread(f, stringsAsFactors = T)
    dim(df)

    link <- "https://eichlerlab.gs.washington.edu/help/mvollger/share/20130606_g1k_3202_samples_ped_population.txt"
    pop <- fread(link, stringsAsFactors = T) %>%
        dplyr::select(SampleID, Sex, Population, Superpopulation) %>%
        data.table()
    pop_small <- pop %>%
        dplyr::select(SampleID, Superpopulation) %>%
        data.table()
    df <- df %>%
        filter(status != "Donor" & sample != "CHM1_2") %>%
        filter(sample != "GRCh38_2") %>%
        filter(reference_name != "chrY") %>%
        filter(`#reference_name.liftover` != "chrY") %>%
        group_by(sample) %>%
        mutate(Sample = gsub("_\\d", "", sample)) %>%
        arrange(-name) %>%
        mutate(
            length = as.numeric(
                reference_end.liftover - reference_start.liftover
            )
        ) %>%
        mutate(cum.bp = cumsum(length)) %>%
        merge(pop, by.x = "Sample", by.y = "SampleID", all.x = TRUE) %>%
        data.table()

    df[is.na(Superpopulation)]$Superpopulation <- ""
    df[Superpopulation == ""]$Superpopulation <- "EUR"

    #
    # Add gene information to the IGC table
    #
    df <- df %>%
        # ADDS THE ACCEPTOR GENES
        add_gene_list_col(
            cols = acceptor_columns, gene_cols = c(2, 3, 4),
            genelist = gene_bed,
            outcol = "acceptor_gene_names", gene_name_col = "gene_name"
        ) %>%
        add_gene_list_col(
            cols = acceptor_columns, gene_cols = c(2, 3, 4),
            genelist = gene_bed,
            outcol = "acceptor_gene_ids", gene_name_col = "gene_id"
        ) %>%
        # ADDS THE DONOR GENES
        add_gene_list_col(
            cols = donor_columns, gene_cols = c(2, 3, 4),
            genelist = gene_bed,
            outcol = "donor_gene_names", gene_name_col = "gene_name"
        ) %>%
        add_gene_list_col(
            cols = donor_columns, gene_cols = c(2, 3, 4),
            genelist = gene_bed,
            outcol = "donor_gene_ids", gene_name_col = "gene_id"
        ) %>%
        # ADDS THE ACCEPTOR EXONS
        add_gene_list_col(
            cols = acceptor_columns, gene_cols = c(2, 3, 4),
            genelist = gene_bed6,
            outcol = "acceptor_exon_names", gene_name_col = "gene"
        ) %>%
        add_gene_list_col(
            cols = acceptor_columns, gene_cols = c(2, 3, 4),
            genelist = gene_bed6,
            outcol = "acceptor_exon_ids", gene_name_col = "gene_id"
        ) %>%
        # ADDS THE DONOR GENES
        add_gene_list_col(
            cols = donor_columns, gene_cols = c(2, 3, 4),
            genelist = gene_bed6,
            outcol = "donor_exon_names", gene_name_col = "gene"
        ) %>%
        add_gene_list_col(
            cols = donor_columns, gene_cols = c(2, 3, 4),
            genelist = gene_bed6,
            outcol = "donor_exon_ids", gene_name_col = "gene_id"
        ) %>%
        data.table()
    #
    # make copies of the dataframe so that if I use "df" else where I will still have these
    #
    igc <- copy(df)
    IGC <- copy(df)
    pre_merged_df <- copy(df)

    #
    # make a copy of total for super pop ecdfs
    #
    dfc <- copy(df)
    dfc$Superpopulation <- "Total"
    df.all <- rbind(df, dfc)
    #
    # Read in merged groups
    #
    # merged_df <-
    length(unique(
        paste0(
            pre_merged_df$group,
            pre_merged_df$reference_name,
            pre_merged_df$`#reference_name.liftover`
        )
    ))
    merged_df <- pre_merged_df %>%
        group_by(group, `#reference_name.liftover`, `reference_name`) %>%
        mutate(
            n = n(),
            igc_count = n(),
            sample_hap = as.character(sample),
            sample = paste(sample[1], n, sep = " X "),
        ) %>%
        group_by(
            group, `#reference_name.liftover`, `reference_name`
        ) %>%
        summarise(
            across(
                everything(),
                ~ if (is.numeric(.)) mean(., na.rm = TRUE) else if (is.character(.)) paste(unique(unlist(strsplit(., ";"))), collapse = ";")
            )
        ) %>%
        data.table()
    dim(merged_df)

    # change coordiantes
    source("utils/gene-conversion-utils.R")
    cols <- c(acceptor_columns, donor_columns, "sample")
    igc_a <- igc[, ..acceptor_columns]
    igc_d <- igc[, ..donor_columns]

    # make converted space
    SDs_to_use <- as.data.table(
        GenomicRanges::reduce(
            toGRanges(rbind(SEDEF_V1.1, SEDEF_V1.1_LOWID))
        )
    )
    # GC_SDs_bool = overlaps(SDs_to_use, igc_a) | overlaps(SDs_to_use, igc_d)
    bedlength(SDs_to_use) / 1e6
    converter <- make_sd_coord_system(sds = SDs_to_use)

    tmp <- convert_to_sd_coords(igc, acceptor_columns, converter)
    dim(tmp)
    igc_sd <- convert_to_sd_coords(
        tmp,
        donor_columns,
        converter
    )
    dim(igc)
    dim(igc_sd)
    mean(igc_sd$reference_end - igc_sd$reference_start) / 1e6
    mean(igc_sd$reference_end.liftover - igc_sd$reference_start.liftover) / 1e6

    ALL_GENES_V1.1.sd <- convert_to_sd_coords(
        ALL_GENES_V1.1, c(1, 2, 3), converter
    )
    GENES_V1.1.sd <- convert_to_sd_coords(GENES_V1.1, c(1, 2, 3), converter)

    f <- "https://eichlerlab.gs.washington.edu/help/mvollger/share/morbity_map_chm13_v1.1.tsv"
    morb <- fread(f) %>%
        mutate(color = case_when(
            DEL && DUP ~ "purple",
            DUP ~ "blue",
            DEL ~ "red",
            TRUE ~ "darkgreen"
        )) %>%
        data.table()
    morb_sd <- convert_to_sd_coords(morb, c(1, 2, 3), converter)

    # drop -> c(
    #    "anno_Sat", "anno_Exons", "anno_InterGenicSD",
    #    "anno_RM", "dist_TSS",
    #    "dist_Exons", "dist_InterGenicSD", "dist_InterGenicUnique"
    # )
    HPRC_SNVS <- fread("data/sd-divergence-results/all_snv_exploded.bed.gz",
        showProgress = TRUE,
        nThread = 8,
        stringsAsFactors = TRUE,
        logical01 = TRUE # , drop = drop
    ) %>%
        merge(pop_small, by.x = "SAMPLE", by.y = "SampleID") %>%
        data.table()
    HPRC_SNVS[is.na(Superpopulation)]$Superpopulation <- ""
    HPRC_SNVS[Superpopulation == ""]$Superpopulation <- "EUR"

    colnames(HPRC_SNVS)
    t_grouped_hprc_snvs <- HPRC_SNVS %>%
        filter((anno_SD == T | anno_Unique == T) & anno_TRF == F) %>%
        dplyr::select(
            SAMPLE, anno_SD, anno_Unique,
            anno_IGC, POS, HAP, `#CHROM`,
            Superpopulation
        ) %>%
        group_by(`#CHROM`, `SAMPLE`, `HAP`, `anno_SD`, `anno_Unique`) %>%
        arrange(`#CHROM`, `POS`) %>%
        mutate(dist_to_next_snv = POS - lag(POS)) %>%
        data.table()
    # rm(HPRC_SNVS)

    t_grouped_hprc_snvs %>%
        drop_na() %>%
        group_by(Superpopulation, anno_SD) %>%
        summarise(
            mean = mean(dist_to_next_snv),
            median = median(dist_to_next_snv)
        ) %>%
        pivot_longer(cols = c("mean", "median"))

    # gc()
    grouped_hprc_snvs <- t_grouped_hprc_snvs %>%
        drop_na() %>%
        filter(dist_to_next_snv < 1e5) %>%
        mutate(
            color = case_when(anno_SD == T ~ "SD", T ~ "Unique"),
            AFR = case_when(Superpopulation == "AFR" ~ "AFR", T ~ "Non-AFR")
        ) %>%
        ungroup() %>%
        data.table()

    rm(t_grouped_hprc_snvs)
    gc()


    ####################################################################################################
    ####################################################################################################
    ####################################################################################################

    source("utils/snv-setup.R")
    infile <- "data/sd-divergence-results/long_windows_with_snv_dist_annotation.bed.gz"
    snv_windows <- read_in_snv_windows(infile)
    chrX <- copy(snv_windows[snv_windows[["#chr"]] == "chrX"])
    chrX$region <- "chrX"
    snv_windows <- rbind(snv_windows, chrX)
    snv_windows$region <- factor(snv_windows$region, levels = names(COLORS))

    ####################################################################################################
    ####################################################################################################
    ####################################################################################################

    con <- pipe("pigz -p32 > Rdata/igc_pigz.Rdata", "wb")
    # same as save.image, but save.image does not work with pipes
    save(list = ls(all.names = TRUE), file = con, envir = .GlobalEnv)
    close(con)
} else {
    load("Rdata/igc_pigz.Rdata", verbose = TRUE)
}
gc()
source("utils/gene-conversion-utils.R")