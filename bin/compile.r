
suppressPackageStartupMessages({
    library(argparse)
    library(tidyr)
    library(dplyr)
    library(purrr)
    library(stringr)
    library(readr)
    library(tximeta)
    library(fishpond)
    library(AnnotationDbi)
    library(BiocFileCache)
    library(SummarizedExperiment)
    library(HDF5Array)
    library(DESeq2)
})

samples = 'quant' %>%
        list.files(full.names = TRUE) %>%
        lapply(function(path) {
            quant <- get_quant_files(path, 'long');
            sample <- basename(path);
            name_list <- tibble::lst("files" = quant, "names" = sample);
        }) %>%
        dplyr::bind_rows() %>%
        as.data.frame()

quants <- tximeta(samples, type='oarfish', dropInfReps = TRUE, txOut = TRUE, skipMeta = TRUE)

txnames <- rownames(quants) %>%
    stringr::str_replace("\\|.*", "")

info = '../../references/txinfo.txt'
txinfo <- readLines(info) %>%
        lapply(function(attributes) {
            strsplit(attributes, "; ") %>%
            unlist() %>%
            gsub("^\\s+|\"", "", .) %>%
            lapply(X=c("transcript_id", "gene_id"), FUN=function(attr, attrs) {
                attrs %>%
                    magrittr::extract(grepl(paste0("^",attr,"\\s+"), .)) %>%
                    magrittr::set_names(strsplit(., "\\s+")[[1]][1]) %>%
                    lapply(function(v) stringr::str_split(v, "\\s+", n=2)[[1]][2]) %>%
                    unlist()
            }, attrs=.) %>%
            unlist()
        }) %>%
        do.call(args=., what='rbind') %>%
        as.data.frame() %>%
        distinct()
tx2g <- txinfo %>%
    dplyr::rename(TXNAME = transcript_id, GENEID = gene_id) %>%
    dplyr::select(TXNAME, GENEID)

gquants <- tximeta(samples, type='oarfish', dropInfReps = TRUE, txOut = FALSE, skipMeta = TRUE, tx2gene = tx2g, ignoreAfterBar = TRUE)

HDF5Array::saveHDF5SummarizedExperiment(gquants, "nanopore_gene_quants")
HDF5Array::saveHDF5SummarizedExperiment(quants, "nanopore_transcript_quants")

short_reads <- HDF5Array::loadHDF5SummarizedExperiment('../../novel')

length(intersect(rownames(short_reads), rownames(gquants)))
sample_map <- read.csv('/private/groups/kimlab/data/cambridge/discover/label.csv', header = FALSE)
colnames(sample_map) <- c('illumina_id', 'sample_id', 'nanopore_id')

nanopore_long <- gquants %>%
    SummarizedExperiment::assays() %>%
    magrittr::extract2('counts') %>%
    as.data.frame() %>%
    dplyr::mutate(gene_id = rownames(.)) %>%
    tidyr::pivot_longer(
        cols = -gene_id,
        names_to = 'nanopore_id',
        values_to = 'count'
    ) %>%
    dplyr::left_join(sample_map, by = 'nanopore_id')
illumina_long <- short_reads %>%
    SummarizedExperiment::assays() %>%
    magrittr::extract2('counts') %>%
    as.data.frame() %>%
    dplyr::mutate(gene_id = rownames(.)) %>%
    tidyr::pivot_longer(
        cols = -gene_id,
        names_to = 'illumina_id',
        values_to = 'count'
    ) %>%
    dplyr::mutate(
        illumina_id = stringr::str_replace(illumina_id, "_S.*$", "")
    ) %>%
    dplyr::left_join(sample_map, by = 'illumina_id')

head(illumina_long)
head(nanopore_long)

comparison <- nanopore_long %>%
    dplyr::left_join(illumina_long, by = c('gene_id', 'sample_id', 'illumina_id', 'nanopore_id')) %>%
    dplyr::rename(nanopore_count = count.x, illumina_count = count.y)

linear_model <- lm(nanopore_count ~ illumina_count, data = comparison)

corrplot <- comparison %>%
    dplyr::filter(nanopore_count > 1 & illumina_count > 1) %>%
    dplyr::bind_rows(list('illumina_count' = 1e6, 'nanopore_count' = 1e6)) %>%
    ggplot2::ggplot(ggplot2::aes(x = illumina_count, y = nanopore_count)) +
    ggplot2::geom_density_2d_filled() +
    ggplot2::scale_x_log10(expand = c(0,0)) + ggplot2::scale_y_log10(expand = c(0,0)) +
    ggplot2::geom_abline(slope = coef(linear_model)[2], intercept = coef(linear_model)[1]) +
    ggplot2::labs(title = paste("Pearson Correlation:", round(summary(linear_model)$r.squared, 3)),
         x = "Illumina Counts",
         y = "Nanopore Counts") +
    creater::theme_create()
ggplot2::ggsave("test.png", corrplot)
