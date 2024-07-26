## Function collection

getstr <-
function(mystring, initial.character, dijige_str_ini = 1,final.character, dijige_str_final = 1 ,plus1 = F)
{

    # check that all 3 inputs are character variables
    if (!is.character(mystring))
    {
        stop('The parent string must be a character variable.')
    }

    if (!is.character(initial.character))
    {
        stop('The initial character must be a character variable.')
    }


    if (!is.character(final.character))
    {
        stop('The final character must be a character variable.')
    }



    # pre-allocate a vector to store the extracted strings
    snippet = rep(0, length(mystring))



    for (i in 1:length(mystring))
    {
        # extract the initial position
        initial.position = gregexpr(initial.character, mystring[i])[[1]][dijige_str_ini] + nchar(initial.character)

        if (plus1 == T) {
            final.position = gregexpr(final.character, mystring[i])[[1]][dijige_str_final] - nchar(final.character) + 1
        }
        else{
            # extract the final position
            final.position = gregexpr(final.character, mystring[i])[[1]][dijige_str_final] - nchar(final.character)
        }

        # extract the substring between the initial and final positions, inclusively
        snippet[i] = substr(mystring[i], initial.position, final.position)
    }

    return(snippet)
}

sc_abs_profile <- function(gr, anno){
    ## center sc annotation
    stop_codon <- anno[anno$type == "stop_codon"]
    stop_codon$exon_number <- as.numeric(stop_codon$exon_number)
    ## keep only one annotation for the splitted stop codon annotation
    stop_codon <- stop_codon %>%
        as.data.frame() %>% group_by(transcript_id) %>%
        filter(exon_number == min(exon_number)) %>%
        makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
    ## set the coordinate of sc to the nucleotide taht most close to 5'
    stop_codon_p <- stop_codon[strand(stop_codon) == "+"]
    end(stop_codon_p) <- start(stop_codon_p)
    stop_codon_n <- stop_codon[strand(stop_codon) == "-"]
    start(stop_codon_n) <- end(stop_codon_n)
    stop_codon <- c(stop_codon_p, stop_codon_n)

    ## find out the position of stop codon on which exon
    exon <- anno[anno$type == "exon"]
    exon <- as.data.frame(exon) %>% group_by(transcript_id) %>%
        mutate(transcript_length = sum(width)) %>%
        makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

    exon$exon_number <- as.numeric(exon$exon_number)
    exon <- exon[exon$transcript_id %in% stop_codon$transcript_id]
    ## calculate the adding length for each exon to the sc
    exon$stop_codon <- start(stop_codon)[match(exon$transcript_id, stop_codon$transcript_id)]
    exon$sc_exon <- ifelse(start(exon) <= exon$stop_codon &
                               end(exon) >= exon$stop_codon, TRUE, FALSE)
    exon$exon_w <- width(exon)

    ## calculate the distance between sc and SS
    exon$sc_dis_start <- 0
    exon$sc_dis_end <- 0
    exon$sc_dis_start[exon$sc_exon == TRUE] <- exon$stop_codon[exon$sc_exon == TRUE] -
        start(exon[exon$sc_exon == TRUE])
    exon$sc_dis_end[exon$sc_exon == TRUE] <- end(exon[exon$sc_exon == TRUE]) -
        exon$stop_codon[exon$sc_exon == TRUE]

    ## find out the sc located exon
    exon <- as.data.frame(exon) %>% group_by(transcript_id) %>%
        mutate(sc_exon_number = exon_number[which(sc_exon)[1]]) %>%
        makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

    ## positive strand
    exon_p <- exon[strand(exon) == "+"]
    ## extract the stop codon exon to prevent them present twice in the data
    exon_p_true <- exon_p[exon_p$sc_exon == TRUE]
    exon_p_true$Add_3 <- exon_p_true$Add_5 <- 0

    ## first calculate the add value before the sc
    exon_p_5 <- exon_p[exon_p$exon_number <= exon_p$sc_exon_number] %>% sort(decreasing = TRUE)
    exon_p_5$exon_w[exon_p_5$sc_exon == TRUE] <- exon_p_5$sc_dis_start[exon_p_5$sc_exon == TRUE]
    exon_p_5 <- as.data.frame(exon_p_5) %>% group_by(transcript_id) %>%
        mutate(Add_5 = c(0, cumsum(exon_w)[-length(exon_w)])) %>%
        makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
    exon_p_5$Add_3 <- 0
    ## then calculate the add value after the sc
    exon_p_3 <- exon_p[exon_p$exon_number >= exon_p$sc_exon_number] %>% sort(decreasing = FALSE)
    exon_p_3$Add_5 <- 0
    exon_p_3$exon_w[exon_p_3$sc_exon == TRUE] <- exon_p_3$sc_dis_end[exon_p_3$sc_exon == TRUE]
    exon_p_3 <- as.data.frame(exon_p_3) %>% group_by(transcript_id) %>%
        mutate(Add_3 = c(0, cumsum(exon_w)[-length(exon_w)])) %>%
        makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

    exon_p <- c(exon_p_5, exon_p_3)

    ## make sure the stop codon containing exon only appear for once
    exon_p <- exon_p[exon_p$sc_exon != TRUE]
    exon_p <- c(exon_p_true, exon_p)

    ## negative strand
    exon_n <- exon[strand(exon) == "-"]
    exon_n_true <- exon_n[exon_n$sc_exon == TRUE]
    exon_n_true$Add_3 <- exon_n_true$Add_5 <- 0
    ## first calculate the add value before the sc
    exon_n_5 <- exon_n[exon_n$exon_number <= exon_n$sc_exon_number] %>% sort(decreasing = FALSE)
    exon_n_5$exon_w[exon_n_5$sc_exon == TRUE] <- exon_n_5$sc_dis_end[exon_n_5$sc_exon == TRUE]
    exon_n_5 <- as.data.frame(exon_n_5) %>% group_by(transcript_id) %>%
        mutate(Add_5 = c(0, cumsum(exon_w)[-length(exon_w)])) %>%
        makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
    exon_n_5$Add_3 <- 0

    ## then calculate the add value after the sc
    exon_n_3 <- exon_n[exon_n$exon_number >= exon_n$sc_exon_number] %>% sort(decreasing = TRUE)
    exon_n_3$Add_5 <- 0
    exon_n_3$exon_w[exon_n_3$sc_exon == TRUE] <- exon_n_3$sc_dis_start[exon_n_3$sc_exon == TRUE]
    exon_n_3 <- as.data.frame(exon_n_3) %>% group_by(transcript_id) %>%
        mutate(Add_3 = c(0, cumsum(exon_w)[-length(exon_w)])) %>%
        makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

    exon_n <- c(exon_n_5, exon_n_3)
    exon_n <- exon_n[exon_n$sc_exon != TRUE]
    exon_n <- c(exon_n_true, exon_n)

    ## merge positive and negative strand
    exon <- c(exon_p, exon_n)

    ## rank  annotation
    exon <- exon[order(exon$level, exon$transcript_support_level, -exon$transcript_length)]
    ## trans_type_rank can be use for the ranking

    ## assign sites to exon annotation
    o <- findOverlaps(gr, exon, select = "first")

    gr$sc_coor <- exon$stop_codon[o]
    gr$exon_number <- exon$exon_number[o]
    gr$sc_exon_number <- exon$sc_exon_number[o]
    gr$exon_start <- start(exon)[o]
    gr$exon_end <- end(exon)[o]
    gr$Add_5 <- exon$Add_5[o]
    gr$Add_3 <- exon$Add_3[o]
    gr$Transcript_ID_stop <- exon$transcript_id[o]

    gr$dist_to_sc <- NA
    ## extract NA guys to prevent calculation error
    gr_NA <- gr[is.na(gr$sc_exon_number)]
    gr <- gr[!is.na(gr$sc_exon_number)]
    ## positive strand
    gr_p <- gr[strand(gr) == "+"]
    ## first calculate the m6A site on the sc_exon
    gr_p$dist_to_sc[gr_p$sc_exon_number == gr_p$exon_number] <-
        start(gr_p)[gr_p$sc_exon_number == gr_p$exon_number] -
        gr_p$sc_coor[gr_p$sc_exon_number == gr_p$exon_number]
    ## m6A site on the exon before sc
    gr_p$dist_to_sc[gr_p$exon_number < gr_p$sc_exon_number] <-
        start(gr_p)[gr_p$exon_number < gr_p$sc_exon_number] -
        gr_p$exon_end[gr_p$exon_number < gr_p$sc_exon_number] -
        gr_p$Add_5[gr_p$exon_number < gr_p$sc_exon_number]
    ## m6A sites on the exon after sc
    gr_p$dist_to_sc[gr_p$exon_number > gr_p$sc_exon_number] <-
        start(gr_p)[gr_p$exon_number > gr_p$sc_exon_number] -
        gr_p$exon_start[gr_p$exon_number > gr_p$sc_exon_number] +
        gr_p$Add_3[gr_p$exon_number > gr_p$sc_exon_number]

    ## negative strand
    gr_n <- gr[strand(gr) == "-"]
    ## first calculate the m6A site on the sc_exon
    gr_n$dist_to_sc[gr_n$exon_number == gr_n$sc_exon_number] <-
        gr_n$sc_coor[gr_n$exon_number == gr_n$sc_exon_number] -
        start(gr_n)[gr_n$exon_number == gr_n$sc_exon_number]
    ## m6A site on the exon before sc
    gr_n$dist_to_sc[gr_n$exon_number < gr_n$sc_exon_number] <-
        gr_n$exon_start[gr_n$exon_number < gr_n$sc_exon_number] -
        start(gr_n)[gr_n$exon_number < gr_n$sc_exon_number] -
        gr_n$Add_5[gr_n$exon_number < gr_n$sc_exon_number]
    ## m6A sites on the exon after sc
    gr_n$dist_to_sc[gr_n$exon_number > gr_n$sc_exon_number] <-
        gr_n$exon_end[gr_n$exon_number > gr_n$sc_exon_number] -
        start(gr_n)[gr_n$exon_number > gr_n$sc_exon_number] +
        gr_n$Add_3[gr_n$exon_number > gr_n$sc_exon_number]

    gr <- c(gr_p, gr_n, gr_NA) %>% sort()
    return(gr)
}

## function to make the sc absolute profile without intron
start_codon_abs_profile <- function(gr, anno){
    ## center sc annotation
    start_codon <- anno[anno$type == "start_codon"]
    start_codon$exon_number <- as.numeric(start_codon$exon_number)
    ## keep only one annotation for the splited start codon annotation
    start_codon <- start_codon %>%
        as.data.frame() %>% group_by(transcript_id) %>%
        filter(exon_number == min(exon_number)) %>%
        makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
    ## set the coordinate of sc to the most 5' nucleotide
    start_codon_p <- start_codon[strand(start_codon) == "+"]
    end(start_codon_p) <- start(start_codon_p)
    start_codon_n <- start_codon[strand(start_codon) == "-"]
    start(start_codon_n) <- end(start_codon_n)
    start_codon <- c(start_codon_p, start_codon_n)

    exon <- anno[anno$type == "exon"]
    exon$exon_number <- as.numeric(exon$exon_number)
    exon <- exon[exon$transcript_id %in% start_codon$transcript_id]
    ## calculate the adding length for each exon to the sc
    exon$start_codon <- start(start_codon)[match(exon$transcript_id, start_codon$transcript_id)]
    exon$start_co_exon <- ifelse(start(exon) <= exon$start_codon &
                                     end(exon) >= exon$start_codon, TRUE, FALSE)
    exon$exon_w <- width(exon)

    ## calculate the distance between sc and SS
    exon$start_co_dis_start <- 0
    exon$start_co_dis_end <- 0
    exon$start_co_dis_start[exon$start_co_exon == TRUE] <-
        exon$start_codon[exon$start_co_exon == TRUE] - start(exon[exon$start_co_exon == TRUE])

    exon$start_co_dis_end[exon$start_co_exon == TRUE] <- end(exon[exon$start_co_exon == TRUE]) -
        exon$start_codon[exon$start_co_exon == TRUE]

    ## find out the sc located exon
    exon <- as.data.frame(exon) %>% group_by(transcript_id) %>%
        mutate(start_co_exon_number = exon_number[which(start_co_exon)[1]]) %>%
        makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

    ## positive strand
    exon_p <- exon[strand(exon) == "+"]
    ## extract the stop codon exon to prevent them present twice in the data
    exon_p_true <- exon_p[exon_p$start_co_exon == TRUE]
    exon_p_true$Add_3 <- exon_p_true$Add_5 <- 0
    ## first calculate the add value before the sc
    exon_p_5 <- exon_p[exon_p$exon_number <= exon_p$start_co_exon_number] %>%
        sort(decreasing = TRUE)
    exon_p_5$exon_w[exon_p_5$start_co_exon == TRUE] <-
        exon_p_5$start_co_dis_start[exon_p_5$start_co_exon == TRUE]
    exon_p_5 <- as.data.frame(exon_p_5) %>% group_by(transcript_id) %>%
        mutate(Add_5 = c(0, cumsum(exon_w)[-length(exon_w)])) %>%
        makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
    exon_p_5$Add_3 <- 0
    ## then calculate the add value after the sc
    exon_p_3 <- exon_p[exon_p$exon_number >= exon_p$start_co_exon_number] %>%
        sort(decreasing = FALSE)
    exon_p_3$Add_5 <- 0
    exon_p_3$exon_w[exon_p_3$start_co_exon == TRUE] <-
        exon_p_3$start_co_dis_end[exon_p_3$start_co_exon == TRUE]
    exon_p_3 <- as.data.frame(exon_p_3) %>% group_by(transcript_id) %>%
        mutate(Add_3 = c(0, cumsum(exon_w)[-length(exon_w)])) %>%
        makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

    exon_p <- c(exon_p_5, exon_p_3)
    exon_p <- exon_p[exon_p$start_co_exon != TRUE]
    exon_p <- c(exon_p_true, exon_p)

    ## negative strand
    exon_n <- exon[strand(exon) == "-"]
    exon_n_true <- exon_n[exon_n$start_co_exon == TRUE]
    exon_n_true$Add_3 <- exon_n_true$Add_5 <- 0
    ## first calculate the add value before the sc
    exon_n_5 <- exon_n[exon_n$exon_number <= exon_n$start_co_exon_number] %>%
        sort(decreasing = FALSE)
    exon_n_5$exon_w[exon_n_5$start_co_exon == TRUE] <-
        exon_n_5$start_co_dis_end[exon_n_5$start_co_exon == TRUE]
    exon_n_5 <- as.data.frame(exon_n_5) %>% group_by(transcript_id) %>%
        mutate(Add_5 = c(0, cumsum(exon_w)[-length(exon_w)])) %>%
        makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
    exon_n_5$Add_3 <- 0

    ## then calculate the add value after the start codon
    exon_n_3 <- exon_n[exon_n$exon_number >= exon_n$start_co_exon_number] %>%
        sort(decreasing = TRUE)
    exon_n_3$Add_5 <- 0
    exon_n_3$exon_w[exon_n_3$start_co_exon == TRUE] <-
        exon_n_3$start_co_dis_start[exon_n_3$start_co_exon == TRUE]
    exon_n_3 <- as.data.frame(exon_n_3) %>% group_by(transcript_id) %>%
        mutate(Add_3 = c(0, cumsum(exon_w)[-length(exon_w)])) %>%
        makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

    exon_n <- c(exon_n_5, exon_n_3)
    exon_n <- exon_n[exon_n$start_co_exon != TRUE]
    exon_n <- c(exon_n_true, exon_n)

    ## merge positive and negative strand
    exon <- c(exon_p, exon_n)

    ## assign the exon base on the rule
    exon <- as.data.frame(exon) %>% group_by(transcript_id) %>%
        mutate(transcript_length = sum(width)) %>%
        mutate(maxExon = max(exon_number)) %>%
        makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

    exon <- exon[order(exon$level, exon$transcript_support_level,
                       -exon$transcript_length)]
    ## assign sites to exon annotation
    o <- findOverlaps(gr, exon, select = "first")

    gr$start_co_coor <- exon$start_codon[o]
    gr$exon_number <- exon$exon_number[o]
    gr$start_co_exon_number <- exon$start_co_exon_number[o]
    gr$exon_start <- start(exon)[o]
    gr$exon_end <- end(exon)[o]
    gr$Add_5 <- exon$Add_5[o]
    gr$Add_3 <- exon$Add_3[o]
    gr$Transcript_ID_start <- exon$transcript_id[o]

    gr$dist_to_start_co <- NA
    ## extract NA guys to prevent calculation error
    gr_NA <- gr[is.na(gr$start_co_exon_number)]
    gr <- gr[!is.na(gr$start_co_exon_number)]
    ## positive strand
    gr_p <- gr[strand(gr) == "+"]
    ## first calculate the m6A site on the start_co_exon
    gr_p$dist_to_start_co[gr_p$start_co_exon_number == gr_p$exon_number] <-
        start(gr_p)[gr_p$start_co_exon_number == gr_p$exon_number] -
        gr_p$start_co_coor[gr_p$start_co_exon_number == gr_p$exon_number]
    ## m6A site on the exon before sc
    gr_p$dist_to_start_co[gr_p$exon_number < gr_p$start_co_exon_number] <-
        start(gr_p)[gr_p$exon_number < gr_p$start_co_exon_number] -
        gr_p$exon_end[gr_p$exon_number < gr_p$start_co_exon_number] -
        gr_p$Add_5[gr_p$exon_number < gr_p$start_co_exon_number]
    ## m6A sites on the exon after sc
    gr_p$dist_to_start_co[gr_p$exon_number > gr_p$start_co_exon_number] <-
        start(gr_p)[gr_p$exon_number > gr_p$start_co_exon_number] -
        gr_p$exon_start[gr_p$exon_number > gr_p$start_co_exon_number] +
        gr_p$Add_3[gr_p$exon_number > gr_p$start_co_exon_number]

    ## negative strand
    gr_n <- gr[strand(gr) == "-"]
    ## first calculate the m6A site on the start_co_exon
    gr_n$dist_to_start_co[gr_n$exon_number == gr_n$start_co_exon_number] <-
        gr_n$start_co_coor[gr_n$exon_number == gr_n$start_co_exon_number] -
        start(gr_n)[gr_n$exon_number == gr_n$start_co_exon_number]
    ## m6A site on the exon before sc
    gr_n$dist_to_start_co[gr_n$exon_number < gr_n$start_co_exon_number] <-
        gr_n$exon_start[gr_n$exon_number < gr_n$start_co_exon_number] -
        start(gr_n)[gr_n$exon_number < gr_n$start_co_exon_number] -
        gr_n$Add_5[gr_n$exon_number < gr_n$start_co_exon_number]
    ## m6A sites on the exon after sc
    gr_n$dist_to_start_co[gr_n$exon_number > gr_n$start_co_exon_number] <-
        gr_n$exon_end[gr_n$exon_number > gr_n$start_co_exon_number] -
        start(gr_n)[gr_n$exon_number > gr_n$start_co_exon_number] +
        gr_n$Add_3[gr_n$exon_number > gr_n$start_co_exon_number]

    gr <- c(gr_p, gr_n, gr_NA) %>% sort()
    gr$dist_to_start_co[which(gr$exon_number > gr$start_co_exon_number)] <-
        gr$dist_to_start_co[which(gr$exon_number > gr$start_co_exon_number)] + 1
    return(gr)
}

readBWtoRLE <- function(x, mainchr = TRUE) {
    rle <- rtracklayer::import.bw(x, as = "Rle")
    if (mainchr == TRUE) {
        rle <- rle[names(rle) %in% standardChromosomes(rle)]
    }
    return(rle)
}

extendposition <- function(gr){
    gr$Position[gr$location == "CDS"]  <- gr$Position[gr$location == "CDS"] + 1
    gr$Position[gr$location == "UTR3"] <- gr$Position[gr$location == "UTR3"] + 2
    return(gr)
}

plot_boruta <- function(boruta_output, title = ""){
    boruta_output_importance <-
        boruta_output$finalDecision %>%
        as.data.frame()
    dd <- data.frame(. = c(rep("Shadow", 3)))
    rownames(dd) <- c("shadowMax", "shadowMean", "shadowMin")
    boruta_output_importance <- rbind(boruta_output_importance, dd)
    df <- as.data.frame(boruta_output$ImpHistory)
    boruta_output <-
        df[,colnames(df) %in% c(rownames(boruta_output_importance))]
    k <- boruta_output
    k[k == -Inf] <- 0
    ord <- apply(k, 2, median)
    ord <- ord[order(ord, decreasing = T)]
    df <- data.frame(value = 0,
                     feature = "a")
    for (i in 1:ncol(boruta_output)) {
        df2 <- data.frame(value = boruta_output[,i],
                          feature = colnames(boruta_output)[i])
        df <- rbind(df, df2)
    }

    df <- df[-1,]
    df$feature <- factor(df$feature, levels = names(ord))
    df$Sig <- boruta_output_importance$.[
        match(df$feature, rownames(boruta_output_importance))
    ]
    df_RPE1 <- df
    dd <- df %>% group_by(df$feature) %>%
        summarise(mean = mean(value))
    dd$seq <- substr(dd$`df$feature`, 1, 2)
    a <- ifelse(dd$seq == "n_", "red", "black")

    df$feature <- as.character(df$feature)
    df$feature[df$feature == "n_m6A_CDS"] <- "CDS m6A"
    df$feature[df$feature == "n_m6A_UTR3"] <- "3'UTR m6A"
    df$feature[df$feature == "n_m6A_UTR5"] <- "5'UTR m6A"
    df$feature <- factor(df$feature,
                         levels = c("CDS m6A", "3'UTR m6A", "5'UTR m6A",
                                    "shadowMax", "shadowMean", "shadowMin"))
    df$Sig <- factor(df$Sig, levels = c("Confirmed", "Shadow", "Tentative", "Rejected"))

    ggplot(df, aes(x = feature, y = value, fill = Sig)) +
        geom_boxplot(outlier.shape = NA) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        ggtitle(title) +
        labs(y = "Importance", x = NULL)
}

norm_cover_sites <- function(gr, rle_p, rle_n, group_name = "", window.size = 20,
                             normalization = T){
    gr <- gr+window.size
    gr_n <- gr[strand(gr) == "-"]
    gr_p <- gr[strand(gr) == "+"]

    ## the truncation signals
    cover_p <- as.matrix(rle_p[gr_p])
    cover_n <- as.matrix(rle_n[gr_n] ) %>% .[,ncol(.):1]
    cover_p[is.na(cover_p)] <- 0
    cover_n[is.na(cover_n)] <- 0

    cover_WT <- rbind(cover_p, cover_n)
    if (normalization == T) {
        tmp <- apply(cover_WT, 1, function(x){
            minx <- min(x)
            maxx <- max(x)
            out <- (x - min(x))/(max(x) - min(x))
            return(out)
        }) %>% t()
        tmp[is.nan(tmp)] <- 0
        df <- data.frame(Signal = colSums(tmp)/length(gr),
                         Position = -window.size:window.size,
                         Group = group_name)
    }
    if (normalization == F) {
        tmp <- apply(cover_WT, 1, function(x){
            minx <- min(x)
            maxx <- max(x)
            out <- (x - min(x))/(max(x) - min(x))
            return(out)
        }) %>% t()
        tmp[is.nan(tmp)] <- 0

        df <- data.frame(Signal = colSums(cover_WT),
                         Position = -window.size:window.size,
                         Group = group_name)
    }
    return(df)
}

cover_sites <- function(gr, rle_p, rle_n, group_name = "", window.size = 20,
         normalization = T, threshold = 4){
    gr <- gr+window.size
    gr_n <- gr[strand(gr) == "-"]
    gr_p <- gr[strand(gr) == "+"]

    ## the truncation signals
    cover_p <- as.matrix(rle_p[gr_p])
    cover_n <- as.matrix(rle_n[gr_n] ) %>% .[,ncol(.):1]
    cover_p[is.na(cover_p)] <- 0
    cover_n[is.na(cover_n)] <- 0

    cover_WT <- rbind(cover_p, cover_n)
    if (normalization == T) {
        df <- data.frame(Signal = colSums(cover_WT >= threshold)/length(gr),
                         Position = -window.size:window.size,
                         Group = group_name)
    }
    if (normalization == F) {
        df <- data.frame(Signal = colSums(cover_WT >= threshold),
                         Position = -window.size:window.size,
                         Group = group_name)
    }
    return(df)
}

cover_sites_v2 <- function(gr, rle_p, rle_n, group_name = "", window.size = 20,
         normalization = T, threshold = 3){
    gr <- gr+window.size
    gr_n <- gr[strand(gr) == "-"]
    gr_p <- gr[strand(gr) == "+"]

    ## the truncation signals
    cover_p <- as.matrix(rle_p[gr_p])
    cover_n <- as.matrix(rle_n[gr_n] ) %>% .[,ncol(.):1]
    cover_p[is.na(cover_p)] <- 0
    cover_n[is.na(cover_n)] <- 0

    cover_WT <- rbind(cover_p, cover_n)
    if (normalization == T) {
        df <- data.frame(Truncation_sites = colSums(cover_WT >= threshold)/length(gr),
                         Position = -window.size:window.size,
                         Group = group_name)
    }
    if (normalization == F) {
        df <- data.frame(Truncation_sites = colSums(cover_WT >= threshold),
                         Position = -window.size:window.size,
                         Group = group_name)
    }
    return(df)
}

last_ss_profile <- function(gr, anno = anno){
    ## PART1 find out the last SS in last CDS exon
    ## make new id to avoid weird sub-transcript, e.g. ENST00000302805.7_PAR_Y
    anno$t_id <- paste0(anno$transcript_id, ",", anno$Parent)
    cds <- anno[anno$type == "CDS"]
    exon <- anno[anno$type == "exon"]
    exon$exon_number <- as.numeric(exon$exon_number)

    ## calculate transcript length for the assignment
    exon <- as.data.frame(exon) %>% group_by(t_id) %>%
        mutate(transcript_length = sum(width)) %>%
        makeGRangesFromDataFrame(., keep.extra.columns = T)

    ## find out the last exon in CDS
    cds$exon_number <- as.numeric(cds$exon_number)
    cds <- as.data.frame(cds) %>% group_by(t_id) %>%
        mutate(max_exon = max(exon_number),
               total_cds = n()) %>%
        makeGRangesFromDataFrame(., keep.extra.columns = T)
    last_cds_exon <- cds[cds$max_exon == cds$exon_number]
    last_cds_exon <- last_cds_exon[last_cds_exon$total_cds != 1]
    last_ss_p <- last_cds_exon[strand(last_cds_exon) == "+"]
    last_ss_n <- last_cds_exon[strand(last_cds_exon) == "-"]
    ## take the 3'SS for the last exon as reference point
    end(last_ss_p) <- start(last_ss_p)
    start(last_ss_n) <- end(last_ss_n)

    last_ss <- c(last_ss_p, last_ss_n)

    ## extract the exon with the annotated CDS and have the intron in CDS
    exon <- exon[exon$t_id %in% last_ss$t_id]
    exon$exon_w <- width(exon)

    exon$last_ss <- start(last_ss)[match(exon$t_id, last_ss$t_id)]
    exon_p <- exon[strand(exon) == "+"]
    exon_n <- exon[strand(exon) == "-"]

    exon_p$last_ss_exon <- ifelse(start(exon_p) == exon_p$last_ss, TRUE, FALSE)
    exon_n$last_ss_exon <- ifelse(end(exon_n)   == exon_n$last_ss, TRUE, FALSE)

    ## find out adding for exons
    ## plus strand
    exon_p_true <- exon_p[exon_p$last_ss_exon == TRUE]
    exon_p_true$Add_3 <- exon_p_true$Add_5 <- 0

    exon_p$last_ss_e_number <- exon_p_true$exon_number[
        match(exon_p$t_id, exon_p_true$t_id)
    ]
    ## first calculate the add value before the last SS
    exon_p_5 <- exon_p[exon_p$exon_number < exon_p$last_ss_e_number] %>% sort(decreasing = TRUE)
    ## assign the distance between stop codon and SS to the stop codon containing exon

    exon_p_5 <- as.data.frame(exon_p_5) %>% group_by(t_id) %>%
        mutate(Add_5 = c(0, cumsum(exon_w)[-length(exon_w)])) %>%
        makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
    exon_p_5$Add_3 <- 0

    ## then calculate the add value after the SS
    exon_p_3 <- exon_p[exon_p$exon_number >= exon_p$last_ss_e_number] %>% sort(decreasing = FALSE)
    exon_p_3$Add_5 <- 0
    exon_p_3 <- as.data.frame(exon_p_3) %>% group_by(t_id) %>%
        mutate(Add_3 = c(0, cumsum(exon_w)[-length(exon_w)])) %>%
        makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

    exon_p <- c(exon_p_5, exon_p_3)

    ## negative strand
    exon_n_true <- exon_n[exon_n$last_ss_exon == TRUE]
    exon_n_true$Add_3 <- exon_n_true$Add_5 <- 0
    exon_n$last_ss_e_number <- exon_n_true$exon_number[
        match(exon_n$t_id, exon_n_true$t_id)
    ]
    ## first calculate the add value before the sc
    exon_n_5 <- exon_n[exon_n$exon_number < exon_n$last_ss_e_number] %>% sort(decreasing = FALSE)
    exon_n_5 <- as.data.frame(exon_n_5) %>% group_by(t_id) %>%
        mutate(Add_5 = c(0, cumsum(exon_w)[-length(exon_w)])) %>%
        makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
    exon_n_5$Add_3 <- 0

    ## then calculate the add value after the sc
    exon_n_3 <- exon_n[exon_n$exon_number >= exon_n$last_ss_e_number] %>% sort(decreasing = TRUE)
    exon_n_3$Add_5 <- 0
    exon_n_3 <- as.data.frame(exon_n_3) %>% group_by(t_id) %>%
        mutate(Add_3 = c(0, cumsum(exon_w)[-length(exon_w)])) %>%
        makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

    exon_n <- c(exon_n_5, exon_n_3)
    exon <- c(exon_p, exon_n)
    exon <- exon[order(exon$level, exon$transcript_support_level, -exon$transcript_length)]

    ## assign sites to exon annotation
    o <- findOverlaps(gr, exon, select = "first")

    gr$exon_number <- exon$exon_number[o]
    gr$last_ss_e_number <- exon$last_ss_e_number[o]
    gr$exon_start <- start(exon)[o]
    gr$exon_end <- end(exon)[o]
    gr$Add_5 <- exon$Add_5[o]
    gr$Add_3 <- exon$Add_3[o]
    gr$Transcript_ID_ss <- exon$transcript_id[o]

    ## extract NA guys to prevent calculation error
    gr_NA <- gr[is.na(gr$last_ss_e_number)]
    gr <- gr[!is.na(gr$last_ss_e_number)]
    ## positive strand
    gr_p <- gr[strand(gr) == "+"]
    ## first calculate the m6A site on the last SS exon
    gr_p$dist_last_ss[gr_p$last_ss_e_number == gr_p$exon_number] <-
        start(gr_p)[gr_p$last_ss_e_number == gr_p$exon_number] -
        gr_p$exon_start[gr_p$last_ss_e_number == gr_p$exon_number]
    ## m6A site on the exon before last SS
    gr_p$dist_last_ss[gr_p$exon_number < gr_p$last_ss_e_number] <-
        start(gr_p)[gr_p$exon_number < gr_p$last_ss_e_number] -
        gr_p$exon_end[gr_p$exon_number < gr_p$last_ss_e_number] -
        gr_p$Add_5[gr_p$exon_number < gr_p$last_ss_e_number]
    ## m6A sites on the exon after last SS
    gr_p$dist_last_ss[gr_p$exon_number >= gr_p$last_ss_e_number] <-
        start(gr_p)[gr_p$exon_number >= gr_p$last_ss_e_number] -
        gr_p$exon_start[gr_p$exon_number >= gr_p$last_ss_e_number] +
        gr_p$Add_3[gr_p$exon_number >= gr_p$last_ss_e_number]

    ## negative strand
    gr_n <- gr[strand(gr) == "-"]
    ## first calculate the m6A site on the last SS exon
    gr_n$dist_last_ss[gr_n$exon_number == gr_n$last_ss_e_number] <-
        gr_n$exon_end[gr_n$exon_number == gr_n$last_ss_e_number] -
        start(gr_n)[gr_n$exon_number == gr_n$last_ss_e_number]
    ## m6A site on the exon before last SS
    gr_n$dist_last_ss[gr_n$exon_number < gr_n$last_ss_e_number] <-
        gr_n$exon_start[gr_n$exon_number < gr_n$last_ss_e_number] -
        start(gr_n)[gr_n$exon_number < gr_n$last_ss_e_number] -
        gr_n$Add_5[gr_n$exon_number < gr_n$last_ss_e_number]
    ## m6A sites on the exon after last SS
    gr_n$dist_last_ss[gr_n$exon_number > gr_n$last_ss_e_number] <-
        gr_n$exon_end[gr_n$exon_number > gr_n$last_ss_e_number] -
        start(gr_n)[gr_n$exon_number > gr_n$last_ss_e_number] +
        gr_n$Add_3[gr_n$exon_number > gr_n$last_ss_e_number]

    gr <- c(gr_p, gr_n, gr_NA) %>% sort()
    return(gr)
}

getres <- function(dds, wt = "WT", treated = "", basemean = 10){
    res <- results(dds, contrast = c("condition", treated, wt))
    res <- as.data.frame(res)
    res <- res[res$baseMean >= basemean,]
    res$gene_id <- substr(rownames(res), 1, 15)
    res$gene_name <- hg38$gene_name[match(rownames(res), hg38$gene_id)]
    return(res)
}

reads.assignment <- function(gr, rle, name = ""){
    reads <- as.matrix(rle[gr])
    mcols(gr)[,ncol(mcols(gr))+1] <- reads
    colnames(mcols(gr))[ncol(mcols(gr))] <- name
    return(gr)
}

read.featureCount <- function(path){
    out <- read.table(path, header = T)
    out <- out[,c(1,6,7)]
    out
}

dis_to_SS <- function(gr, anno){
    gr <- gr[!is.na(gr$Transcript_ID)]
    used_anno <- anno[anno$transcript_id %in% gr$Transcript_ID]
    used_anno <- used_anno[used_anno$type == "exon"]

    o <- findOverlaps(gr, used_anno, select = "first")
    gr$exon_start <- start(used_anno)[o]
    gr$exon_end <- end(used_anno)[o]
    gr$exon_number <- used_anno$exon_number[o]
    gr$exon_max <- used_anno$exon_max[o]
    ## calculate the shortest distance to the splicing junction
    gr$dis_1 <- abs(start(gr) - gr$exon_start)
    gr$dis_2 <- abs(start(gr) - gr$exon_end)


    mcols(gr) <- mcols(gr) %>% as.data.frame %>% rowwise() %>%
        mutate(dis_SS = min(dis_1, dis_2))

    ## assign the distance to 3' and 5' SS
    gr$dis_3SS <- 0
    gr$dis_3SS[as.character(strand(gr)) == "+"] <-
        gr$dis_1[as.character(strand(gr)) == "+"]
    gr$dis_3SS[as.character(strand(gr)) == "-"] <-
        gr$dis_2[as.character(strand(gr)) == "-"]

    gr$dis_5SS <- 0
    gr$dis_5SS[as.character(strand(gr)) == "+"] <-
        gr$dis_2[as.character(strand(gr)) == "+"]
    gr$dis_5SS[as.character(strand(gr)) == "-"] <-
        gr$dis_1[as.character(strand(gr)) == "-"]

    ## for distance to the SS of sites on the first and last exon we need to correct the
    gr_n <- gr[strand(gr) == "-"]
    gr_p <- gr[strand(gr) == "+"]

    #### negative strand
    gr_n$dis_SS[which(gr_n$exon_number == "1")] <-
        gr_n$dis_1[which(gr_n$exon_number == "1")]
    gr_n$dis_SS[which(gr_n$exon_number == gr_n$exon_max)] <-
        gr_n$dis_2[which(gr_n$exon_number == gr_n$exon_max)]

    ##### relative position within exon
    gr_n$exon_relative_position <-
        gr_n$dis_2/(gr_n$exon_end - gr_n$exon_start)

    #### positive strand
    gr_p$dis_SS[which(gr_p$exon_number == "1")] <-
        gr_p$dis_2[which(gr_p$exon_number == "1")]
    gr_p$dis_SS[which(gr_p$exon_number == gr_p$exon_max)] <-
        gr_p$dis_1[which(gr_p$exon_number == gr_p$exon_max)]

    ##### relative position within exon
    gr_p$exon_relative_position <-
        gr_p$dis_1/(gr_p$exon_end - gr_p$exon_start)

    ## merge
    gr <- c(gr_p, gr_n)

    ## assign the max SS distance to the sites which located on the transcript without intron
    max_dis_SS <- max(gr_p$dis_SS, gr_n$dis_SS)

    gr$dis_3SS[gr$exon_number == "1"] <- max_dis_SS
    gr$dis_5SS[gr$exon_number == gr$exon_max] <- max_dis_SS

    gr$dis_SS[gr$exon_max == "1"] <- max_dis_SS
    gr$dis_3SS[gr$exon_max == "1"] <- max_dis_SS
    gr$dis_5SS[gr$exon_max == "1"] <- max_dis_SS
    return(gr)
}

read_Anke_DA <- function(path){
    out <- read.csv(path)
    colnames(out)[7:8] <- c("baseMean", "log2FC")
    return(out)
}


plotMeta <- function(df, title, adjust)
{
    df$Position[df$location == "CDS"] <-
        df$Position[df$location == "CDS"] + 1
    df$Position[df$location == "UTR3"] <-
        df$Position[df$location == "UTR3"] + 2
    p1 <- ggplot(df, aes(x = Position)) +
        geom_density(fill = "white", alpha= 0.05,
                     adjust = adjust, color = "Orange")  +
        ggtitle(title) +
        scale_x_continuous(breaks = c(0.5,1.5,2.5, 5),
                           labels = c("5'UTR","CDS","3'UTR", "No Map")) +
        geom_vline(xintercept=c(1,2), linetype=2, color="cornflowerblue") +
        theme_bw()
    return(p1)
}

make_pie_chart <- function(x, title = ""){
    tmp2 <- table(x) %>% as.data.frame(.)
    colnames(tmp2) <- c("seq", "Number")
    # Compute percentages
    tmp2$fraction = tmp2$Number/sum(tmp2$Number)

    # Compute the cumulative percentages (top of each rectangle)
    tmp2$ymax = cumsum(tmp2$fraction)

    # Compute the bottom of each rectangle
    tmp2$ymin = c(0, head(tmp2$ymax, n=-1))

    # Compute label position
    tmp2$labelPosition <- (tmp2$ymax + tmp2$ymin) / 2

    # Compute a good label
    tmp2$label <- paste0(tmp2$seq, "\n value: ", round(tmp2$fraction,3)*100, "%")

    # Make the plot
    ggplot(tmp2, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=seq)) +
        geom_rect() +
        ggrepel::geom_label_repel(x=3.5, aes(y=labelPosition, label=label), size=6) +
        scale_fill_brewer(palette=4) +
        coord_polar(theta="y") +
        xlim(c(2, 4)) +
        theme_void() +
        theme(legend.position = "none") +
        labs(title = title)
}
