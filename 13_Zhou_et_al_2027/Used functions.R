## Used functions

m6A_assignment_intron_version <- function(gr, anno){
    
    ## processing anno
    exon <- anno[anno$type == "exon"]
    tran_length <- as.data.frame(exon) %>% group_by(transcript_id, gene_id) %>%
        summarise(trans_len = sum(width), .groups = "drop")
    
    exon$trans_len <- tran_length$trans_len[match(exon$transcript_id, tran_length$transcript_id)]
    exon$exon_number <- as.numeric(exon$exon_number)
    
    exon <- as.data.frame(exon) %>% group_by(transcript_id) %>%
        mutate(exon_total = max(exon_number)) %>%
        makeGRangesFromDataFrame(., keep.extra.columns = T)
    exon <- exon[order(exon$level, exon$transcript_support_level, -exon$trans_len)]
    
    gene <- anno[anno$type == "gene"]
    
    ## save the used for later usage
    used <- anno[anno$type %in%c("exon")]
    used$transcript_length <- exon$trans_len[match(used$transcript_id, exon$transcript_id)]
    
    ## get intron annotation
    exon$exonID <- paste0(exon$transcript_id, ",", exon$exon_number)
    introns <- setdiff(gene, exon)
    o <- findOverlaps(introns, gene, select = "first")
    introns$gene_id <- gene$gene_id[o]
    introns$gene_name <- gene$gene_name[o]
    
    ## assgined first to annotations
    used <- used[order(used$level, used$transcript_support_level,
                       -used$transcript_length)]
    introns$type <- "intron"
    introns$transcript_id <- NA
    merge_anno <- c(introns, used)
    
    ## assign the sites based on the hirachy
    o <- findOverlaps(gr, merge_anno, select = "first")
    
    gr$type <- merge_anno$type[o]
    gr$gene_id <- merge_anno$gene_id[o]
    gr$transcript_id <- merge_anno$transcript_id[o]
    gr$gene_name <- merge_anno$gene_name[o]
    return(gr)
}

## function to calculate the distance to SS including introns
dis_to_SS_v2 <- function(gr, anno, width_bs = 5){
    ## define the original width of binding sites
    w <- (width_bs-1)/2
    ## center the bs
    ## Store the original information for the peaks
    gr$oriStart <- start(gr)
    gr$oriEnd <- end(gr)
    
    ## Extract the center of all peaks
    gr$half_length <- width(gr)/2
    ## For the peaks which width = 1, keep them as what they are
    gr$half_length[gr$half_length == 0.5] <- 0
    gr$half_length <- round(gr$half_length)
    GenomicRanges::start(gr) <- GenomicRanges::start(gr) +
        gr$half_length
    GenomicRanges::end(gr) <- start(gr)
    gr$half_length <- NULL
    
    gr$center <- start(gr)
    
    ## identify the 5'ss and 3'ss of intron
    ## calculate the length of transcript for assignment hierarchy
    exon <- anno[anno$type == "exon"]
    tran_length <- as.data.frame(exon) %>% group_by(transcript_id, gene_id) %>%
        summarise(trans_len = sum(width), .groups = "drop")
    
    exon$trans_len <- tran_length$trans_len[match(exon$transcript_id, tran_length$transcript_id)]
    exon$exon_number <- as.numeric(exon$exon_number)
    
    exon <- as.data.frame(exon) %>% group_by(transcript_id) %>%
        mutate(exon_total = max(exon_number)) %>%
        makeGRangesFromDataFrame(., keep.extra.columns = T)
    exon <- exon[order(exon$level, exon$transcript_support_level, -exon$trans_len)]
    
    trans <- anno[anno$type == "transcript"]
    trans$trans_len <- tran_length$trans_len[match(trans$transcript_id, tran_length$transcript_id)]
    trans <- trans[order(trans$level, trans$transcript_support_level, -trans$trans_len)]
    
    ## assign gr to transcript
    o <- findOverlaps(gr, trans, select = "first")
    
    gr$transcript_id <- trans$transcript_id[o]
    
    ## get intron annotation
    exon$exonID <- paste0(exon$transcript_id, ",", exon$exon_number)
    exonsByTranscript <- split(anno[anno$type=="exon"],
                               anno[anno$type=="exon"]$transcript_id)
    trans <- trans[match(names(exonsByTranscript),
                         trans$transcript_id)]
    
    introns <- psetdiff(trans, exonsByTranscript)
    introns <- as(introns, "GRangesList") %>% unlist
    
    introns$transcript_id <- names(introns)
    names(introns) <- NULL
    
    ## make intron number
    introns.n <- introns[strand(introns) == "-"] %>%
        sort(., decreasing = TRUE) %>% as.data.frame() %>%
        group_by(transcript_id) %>%
        mutate(intron_number = seq_len(n())) %>%
        makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
    introns.p <- introns[strand(introns) == "+"] %>%
        sort(.) %>% as.data.frame() %>% group_by(transcript_id) %>%
        mutate(intron_number = seq_len(n())) %>%
        makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
    introns <- c(introns.n, introns.p)
    
    introns$intronID1 <- paste0(introns$transcript_id,
                                ",",introns$intron_number)
    introns$intronID2 <- paste0(introns$transcript_id,
                                ",",introns$intron_number+1)
    ## Assign the length of flanking exons this value will be used in ssprofile
    introns$exon1_l <- width(exon)[match(introns$intronID1,
                                         exon$exonID)]
    introns$exon2_l <- width(exon)[match(introns$intronID2,
                                         exon$exonID)]
    
    introns$transcript_length <-
        width(trans)[match(introns$transcript_id,
                           trans$transcript_id)]
    introns$level <-
        trans$level[match(introns$transcript_id,
                          trans$transcript_id)]
    
    ## tsl means transcript_support_level
    introns$tsl <-
        trans$transcript_support_level[match(introns$transcript_id,
                                             trans$transcript_id)]
    
    ## Remove NA avoid warning message
    introns$level[is.na(introns$level)] <- 6
    introns$tsl[is.na(introns$tsl)] <- 6
    introns$tsl[introns$tsl == "NA"] <-6
    
    ## Change the type of data for the ranking step
    introns$level <- as.integer(introns$level)
    introns$transcript_support_level <- as.integer(introns$tsl)
    introns$tsl <- NULL
    introns <- introns[order(introns$level, introns$transcript_support_level,
                             -introns$transcript_length)]
    
    ## overlap to intron
    o <- findOverlaps(gr, introns, select = "first")
    gr$trans_id_intron <- introns$transcript_id[o]
    ## overlap to exon
    o <- findOverlaps(gr, exon, select = "first")
    gr$trans_id_exon <- exon$transcript_id[o]
    
    ## split the data to intron and exon for the later calculation
    gr_intron <- gr[which(gr$transcript_id == gr$trans_id_intron)]
    gr_exon <- gr[which(gr$transcript_id == gr$trans_id_exon)]
    gr_other <- gr[!gr %over% c(gr_intron, gr_exon)]
    
    ## intron grp
    o <- findOverlaps(gr_intron, introns, select = "first")
    gr_intron$intron_start <- start(introns)[o]
    gr_intron$intron_end <- end(introns)[o]
    
    ## split into positive and negative strand
    gr_intron_p <- gr_intron[strand(gr_intron) == "+"]
    gr_intron_n <- gr_intron[strand(gr_intron) == "-"]
    
    gr_intron_p$dis_5SS <- gr_intron_p$center - gr_intron_p$intron_start
    gr_intron_n$dis_5SS <- gr_intron_n$intron_end - gr_intron_n$center
    
    gr_intron_p$dis_3SS <- gr_intron_p$center - gr_intron_p$intron_end
    gr_intron_n$dis_3SS <- gr_intron_n$intron_start - gr_intron_n$center
    
    gr_intron <- c(gr_intron_p, gr_intron_n)
    
    ## exon group
    o <- findOverlaps(gr_exon, exon, select = "first")
    gr_exon$exon_start <- start(exon)[o]
    gr_exon$exon_end <- end(exon)[o]
    gr_exon$exon_number <- exon$exon_number[o]
    gr_exon$exon_total <- exon$exon_total[o]
    
    ## split into p and n strand
    gr_exon_p <- gr_exon[strand(gr_exon) == "+"]
    gr_exon_n <- gr_exon[strand(gr_exon) == "-"]
    
    ## due to the 1nt difference of boundary of exon and intron here I +/- 1 to adjust the distance
    gr_exon_p$dis_5SS <- gr_exon_p$center - gr_exon_p$exon_end + 1
    gr_exon_n$dis_5SS <- gr_exon_n$exon_start - gr_exon_n$center + 1
    
    gr_exon_p$dis_3SS <- gr_exon_p$center - gr_exon_p$exon_start - 1
    gr_exon_n$dis_3SS <- gr_exon_n$exon_end - gr_exon_n$center - 1
    
    ## make sure it is the distance to intron SS not the exon boundary
    gr_exon <- c(gr_exon_p, gr_exon_n)
    gr_exon$dis_3SS[gr_exon$exon_number == 1] <- NA
    gr_exon$dis_5SS[gr_exon$exon_number == gr_exon$exon_total] <- NA
    
    gr_out <- c(gr_intron, gr_exon, gr_other)
    mcols(gr_out)$dis_SS <- pmin(abs(mcols(gr_out)$dis_3SS),
                                 abs(mcols(gr_out)$dis_5SS),
                                 na.rm = TRUE)
    
    return(gr_out+w)
}

## IRtracker
IRTracker_nanopore <- function(anno, file_bw_path,
                               level = c("1","2"), 
                               tsl = c("1","2","3"),
                               Retainable = TRUE){
    ## The used annotation file should be the gff3 file which loaded with import.gff3 function and 
    ## in a GRange object.
    ## Annotation file filtering
    exon <- anno[anno$type == "exon"]
    transcript <- anno[anno$type == "transcript"]
    
    if (!any(is.na(level)) | !any(is.na(tsl))) {
        if(!any(is.na(level)) & !any(is.na(tsl))){
            exon <- exon[exon$level %in% level & exon$transcript_support_level %in% tsl]
            transcript <- transcript[transcript$level %in% level & 
                                         transcript$transcript_support_level %in% tsl]
        }
        if(!any(is.na(level)) & any(is.na(tsl))){
            exon <- exon[exon$level %in% level]
            transcript <- transcript[transcript$level %in% level]
        }
        if(any(is.na(level)) & !any(is.na(tsl))){
            exon <- exon[exon$transcript_support_level %in% tsl]
            transcript <- transcript[transcript$transcript_support_level %in% tsl]
        }
    } 
    
    ## get intron annotation
    exonsByTranscript <- split(exon,
                               exon$transcript_id)
    transcript <- transcript[match(names(exonsByTranscript),
                                   transcript$transcript_id)]
    
    introns <- psetdiff(transcript, exonsByTranscript)
    introns <- as(introns, "GRangesList") %>% unlist(.)
    
    introns$transcript_id <- names(introns)
    names(introns) <- NULL
    
    ## assign back the annotation information
    introns$level <-
        transcript$level[match(introns$transcript_id,
                               transcript$transcript_id)]
    introns$transcript_support_level <-
        transcript$transcript_support_level[match(introns$transcript_id,
                                                  transcript$transcript_id)]
    introns$gene_type <-
        transcript$gene_type[match(introns$transcript_id,
                                   transcript$transcript_id)]
    introns$gene_id <-
        transcript$gene_id[match(introns$transcript_id,
                                 transcript$transcript_id)] 
    introns$gene_name <-
        transcript$gene_name[match(introns$transcript_id,
                                   transcript$transcript_id)] 
    ## Remove NA avoid warning message
    introns$level[is.na(introns$level)] <- 6
    introns$transcript_support_level[is.na(introns$transcript_support_level)] <- 6
    introns$transcript_support_level[introns$transcript_support_level == "NA"] <-6
    
    ## Define the intron groups by whether it is fully within an exon annotation
    o <- findOverlaps(introns, exon, type = "within")
    introns$group <- "clean"
    introns$group[unique(queryHits(o))] <- "Retainable"
    introns$Intron_ID <- paste0(seqnames(introns), ":",start(introns), "-", 
                                end(introns),"_", strand(introns))
    ## remove duplicate introns and keep the one with lowest level and TSL
    introns <- introns[order(introns$level, introns$transcript_support_level)]
    introns <- unique(introns)
    
    if(Retainable == T){
        introns <- introns[introns$group == "Retainable"]
    }
    introns$Intron_ID <- paste0(seqnames(introns), ":",start(introns), "-", 
                                end(introns),"_", strand(introns))
    ## Fetch the first 20nt of the flanking exons.
    exons_up <- flank(introns, width=20)
    exons_down <- flank(introns, width=20, start=FALSE)
    
    ## load bw files
    files <- list.files(file_bw_path, pattern = ".bw")
    
    df_out <- data.frame(Intron_ID = "No",
                         group = "No",
                         IntronDepth = 0,
                         IRratio = 0,
                         SSratio_L = 0, ## splice sites exon:intron ratio on left
                         SSratio_R = 0,
                         Sample = "No")
    
    ## Fast mean coverage of an RleList over GRanges, order-preserving
    .mean_cov_rlelist_gr <- function(cover_rlelist, gr) {
        stopifnot(is(cover_rlelist, "RleList"), is(gr, "GRanges"))
        
        ## If there are seqlevels in gr not present in cover, drop them (keep NA later)
        keep <- as.character(seqnames(gr)) %in% names(cover_rlelist)
        
        out <- rep(NA_real_, length(gr))
        if (!any(keep)) return(out)
        
        gr2 <- gr[keep]
        gr2$.__idx <- which(keep)
        
        gr_by_seq <- split(gr2, as.character(seqnames(gr2)))
        
        res_list <- mapply(
            function(rle, grs) {
                if (length(grs) == 0L) return(list(idx = integer(0), mean = numeric(0)))
                list(
                    idx  = grs$.__idx,
                    mean = viewMeans(Views(rle, start(grs), end(grs)))
                )
            },
            cover_rlelist[names(gr_by_seq)],
            gr_by_seq,
            SIMPLIFY = FALSE
        )
        
        idx <- unlist(lapply(res_list, `[[`, "idx"), use.names = FALSE)
        mu  <- unlist(lapply(res_list, `[[`, "mean"), use.names = FALSE)
        
        out[idx] <- mu
        out
    }
    
    ## do the calculation through all samples
    for (i in 1:length(files)) {
        cover <- import.bw(paste0(file_bw_path, files[[i]]), as = 'RleList')
        
        ## get the mean signal for the three categories from the current sample.
        introns_signal    <- .mean_cov_rlelist_gr(cover, introns)
        exons_up_signal   <- .mean_cov_rlelist_gr(cover, exons_up)
        exons_down_signal <- .mean_cov_rlelist_gr(cover, exons_down)
        
        ## get boundary ratio, remember to remove the strand information for resize function to make 
        ## sure the calculation result is correct
        tmp <- introns
        strand(tmp) <- "*"
        ss_signal_1 <-  cover[resize(tmp, 1, fix = "start") + 3] %>% as.matrix(.)
        ss_signal_2 <-  cover[resize(tmp, 1, fix = "end") + 3] %>% as.matrix(.)
        
        ss_ratio_1 <- rowSums(ss_signal_1[,1:3])/rowSums(ss_signal_1[,4:6])
        ss_ratio_2 <- rowSums(ss_signal_2[,5:7])/rowSums(ss_signal_2[,2:4])
        
        df <- data.frame(Intron_ID = introns$Intron_ID,
                         group = introns$group,
                         IntronDepth = introns_signal,
                         ## for nanopore since read can easily truncate in intron, I use the max value of exon signal 
                         ## for calculate the IRratio
                         IRratio = introns_signal/pmax(exons_up_signal,exons_down_signal),
                         SSratio_L = ss_ratio_1, ## splice sites exon:intron ratio on left
                         SSratio_R = ss_ratio_2,
                         Sample = files[i])
        df_out <- rbind(df_out, df)
        print(paste0(files[i], " finished"))
    }
    ## remove the first empty row
    df_out <- df_out[-1,]
    
    df_out$gene_id <- introns$gene_id[match(df_out$Intron_ID, introns$Intron_ID)]
    df_out$gene_name <- introns$gene_name[match(df_out$Intron_ID, introns$Intron_ID)]
    
    return(df_out)
}



## IR comparison
make_IR_comparison <- function(df, intronDepth = 1, SSratio = 1.2, AveExpr = 0.1,
                               condition = c("DMSO", "PlaB"), maxExon_expression = 10, 
                               n_rep = 4){
    library(limma)
    ## subset the used data
    df <- df[df$condition %in% condition,]
    df$Sample <- paste0(df$condition, "_", df$rep)
    ## filtering
    df_filtered <- df[df$IntronDepth > intronDepth & df$IRratio <= 1 & 
                          (df$SSratio_L > SSratio | df$SSratio_R > SSratio) &
                          df$maxExon_expression >= maxExon_expression,]
    tmp <- table(df_filtered$Intron_ID)
    
    ## keep the ID which happen in at least xx samples
    df_filtered_id <- names(tmp[tmp >= n_rep])
    ## make df for analysis
    df_new <- df_filtered[df_filtered$Intron_ID %in% df_filtered_id,] %>%
        dplyr::select(Intron_ID, Sample, IRratio) %>%
        pivot_wider(names_from = Sample, values_from = IRratio) %>%
        column_to_rownames("Intron_ID") %>%
        as.matrix()
    
    ## adjust this for different data
    sample_info <- data.frame(
        Sample = colnames(df_new),
        condition = sub("_.*", "", colnames(df_new))
    )
    
    sample_info$condition <- factor(sample_info$condition, levels = condition)
    
    ## Design
    design <- model.matrix(~ 0 + condition, data = sample_info)
    colnames(design) <- levels(sample_info$condition)
    
    contrast_name <- paste0(condition[2], "_vs_", condition[1])
    contrast_expr <- paste(condition[2], "-", condition[1])
    
    contrast_mat <- makeContrasts(
        contrasts = contrast_expr,
        levels = design
    )
    
    colnames(contrast_mat) <- contrast_name
    
    fit <- lmFit(df_new, design)
    fit <- contrasts.fit(fit, contrast_mat)
    fit <- eBayes(fit)
    
    res <- topTable(fit,
                    coef = contrast_name,
                    number = Inf,
                    sort.by = "P"
    ) %>%
        dplyr::rename(delta_IR_lm = logFC)
    
    out <- merge(df_new, res, by=0, all=TRUE) %>%
        dplyr::rename(Intron_ID = Row.names)
    out <- out[out$AveExpr >= AveExpr,]
    ## add mean IRratio for each condition
    out[[paste0(condition[[1]], "_mean")]] <- 
        rowMeans(out[, grep(condition[[1]], colnames(out))])
    out[[paste0(condition[[2]], "_mean")]] <- 
        rowMeans(out[, grep(condition[[2]], colnames(out))])
    ## LFC calculation
    out$lgo2FC_IRratio <- 
        log2(out[[paste0(condition[[2]], "_mean")]]/out[[paste0(condition[[1]], "_mean")]])
    out$maxExon_expression <- df$maxExon_expression[match(
        out$Intron_ID, df$Intron_ID
    )]
    return(out)
}