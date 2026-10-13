#' Motif scanning with a corrected background conversion
#'
#' This file vendors the parts of motifmatchr (1.32.0, GPL-3, Alicia Schep)
#' and Signac (1.17.1, MIT, Tim Stuart) that Pando uses to scan candidate
#' regions for motifs (`Signac::AddMotifs()` -> `Signac::CreateMotifMatrix()`
#' -> `motifmatchr::matchMotifs()`), with one bug fixed.
#'
#' The bug: `motifmatchr:::convert_pwm()` should re-express a PWM stored as
#' log2(p / B) (B = `TFBSTools::bg()` of the PWM, 0.25 for Pando's motifs) as
#' log2(p / bg), where bg is the nucleotide composition used for scanning
#' ("subject" by default). It subtracts log2(B) - log2(bg) instead of adding
#' it, so MOODS scans log2(p * bg / B^2) = log2(p * bg / 0.0625). Bases that
#' are common in the background are rewarded instead of penalised. The same
#' sign error affects the natural-log branch, and frequency PWMs are passed
#' through unconverted because `pwmType()` returns 'frequency' while
#' `convert_pwm()` tests for 'prob'.
#'
#' One further deliberate difference from Signac: ranges on contigs missing from
#' the genome get an all-FALSE row. Signac's CreateMotifMatrix() builds that
#' filler with sparseMatrix(i = n_missing, j = n_motifs), which marks the last
#' missing range as matching the last motif.
#'
#' Only the R layer is vendored. The MOODS C++ scanner and its p-value
#' thresholding are unaffected by the bug (they already receive the correct
#' background) and are still called from the installed motifmatchr.
#'
#' @keywords internal
#' @name motifmatchr_patched
NULL


# Compiled MOODS entry points of the installed motifmatchr ---------------------

mm_cpp <- function(name){
    utils::getFromNamespace(name, ns = 'motifmatchr')
}


# PWM conversion (motifmatchr R/utils.R, fixed) --------------------------------

#' Convert a PWM to log-odds against the scanning background
#'
#' Fixed version of `motifmatchr:::convert_pwm()`. Returns
#' log(p / bg_freqs) in the base the PWM is stored in (log2 for the 'log2'
#' and 'frequency' types).
#'
#' @param pwm A \code{PWMatrix}.
#' @param bg_freqs Numeric vector of A, C, G, T background frequencies.
#'
#' @return A 4 x motif length numeric matrix.
#'
#' @keywords internal
convert_pwm_fixed <- function(pwm, bg_freqs){
    type <- motifmatchr::pwmType(pwm)
    mat <- as.matrix(pwm)
    bg_mat <- matrix(bg_freqs, nrow = 4, ncol = ncol(mat), byrow = FALSE)
    pwm_bg <- matrix(TFBSTools::bg(pwm), nrow = 4, ncol = ncol(mat), byrow = FALSE)
    out <- switch(type,
        # Stored log2(p / B): log2(p / bg) = log2(p / B) + log2(B) - log2(bg)
        log2 = mat + log2(pwm_bg) - log2(bg_mat),
        log = mat + log(pwm_bg) - log(bg_mat),
        frequency = log2(mat / bg_mat)
    )
    dimnames(out) <- dimnames(mat)
    return(out)
}

convert_pwms_fixed <- function(pwms, bg_freqs){
    stopifnot(inherits(pwms, 'PWMatrixList'))
    lapply(pwms, convert_pwm_fixed, bg_freqs)
}


# Background frequencies (motifmatchr R/utils.R) -------------------------------

check_bg_fixed <- function(bg_freqs){
    if (length(bg_freqs) != 4){
        stop('Invalid background frequencies -- should be length 4')
    }
    if (!isTRUE(all.equal(sum(bg_freqs), 1)) || min(bg_freqs) <= 0){
        stop('Invalid background frequencies. Should be positive and sum to 1')
    }
    if (!is.null(names(bg_freqs))){
        if (!all(c('A', 'C', 'G', 'T') %in% names(bg_freqs))){
            stop('Background nucleotide frequencies have names that ',
                 "don't match nucleotides! (A,C,G,T)")
        }
        bg_freqs <- bg_freqs[c('A', 'C', 'G', 'T')]
    }
    return(bg_freqs)
}

nuc_freqs <- function(seqs){
    freqs <- Biostrings::letterFrequency(seqs, letters = c('A', 'C', 'G', 'T'))
    if (is.matrix(freqs)) freqs <- colSums(freqs)
    return(freqs / sum(freqs))
}

get_bg_fixed <- function(bg, seqs, genome){
    if (is.numeric(bg)){
        return(check_bg_fixed(bg))
    }
    bg <- match.arg(bg, c('subject', 'genome', 'even'))
    if (bg == 'subject'){
        return(nuc_freqs(seqs))
    }
    if (bg == 'genome'){
        if (!inherits(genome, 'BSgenome')){
            stop("bg = 'genome' requires a BSgenome object")
        }
        # Same as motifmatchr's getNucFreqs() for BSgenome objects
        param <- methods::new('BSParams', X = genome, FUN = Biostrings::letterFrequency)
        freqs <- colSums(do.call(rbind, BSgenome::bsapply(param, letters = c('A', 'C', 'G', 'T'))))
        return(freqs / sum(freqs))
    }
    return(c(A = 0.25, C = 0.25, G = 0.25, T = 0.25))
}


# Scanning (motifmatchr R/match_motifs.R, GRanges/PFMatrixList path) -----------

#' Scan genomic ranges for motifs with the corrected background conversion
#'
#' Replacement for \code{motifmatchr::matchMotifs(pwms, subject, genome,
#' out = 'matches' or 'positions')} for a \code{PFMatrixList} or
#' \code{PWMatrixList} and a \code{GRanges}.
#'
#' @param pwms A \code{PFMatrixList} or \code{PWMatrixList}.
#' @param subject A \code{GRanges} object.
#' @param genome A \code{BSgenome} object.
#' @param bg Background nucleotide frequencies: 'subject' (composition of the
#' scanned ranges, motifmatchr's default), 'genome', 'even' or a numeric
#' vector of A, C, G, T frequencies.
#' @param out 'matches' for a sparse logical ranges x motifs matrix,
#' 'positions' for a \code{GRangesList} of motif hits.
#' @param p.cutoff P-value cutoff for motif matches.
#' @param w Window size for MOODS filtration.
#'
#' @keywords internal
match_motifs_fixed <- function(
    pwms,
    subject,
    genome,
    bg = 'subject',
    out = c('matches', 'positions'),
    p.cutoff = 5e-05,
    w = 7
){
    out <- match.arg(out)
    if (inherits(pwms, 'PFMatrixList')){
        pwms <- do.call(TFBSTools::PWMatrixList, lapply(pwms, TFBSTools::toPWM))
    }
    GenomicRanges::strand(subject) <- '+'
    seqs <- Biostrings::getSeq(genome, subject)
    bg_freqs <- get_bg_fixed(bg, seqs, genome)
    seqs <- as.character(seqs)
    motif_mats <- convert_pwms_fixed(pwms, bg_freqs)

    if (out == 'matches'){
        hits <- mm_cpp('get_motif_ix')(motif_mats, seqs, bg_freqs, p.cutoff, w)
        hits <- methods::as(methods::as(hits, 'lMatrix'), 'CsparseMatrix')
        colnames(hits) <- names(pwms)
        return(hits)
    }

    hits <- mm_cpp('get_motif_positions')(motif_mats, seqs, bg_freqs, p.cutoff, w)
    positions <- lapply(seq_along(motif_mats), function(x){
        m_ix <- which(hits$motif_ix == x - 1)
        GenomicRanges::GRanges(
            GenomicRanges::seqnames(subject)[hits$seq_ix[m_ix] + 1],
            IRanges::IRanges(
                start = GenomicRanges::start(subject[hits$seq_ix[m_ix] + 1]) + hits$pos[m_ix],
                width = ncol(motif_mats[[x]])
            ),
            strand = hits$strand[m_ix],
            score = hits$score[m_ix]
        )
    })
    names(positions) <- names(pwms)
    return(GenomicRanges::GRangesList(positions))
}


# Motif object (Signac AddMotifs.default + CreateMotifMatrix) ------------------

#' Build a Signac \code{Motif} object with the corrected background conversion
#'
#' Drop-in replacement for \code{Signac::AddMotifs()} on a \code{GRanges}
#' that scans with \code{match_motifs_fixed()}.
#'
#' @param object A \code{GRanges} object with the regions to scan.
#' @param genome A \code{BSgenome} object or a string recognised by
#' \code{BSgenome::getBSgenome()}.
#' @param pfm A \code{PFMatrixList}.
#' @param bg Background nucleotide frequencies, see \code{match_motifs_fixed()}.
#' @param p.cutoff P-value cutoff for motif matches.
#' @param verbose Display messages.
#'
#' @return A \code{Motif} object.
#'
#' @keywords internal
add_motifs_fixed <- function(
    object,
    genome,
    pfm,
    bg = 'subject',
    p.cutoff = 5e-05,
    verbose = TRUE
){
    if (is.null(names(pfm))){
        warning("No 'names' attribute found in PFMatrixList. ",
                'Extracting names from individual entries.', immediate. = TRUE)
        names(pfm) <- vapply(pfm, methods::slot, FUN.VALUE = 'character', 'name')
    }
    if (is.character(genome)){
        genome <- BSgenome::getBSgenome(genome = genome)
    }

    # Ranges on contigs missing from the genome get no matches
    keep <- as.character(GenomicRanges::seqnames(object)) %in% GenomeInfoDb::seqlevels(genome)
    if (!all(keep)){
        warning('Not all seqlevels present in supplied genome', immediate. = TRUE)
    }
    scan_ranges <- object[keep]

    log_message('Building motif matrix', verbose = verbose)
    motif_matrix <- match_motifs_fixed(
        pwms = pfm,
        subject = scan_ranges,
        genome = genome,
        bg = bg,
        out = 'matches',
        p.cutoff = p.cutoff
    )
    if (!all(keep)){
        hits <- Matrix::summary(motif_matrix)
        motif_matrix <- Matrix::sparseMatrix(
            i = which(keep)[hits$i], j = hits$j, x = rep(TRUE, nrow(hits)),
            dims = c(length(object), ncol(motif_matrix)),
            dimnames = list(NULL, colnames(motif_matrix))
        )
    }
    rownames(motif_matrix) <- Signac::GRangesToString(object)

    log_message('Finding motif positions', verbose = verbose)
    motif_positions <- match_motifs_fixed(
        pwms = pfm,
        subject = scan_ranges,
        genome = genome,
        bg = bg,
        out = 'positions',
        p.cutoff = p.cutoff
    )
    GenomeInfoDb::seqinfo(motif_positions) <-
        GenomeInfoDb::seqinfo(genome)[GenomeInfoDb::seqlevels(motif_positions)]

    log_message('Creating Motif object', verbose = verbose)
    motif <- Signac::CreateMotifObject(
        data = motif_matrix,
        positions = motif_positions,
        pwm = pfm
    )
    return(motif)
}
