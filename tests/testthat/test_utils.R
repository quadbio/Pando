context('Tests for utility functions')

library(testthat)
library(GenomicRanges)
library(matrixStats)
library(doParallel)
library(Pando)



#### find_peaks_near_genes ####

test_peaks <- GRanges(
    ranges = IRanges(
        start = c(10, 20, 30, 40, 45),
        end = c(15, 22, 33, 45, 46)),
    seqnames = 'chr1'
)

test_genes1 <- GRanges(
    ranges = IRanges(
        start = c(1, 25, 30),
        end = c(8, 30, 40)),
    strand = '+',
    seqnames = 'chr1',
    gene_name = LETTERS[1:3]
)

test_genes2 <- GRanges(
    ranges = IRanges(
        start = c(1, 20, 30),
        end = c(8, 30, 40)),
    strand = '-',
    seqnames = 'chr1',
    gene_name = LETTERS[1:3]
)

test_that('find_peaks_near_genes works.', {
    expect_true(
        all(find_peaks_near_genes(test_peaks, test_genes2, upstream=5, downstream=5) == c(1,0,0,0,0,1,1,1,0,0,0,0,1,1,1))
    )
    expect_true(
        all(find_peaks_near_genes(test_peaks, test_genes1, upstream=5, downstream=5) == c(1,0,0,0,0,0,1,1,0,0,0,0,1,1,1))
    )
    expect_true(
        all(find_peaks_near_genes(test_peaks, test_genes1, upstream=5, downstream=5, only_tss=T) ==c(0,0,0,0,0,0,1,1,0,0,0,0,1,0,0))
    )
    expect_true(
        all(find_peaks_near_genes(test_peaks, test_genes2, upstream=5, downstream=5, only_tss=T) == c(1,0,0,0,0,0,0,1,0,0,0,0,0,1,1))
    )
})



#### aggregate_matrix ####

test_mat <- t(matrix(c(
    1,1,1,1,
    2,3,4,5,
    6,7,8,9,
    2,4,6,8
), nrow=4))


test_that('aggregate_matrix works without groups.', {
    expect_true(all(aggregate_matrix(test_mat) == colMeans2(test_mat)))
    expect_true(all(aggregate_matrix(test_mat, fun='sd') == colSds(test_mat)))
})

test_that('aggregate_matrix works with groups.', {
    expect_true(
        all(aggregate_matrix(test_mat, groups=c('A', 'A', 'B', 'C'))[1, ] == colMeans2(test_mat[1:2, ]))
    )
    expect_true(
        all(aggregate_matrix(test_mat, groups=c('A', 'A', 'B', 'A'))[1, ] == colMeans2(test_mat[c(1:2,4), ]))
    )
    expect_true(
        all(aggregate_matrix(test_mat, groups=c('A', 'A', 'B', 'A'), fun='max')[1, ] == colMaxs(test_mat[c(1:2,4), ]))
    )
    expect_true(
        all(aggregate_matrix(test_mat, groups=c('A', 'A', 'B', 'A'), fun='min')[1, ] == colMins(test_mat[c(1:2,4), ]))
    )
    expect_true(
        all(aggregate_matrix(test_mat, groups=c('A', 'A', 'B', 'A'), fun='sd')[1, ] == colSds(test_mat[c(1:2,4), ]))
    )
})

test_that('aggregate_matrix works with function input.', {
    expect_true(
        all(aggregate_matrix(test_mat, groups=c('A', 'A', 'B', 'A'), fun=colMeans2)[1, ] == colMeans2(test_mat[c(1:2,4), ]))
    )
    expect_true(
        all(aggregate_matrix(test_mat, groups=c('A', 'A', 'B', 'A'), fun=colMedians)[1, ] == colMedians(test_mat[c(1:2,4), ]))
    )
})

test_that('aggregate_matrix keeps names.', {
    expect_true(
        all(colnames(aggregate_matrix(test_mat, groups=c('A', 'A', 'B', 'A'), fun=colMeans2)) == colnames(test_mat))
    )
})


#### sparse_cor ####

colnames(test_mat) <- LETTERS[1:ncol(test_mat)]

test_that('sparse_cor gives the same result as cor.', {
    expect_equal(
        as.matrix(sparse_cor(test_mat, method='pearson')), as.matrix(cor(test_mat, method='pearson'))
    )
    expect_equal(
        as.matrix(sparse_cor(test_mat, method='spearman')), as.matrix(cor(test_mat, method='spearman'))
    )
    expect_equal(
        as.matrix(sparse_cor(test_mat, method='kendall')), as.matrix(cor(test_mat, method='kendall'))
    )
})



#### map_par ####

test_list <- list(
    A = c(1,3,4),
    B = c(5,6,8),
    C = c(2,4,5)
)

test_vec <- c(1,2,5,6,3,9)
names(test_vec) <- LETTERS[1:length(test_vec)]

test_that('map_par works on lists.', {
    expect_mapequal(map_par(test_list, max, verbose=F), lapply(test_list, max))
})

test_that('map_par works on vectors.', {
    expect_mapequal(map_par(test_vec, sqrt, verbose=F), lapply(test_vec, sqrt))
})

registerDoParallel(cores=2)

test_that('map_par works on in parallel', {
    expect_mapequal(map_par(test_vec, sqrt, verbose=F, parallel=T), lapply(test_vec, sqrt))
    expect_mapequal(map_par(test_list, sqrt, verbose=F, parallel=T), lapply(test_list, sqrt))
})






