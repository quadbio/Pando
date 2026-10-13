context('Tests for the patched motifmatchr background conversion')

library(testthat)
library(TFBSTools)
library(Pando)


pfm <- PFMatrix(profileMatrix = matrix(
    c(8, 1, 1, 0,  1, 8, 0, 1,  0, 1, 8, 1), nrow = 4,
    dimnames = list(c('A', 'C', 'G', 'T'))
))
uniform <- c(A = 0.25, C = 0.25, G = 0.25, T = 0.25)
bg <- c(A = 0.2, C = 0.3, G = 0.3, T = 0.2)


test_that('log2 PWMs are converted to log2(p / bg)', {
    for (pc in c(0, 0.8)){
        pwm <- toPWM(pfm, type = 'log2probratio', bg = uniform, pseudocounts = pc)
        p <- 2^as.matrix(pwm) * 0.25
        expect_equal(unname(Pando:::convert_pwm_fixed(pwm, bg)), unname(log2(p / bg)))
    }
})

test_that('natural log and frequency PWMs are converted', {
    pwm <- toPWM(pfm, type = 'log2probratio', bg = uniform)
    p <- 2^as.matrix(pwm) * 0.25
    pwm_ln <- pwm
    pwm_ln@profileMatrix <- log(p / 0.25)
    expect_equal(unname(Pando:::convert_pwm_fixed(pwm_ln, bg)), unname(log(p / bg)))

    pwm_prob <- toPWM(pfm, type = 'prob')
    expect_equal(unname(Pando:::convert_pwm_fixed(pwm_prob, bg)),
                 unname(log2(as.matrix(pwm_prob) / bg)))
})

test_that('uniform background leaves the PWM unchanged, as in motifmatchr', {
    pwm <- toPWM(pfm, type = 'log2probratio', bg = uniform)
    expect_equal(Pando:::convert_pwm_fixed(pwm, uniform), as.matrix(pwm))
    expect_equal(Pando:::convert_pwm_fixed(pwm, uniform),
                 motifmatchr:::convert_pwm(pwm, uniform))
})

test_that('background frequencies are validated and ordered by name', {
    expect_equal(Pando:::check_bg_fixed(c(T = 0.2, G = 0.3, C = 0.3, A = 0.2)), bg)
    expect_error(Pando:::check_bg_fixed(c(0.5, 0.5, 0.5, 0.5)))
    expect_error(Pando:::check_bg_fixed(c(X = 0.25, C = 0.25, G = 0.25, T = 0.25)))
})
