# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

# check whether bootstrap samples have all observations in the bag
whichAllInBag <- function(n, samples) {
    seqN <- seq_len(n)
    m <- apply(samples, 2, function(i) length(seqN[-i]))  # number of out-of-bag observations
    which(m == 0)
}
