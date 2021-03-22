q <- 71
m <- 120
n <- 2617
k <- 328

down <- phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE) 
print(down)


q <- 99
m <- 147
n <- 2590
k <- 399
down_small_up_total <- phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE) 
print(down_small_up_total)
