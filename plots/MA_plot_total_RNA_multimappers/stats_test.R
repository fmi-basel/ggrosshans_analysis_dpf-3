# Dpf-3 null vs mut-2

## downregulated
q <- 835
m <- 1855
n <- 16041
k <- 4572
down <- phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE) 
print(down)


## upregulated
q <- 647
m <- 1264
n <- 16632
k <- 3454
up <- phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE) 
print(up)

# Dpf-3 null vs mut-7

## downregulated
q <- 752
m <- 1851
n <- 15708
k <- 4365
down <- phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE)
print(down)

## upregulated
q <- 644
m <- 1268
n <- 16291
k <- 3436
up <- phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE) 
print(up)

# mut-2 vs mut-7

## downregulated
q <- 2475
m <- 4550
n <- 13780
k <- 4365
down <- phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE)
print(down)

## upregulated
q <- 2760
m <- 4120
n <- 14210
k <- 4515
up <- phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE) 
print(up)

