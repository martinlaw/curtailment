library(librarian)
shelf(devtools, flextable)
load_all()

#d1 <- singlearmDesign(n.max=c(48, 50), p0=0.1, p1=0.242, alpha=0.05, power=0.8)
# More aggressive stopping:
# futility stopping: increase maxthetaF
# efficacy stopping: decrease minthetaE

#d2 <- singlearmDesign(n.max=c(48, 50), p0=0.1, p1=0.242, alpha=0.05, power=0.8, maxthetaF=1)


two.stage <- find2stageDesigns(n.max = 47:48, # No thetaF constraint
                               p0 = 0.1,
                               p1 = 0.242,
                               power = 0.8,
                               alpha = 0.05,
                               benefit=TRUE)
two.stage$all.des[6, ] # n1=20, r1=2
chosen.des <- two.stage$all.des[6, ]
n1 <- chosen.des$n1
n <- chosen.des$n

# two.stage.thetaF <- find2stageDesigns(n.max = 48:49,
#                                p0 = 0.1,
#                                p1 = 0.242,
#                                power = 0.8,
#                                alpha = 0.05,
#                                benefit=TRUE,
#                                maxthetaF=thetaF)
# two.stage.thetaF$all.des


conts.monitoring <- singlearmDesign(n.max = c(48,49),
                                    p0 = 0.1,
                                    p1 = 0.242,
                                    power = 0.8,
                                    alpha = 0.05)
conts.monitoring$all.des
thetaF <- conts.monitoring$all.des[, "thetaF"]

bounds <- drawDiagram(conts.monitoring)$bounds.mat
fail.bounds <- bounds$fail
# Replace fully sequential futility boundaries at n1 and n with two-stage boundaries:
fail.bounds[n1] <- chosen.des$r1
fail.bounds[n] <- chosen.des$r



# Test of "fully sequential in blocks" idea:
set.seed(21052026)
nsim <- 10000
p0 <- 0.1
p1 <- 0.242

h0 <- matrix(data=rbinom(n*nsim, 1, prob=p0), nrow=nsim)
h1 <- matrix(data=rbinom(n*nsim, 1, prob=p1), nrow=nsim)

p_reject <- function(y, n1, r1, e1, n, r, cts.fail.bounds){
  s <- t(apply(y, 1, cumsum)) # Cumulative successes

  # Replace fully sequential futility boundaries at n1 and n with two-stage boundaries:
  fail.bounds[n1] <- r1
  fail.bounds[n] <- r

  # Two stage:
  stage1.fail <- s[,n1]<=r1 # Fail at interim
  stage1.reject <- s[,n1]>e1 # Stop for efficacy at interim
  stage1.continue <- !stage1.fail & !stage1.reject # Do not stop at interim
  stage2.fail <- stage1.continue & s[,n]<=r # Fail at end (after not stopping at interim)
  stage2.reject <- stage1.continue & s[,n]>r # Reject H0 at end (after not stopping at interim)

  # Sequential in blocks:
  # What to do when futility boundary is crossed but we would reject H0 under two-stage boundaries? For now, reject:
  conts.stage1.reject <- stage1.reject
  conts.stage1.fail <- !stage1.reject & apply(s[, 1:n1], 1, function(x) any(x <= fail.bounds[1:n1])) # S1 failure under "conts monitoring"
  conts.stage1.continue <- !conts.stage1.fail & !conts.stage1.reject # Do not stop at interim (for failure or efficacy)
  conts.stage2.reject <- conts.stage1.continue & s[,n]>r # Reject H0 regardless of of any crossing of "conts" futility boundaries.
  conts.stage2.fail <- !conts.stage2.reject & conts.stage1.continue & apply(s[,(n1+1):n], 1, function(x) any(x <= fail.bounds[(n1+1):n])) # S2 failure under "conts monitoring".

  reject.2stage <- stage1.reject | stage2.reject
  reject.conts <- conts.stage1.reject | conts.stage2.reject

  ess.2stage <- (sum(n*stage1.continue) + sum(n1*!stage1.continue))/nsim
  ess.conts <- (sum(n*conts.stage1.continue) + sum(n1*!conts.stage1.continue))/nsim

  # Sanity checking:
  two.stage.results <- cbind(stage1.fail, stage1.reject, stage2.fail, stage2.reject)
  if(any(rowSums(two.stage.results)!=1)) stop("Ambiguous result in two-stage design")
  conts.results <- cbind(conts.stage1.reject, conts.stage1.fail, conts.stage2.fail, conts.stage2.reject)
  if(any(rowSums(conts.results)!=1)) stop("Ambiguous result in sequential in blocks design")

  return(data.frame(reject=c(mean(reject.2stage), mean(reject.conts)), ESS=c(ess.2stage, ess.conts)))
}
h0 <- p_reject(y=h0, n1=chosen.des$n1, r1=chosen.des$r1, e1=chosen.des$e1, n=chosen.des$n, r=chosen.des$r, cts.fail.bounds=fail.bounds)
h1 <- p_reject(y=h1, n1=chosen.des$n1, r1=chosen.des$r1, e1=chosen.des$e1, n=chosen.des$n, r=chosen.des$r, cts.fail.bounds=fail.bounds)

results <- cbind(h0, h1)
results <- round(results, 3)
rownames(results) <- c("Two stage", "Seq. in blocks")
results <- cbind(rownames(results), results)
names(results) <- c("Design", "T1ER", "ESSH0", "Power", "ESSH1")

qflextable(results)

# oc <- cbind(stage1.fail, stage1.reject, stage1.continue, stage2.fail, stage2.reject,
#             conts.stage1.fail, conts.stage1.reject, conts.stage1.continue, conts.stage2.fail, conts.stage2.reject,
#             t1, conts.t1)



# Next step:
#
# Quantify the proportion of failures (at stage 1 and at stage 2) that would occur under conts.
# monitoring but not under the two-stage design. Is there a benefit in terms of ESSH0?

# T1ER
t1.conts <- !conts.stage1.fail & !conts.stage2.fail
t1.two.stage <- !stage1.fail & !stage2.fail
