context("ffmanova")

data(dressing)

# add character coded variables (categorical)

dressing$pressCat <- as.character(dressing$press)
dressing$stabCat <- as.character(dressing$stab)
dressing$emulCat <- as.character(dressing$emul)
dressing$to <- 2
dressing$dayNum <- as.numeric(dressing$day)


test_that("Version_1_0_0", {
  mod1visc <- ffmanova(visc ~ (pressCat + stabCat + emulCat)^2 + day, data = dressing, returnModel = FALSE)
  mod2visc <- ffmanova(visc ~ (press + stab + emul)^2 + I(press^2) + I(stab^2) + I(emul^2) + day, data = dressing, returnModel = FALSE)
  mod1rheo <- ffmanova(rheo ~ (pressCat + stabCat + emulCat)^2 + day, data = dressing, stand = FALSE, returnModel = FALSE)
  mod2rheo <- ffmanova(rheo ~ (press + stab + emul)^2 + I(press^2) + I(stab^2) + I(emul^2) + day, data = dressing, returnModel = FALSE)
  
  mod1visc1 <- ffmanova:::ffmanova_Version_1_0_0(visc ~ (pressCat + stabCat + emulCat)^2 + day, data = dressing)
  mod2visc1 <- ffmanova:::ffmanova_Version_1_0_0(visc ~ (press + stab + emul)^2 + I(press^2) + I(stab^2) + I(emul^2) + day, data = dressing)
  mod1rheo1 <- ffmanova:::ffmanova_Version_1_0_0(rheo ~ (pressCat + stabCat + emulCat)^2 + day, data = dressing, stand = FALSE)
  mod2rheo1 <- ffmanova:::ffmanova_Version_1_0_0(rheo ~ (press + stab + emul)^2 + I(press^2) + I(stab^2) + I(emul^2) + day, data = dressing)
  
  expect_equivalent(mod1visc, mod1visc1)  # equivalent since missing y name fixed in new version
  expect_equivalent(mod1visc, mod1visc1)  # 
  expect_equal(mod1rheo, mod1rheo1)
  expect_equal(mod2rheo, mod2rheo1)
})

if (require(car)) test_that("car::Anova", {
  ff1 <- ffAnova(visc ~ (pressCat + stabCat + emulCat)^2 + day, data = dressing)
  ff2 <- ffAnova(visc ~ (press + stab + emul)^2 + day, data = dressing)
  ff0 <- ffAnova(visc ~ (press + stab + emul)^2 + dayNum - 1, data = dressing)
  lm1 <- lm(visc ~ (pressCat + stabCat + emulCat)^2 + day, data = dressing)
  lm2 <- lm(visc ~ (press + stab + emul)^2 + day, data = dressing)
  lm0 <- lm(visc ~ (press + stab + emul)^2 + dayNum - 1, data = dressing)
  ff1b <- ffAnova(lm1)
  ff2b <- ffAnova(lm2)
  car0 <- Anova(lm0)
  car1 <- Anova(lm1)
  car2 <- Anova(lm2)
  expect_equal(ff1, ff1b)
  expect_equal(ff2, ff2b)
  expect_equivalent(ff0[, names(car0)], car0)
  expect_equivalent(ff1[, names(car1)], car1)
  expect_equivalent(ff2[, names(car2)], car2)
})


test_that("lm predictions", {
  lm1 <- lm(visc ~ pressCat + (stabCat + emul + I(emul^2))^2 + dayNum, data = dressing)
  lm0 <- lm(visc ~ press + (stab + I(emul^2) + emul + I(emul^2))^2 + dayNum - 1, data = dressing)
  
  f1A <- ffmanova(cbind(visc, rheo, visc) ~ pressCat + (stabCat + emul + I(emul^2))^2 + dayNum, data = dressing)
  f1a <- ffmanova(visc ~ pressCat + (stabCat + emul + I(emul^2))^2 + dayNum, data = dressing)
  f1b <- ffmanova(lm1, stand = TRUE, returnYhat = TRUE, returnYhatStd = TRUE)
  f1c <- ffmanova(lm1, stand = FALSE)
  f0a <- ffmanova(visc ~ press + (stab + I(emul^2) + emul + I(emul^2))^2 + dayNum - 1, data = dressing)
  f0b <- ffmanova(lm0, stand = TRUE, returnYhat = TRUE, returnYhatStd = TRUE)
  f0c <- ffmanova(lm0, stand = FALSE)
  
  p1 <- predict(lm1, se.fit = TRUE)
  p0 <- predict(lm0, se.fit = TRUE)
  
  
  p1A <- predict(f1A)
  p1a <- predict(f1a)
  p1b <- predict(f1b)
  p1c <- predict(f1c)
  p0a <- predict(f0a)
  p0b <- predict(f0b)
  p0c <- predict(f0c)
  
  expect_equivalent(p1$fit, p1A$YnewPred[, 1])
  expect_equivalent(p1$fit, p1A$YnewPred[, 11])
  
  expect_equivalent(p1$fit, f1b$Yhat)
  expect_equivalent(p1$fit, p1a$YnewPred)
  expect_equivalent(p1$fit, p1b$YnewPred)
  expect_equivalent(p1$fit, p1c$YnewPred)
  
  expect_equivalent(p0$fit, f0b$Yhat)
  expect_equivalent(p0$fit, p0a$YnewPred)
  expect_equivalent(p0$fit, p0b$YnewPred)
  expect_equivalent(p0$fit, p0c$YnewPred)
  
  expect_equivalent(p1$se.fit, p1A$YnewStd[, 1])
  expect_equivalent(p1$se.fit, p1A$YnewStd[, 11])
  
  
  expect_equivalent(p1$se.fit, f1b$YhatStd)
  expect_equivalent(p1$se.fit, p1a$YnewStd)
  expect_equivalent(p1$se.fit, p1b$YnewStd)
  expect_equivalent(p1$se.fit, p1c$YnewStd)
  
  expect_equivalent(p0$se.fit, f0b$YhatStd)
  expect_equivalent(p0$se.fit, p0a$YnewStd)
  expect_equivalent(p0$se.fit, p0b$YnewStd)
  expect_equivalent(p0$se.fit, p0c$YnewStd)
})


test_that("simple model predictions with linComb", {
  day <- dressing$day
  dayC <- as.character(dressing$day)
  visc <- dressing$visc
  emul <- dressing$emul
  stab <- dressing$stab
  
  dayF1 <- data.frame(day = unique(day), dayC = unique(dayC))
  dayF2 <- dayF1
  dayF2$emul <- mean(emul)
  dayF2$stab <- mean(stab)
  dayF3 <- dayF2
  dayF3$stab <- 0.3
  dayF3$emul <- NA
  
  ffa <- ffmanova(visc ~ emul + stab + day, stand = FALSE, data = dressing)
  ffb <- ffmanova(lm(cbind(visc, visc) ~ emul + stab + day))
  ffc <- ffmanova(visc ~ emul + I(stab) + dayC)
  
  p1a <- predict(ffa, newdata = dayF1)
  p1b <- predict(ffb, newdata = dayF1)
  p1c <- predict(ffc, newdata = dayF1)
  p2a <- predict(ffa, newdata = dayF2)
  p2b <- predict(ffb, newdata = dayF2)
  p2c <- predict(ffc, newdata = dayF2)
  
  
  p3a <- predict(ffa, newdata = dayF3)
  p3c <- predict(ffc, newdata = dayF3)
  
  expect_equal(p1a, p2a)
  expect_equal(p1a, p1c)
  expect_equal(p1a, p2c)
  expect_equal(p1b, p2b)
  
  expect_equal(p3a, p3c)
  
  p1b23 <- predict(ffb, newdata = dayF1[2:3, ])
  expect_equal(p1b23[[1]], p1b[[1]][2:3, , drop = FALSE])
  expect_equal(p1b23[[2]], p1b[[2]][2:3, , drop = FALSE])
  
  p1b[[1]] <- p1b[[1]][, 2, drop = FALSE]
  p1b[[2]] <- p1b[[2]][, 1, drop = FALSE]
  expect_equal(p1b, p1a)
  
  p1c23 <- predict(ffc, newdata = dayF1[2:3, ])
  expect_equal(p1c23[[1]], p1c[[1]][2:3, , drop = FALSE])
  expect_equal(p1c23[[2]], p1c[[2]][2:3, , drop = FALSE])
  
  emul_ <- dressing$emul
  stab_ <- dressing$stab
  lma_ <- lm(visc ~ emul_ + stab_ + day, data = dressing)
  ffa_ <- ffmanova(lma_)
  emulF <- data.frame(emul = c(0.2, 0.1), stab = NA, emul_ = c(0.2, 0.1), stab_ = 0.3)
  stabF <- data.frame(emul = NA, stab = unique(dressing$stab), emul_ = 0.1, stab_ = unique(dressing$stab))
  esF <- rbind(emulF, stabF)
  
  cC <- t(matrix(c(10, -10, 0, 0, 0, 0, 0, -5, 0, 5), 5, 2))
  rownames(cC) <- c("emul", "stab")
  
  expect_error(predict(ffa, newdata = esF, cC))
  
  pcC <- predict(ffa_, newdata = esF, cC)
  expect_equal(cC %*% predict(ffa, newdata = esF)$YnewPred, pcC$YnewPred)
  
  expect_equivalent(as.data.frame(summary(lma_)["coefficients"][[1]][2:3, 1:2]), as.data.frame(pcC))
  
  
})



test_that("Advanced model predictions", {
  dayC <- as.character(dressing$day)
  
  l1 <- lm(visc ~ (pressCat + stabCat + emulCat)^2 + day, data = dressing)
  m1 <- ffmanova(visc ~ (I(pressCat) + as.character(stabCat) + emulCat)^2 + dayC, data = dressing)
  
  
  pl1 <- predict(l1, se.fit = TRUE)
  pm1 <- predict(m1)
  
  expect_equivalent(pl1$fit, pm1$YnewPred)
  expect_equivalent(pl1$se.fit, pm1$YnewStd)
  
  
  m2 <- ffmanova(visc ~ day + (pressCat + stabCat + emulCat)^2 - 1, data = dressing)
  pm2 <- predict(m2)
  expect_equal(pm1, pm2)
  
  
  m3 <- ffmanova(visc ~ to + to:((pressCat + stabCat + emulCat)^2 + day) - 1, data = dressing, stand = FALSE)
  pm3 <- predict(m3)
  expect_equal(pm1, pm3)
  
  m4 <- ffmanova(visc ~ (press + stab + emul)^2 + (I(press^2) + I(stab^2) + I(emul^2))^2 + (press + stab):I(emul^2) + (press + emul):I(stab^2) + (stab + emul):I(press^2) + day, data = dressing)
  pm4 <- predict(m4)
  expect_equal(pm1, pm4)
  
  
  m5 <- ffmanova(visc ~ day + (press + stab + emul)^2 + (I(press^2) + I(stab^2) + I(emul^2))^2 + (press + stab):I(emul^2) + (press + emul):I(stab^2) + (stab + emul):I(press^2) - 1, 
                 data = dressing)
  pm5 <- predict(m5)
  expect_equal(pm1, pm5)
  
  
  dayF <- data.frame(day = dressing$day[2:4], dayC = dayC[2:4], to = 2)
  
  
  p_m1_dayF <- predict(m1, dayF)
  expect_equal(predict(m2, dayF), p_m1_dayF)
  expect_equal(predict(m3, dayF), p_m1_dayF)
  expect_equal(predict(m4, dayF), p_m1_dayF)
  expect_equal(predict(m5, dayF), p_m1_dayF)
  
  
  emulF <- data.frame(emul = unique(dressing$emul), emulCat = unique(dressing$emulCat), to = 2)
  emulDayF <- emulF
  emulDayF$day <- dressing$day[2:4]
  emulDayF$dayC <- dayC[2:4]
  
  expect_equal(predict(m1, emulF), predict(m3, emulF))
  expect_equal(predict(m1, emulDayF), predict(m3, emulDayF))
  expect_false(anyNA(predict(m2, dayF)))
  expect_false(anyNA(predict(m2, emulDayF)))
  expect_false(anyNA(predict(m4, dayF)))
  expect_false(anyNA(predict(m4, emulDayF)))
  expect_false(anyNA(predict(m5, dayF)))
  expect_false(anyNA(predict(m5, emulDayF)))
  
  
  expect_equal(length(sort((unique(round(c(ffAnova(m1)[, 4], ffAnova(m2)[, 4], ffAnova(m3)[, 4]), 6))))), 9)
  expect_equal(length(sort((unique(round(c(ffAnova(m4)[, 4], ffAnova(m5)[, 4]), 6))))), 20)
  
  
})


test_that("Embedded matrix as regressor", {
  dayC <- as.character(dressing$day)
  rhe1 <- dressing$rhe[, 1]
  rhe2 <- dressing$rhe[, 2]
  rhe12 <- cbind(rhe1, rhe2)
  
  l1 <- lm(visc ~ day + rhe12 + I(rhe12^2), data = dressing)
  m1 <- expect_warning(ffmanova(visc ~ day + rhe12 + I(rhe12^2), data = dressing))
  m2 <- ffmanova(visc ~ day + rhe1 + rhe2 + I(rhe1^2) + I(rhe2^2), data = dressing)
  m1b <- expect_warning(ffmanova(l1))
  expect_identical(m1, m1b)
  
  expect_equivalent(anova(l1)[2:3, ], ffAnova(m1)[2:3, ])
  if (require(car)) {
    expect_equivalent(Anova(l1)[c(1, 3, 4), ], ffAnova(m1)[c(1, 3, 4), c(2, 1, 4, 5)])
  }
  
  r12 <- data.frame(rhe1 = m2$ffModel$xlev$rhe1, rhe2 = m2$ffModel$xlev$rhe2)
  
  expect_equal(predict(m1, data.frame(rhe12 = I(as.matrix(r12)))), predict(m2, r12))
  r12$rhe2 = rev(r12$rhe2)
  expect_equal(predict(m1, data.frame(rhe12 = I(m1$ffModel$xlev$rhe12))), predict(m2, r12))
  
  
  
})