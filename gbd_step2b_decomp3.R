## gbd_step2b_decomp3.R —— Das Gupta 对称分解: 3因子(人口/年龄结构/率) + 4因子(+PAF)
## 输入: gbd_age.rds (step0b 输出) + gbd_paf.rds (含 1990/2023 PAF)
## 逻辑: Number = P × Σ_a (S_a × R_a) [× PAF]；对每个国家做因子置换平均, 再全球加总
library(data.table)

wk <- "/mnt/f/HRS/GBD/work"; setwd(wk)
age <- as.data.table(readRDS("gbd_age.rds"))
paf <- as.data.table(readRDS("gbd_paf.rds"))

bands <- c("15-49 years", "50-69 years", "70+ years")

## ---- 组装国家-年份因子表 ----
pop_t <- age[grepl("^dalys", measure) & age %in% bands & !is.na(pop),
             .(P = sum(pop)), by = .(location, year)]
sh <- age[grepl("^dalys", measure) & age %in% bands & !is.na(pop),
          .(location, year, age, S = pop)]
rt <- age[grepl("^dalys", measure) & age %in% bands & !is.na(rate),
          .(location, year, age, R = rate)]
sh[, S := S / sum(S), by = .(location, year)]

sw <- dcast(sh, location + year ~ age, value.var = "S")
setnames(sw, bands, paste0("S", 1:3))
rw <- dcast(rt, location + year ~ age, value.var = "R")
setnames(rw, bands, paste0("R", 1:3))

pafcol <- grep("paf", names(paf), value = TRUE)[1]
pw <- paf[year %in% c(1990, 2023), .(location, year, paf = get(pafcol))]

dt <- Reduce(function(x, y) merge(x, y, by = c("location", "year")),
             list(pop_t, sw, rw, pw))
dt <- dt[complete.cases(dt) & P > 0]
stopifnot("国家数≠204" = uniqueN(dt$location) == 204)

## ---- 内部一致性校验: f3 应复现 All Ages DALYs Number ----
chk <- merge(dt[year == 2023],
             age[grepl("^dalys", measure) & age == "All ages", .(location, year, daly_n = number)],
             by = c("location", "year"))
chk[, recon := P * (S1*R1 + S2*R2 + S3*R3) / 1e5]
cat("== 重构校验 (2023, 万DALY) ==\n")
cat(sprintf("官方 Number 合计: %.0f  重构: %.0f  相对差: %.2f%%\n",
            sum(chk$daly_n), sum(chk$recon),
            100 * (sum(chk$recon) - sum(chk$daly_n)) / sum(chk$daly_n)))

## ---- Das Gupta 对称分解（自写置换, 不依赖 combinat） ----
permn <- function(v) {
  if (length(v) == 1) return(list(v))
  out <- list()
  for (i in seq_along(v)) {
    sub <- permn(v[-i])
    for (s in sub) out[[length(out) + 1]] <- c(v[i], s)
  }
  out
}

dg <- function(fac_cols, dat) {
  d0 <- as.matrix(dat[year == 1990, ..fac_cols])
  d1 <- as.matrix(dat[year == 2023, ..fac_cols])
  k <- length(fac_cols)
  pm <- permn(1:k)
  eff <- matrix(0, nrow = length(pm), ncol = k)
  f <- function(m) m[, "P"] * (m[, "S1"]*m[, "R1"] + m[, "S2"]*m[, "R2"] + m[, "S3"]*m[, "R3"]) / 1e5 *
         if ("paf" %in% colnames(m)) m[, "paf"] else 1
  for (p in seq_along(pm)) {
    cur0 <- d0; cur1 <- d0
    for (i in pm[[p]]) {
      cur1[, i] <- d1[, i]
      eff[p, i] <- sum(f(cur1)) - sum(f(cur0))
      cur0 <- cur1
    }
  }
  colMeans(eff)
}

res3 <- dg(c("P", "S1", "S2", "S3", "R1", "R2", "R3"), dt)
## 注意: 3因子=人口P/结构S(3年龄合并)/率R(3年龄合并)
agg3 <- c(Population = sum(res3[1]), AgeStructure = sum(res3[2:4]), Rate = sum(res3[5:7]))
tot3 <- sum(agg3)
obs3 <- sum(dt[year == 2023, P * (S1*R1 + S2*R2 + S3*R3) / 1e5]) -
        sum(dt[year == 1990, P * (S1*R1 + S2*R2 + S3*R3) / 1e5])
cat("\n== 3因子分解: 痴呆总 DALY 变化 (1990→2023) ==\n")
print(round(agg3))
cat(sprintf("分解合计 %.0f vs 实际差 %.0f (残差 %.4f%%)\n",
            tot3, obs3, 100 * (tot3 - obs3) / obs3))
cat(sprintf("占比: 人口 %.1f%% / 结构 %.1f%% / 率 %.1f%%\n",
            100 * agg3[1] / tot3, 100 * agg3[2] / tot3, 100 * agg3[3] / tot3))

res4 <- dg(c("P", "S1", "S2", "S3", "R1", "R2", "R3", "paf"), dt)
agg4 <- c(Population = res4[1], AgeStructure = sum(res4[2:4]),
          Rate = sum(res4[5:7]), PAF = res4[8])
tot4 <- sum(agg4)
obs4 <- sum(dt[year == 2023, P * (S1*R1 + S2*R2 + S3*R3) * paf / 1e5]) -
        sum(dt[year == 1990, P * (S1*R1 + S2*R2 + S3*R3) * paf / 1e5])
cat("\n== 4因子分解: 抑郁归因痴呆 DALY 变化 (1990→2023) ==\n")
print(round(agg4))
cat(sprintf("分解合计 %.0f vs 实际差 %.0f (残差 %.4f%%)\n",
            tot4, obs4, 100 * (tot4 - obs4) / obs4))
cat(sprintf("占比: 人口 %.1f%% / 结构 %.1f%% / 率 %.1f%% / PAF %.1f%%\n",
            100 * agg4[1] / tot4, 100 * agg4[2] / tot4,
            100 * agg4[3] / tot4, 100 * agg4[4] / tot4))

saveRDS(list(dt = dt, agg3 = agg3, agg4 = agg4),
        "gbd_step2b_decomp3.rds")
cat("\n已保存 gbd_step2b_decomp3.rds\n")
