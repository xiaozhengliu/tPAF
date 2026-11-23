## gbd_step5_rr_sensitivity.R —— RR 敏感性（审稿必问项）v1.1
## 逻辑：pe（GBD MDD 患病率）由 step3b 用 RR=1.31 从锁版 untreated PAF 反推，
##       该值不随 RR 变化；RR 只进入下游 paf/tpaf 公式。故对每个 RR' 精确重算，
##       非近似。RR'=1.31 行必须复现锁版数字（脚本内置校验）。
wk <- "/mnt/f/HRS/GBD/work"; setwd(wk)
library(data.table)

RR0  <- 1.31                       # 与 step3b 一致的 pe 反推用 RR
ETA  <- 0.3
rr_grid <- c(1.31, 1.5, 2.0, 2.2)  # 1.31 主分析 + 敏感性（覆盖 Commission 的 2.2）

## ---- 输入 ----
paf <- as.data.table(readRDS("gbd_paf.rds"))[year == 2023]
stopifnot(all(c("location","paf","dalys") %in% names(paf)))
mat <- fread("mdd_mat_204_santomauro2024.csv")
cnow <- mat[tolower(sex) == "both" & year == max(year),
            .(location = location_gbd, c_now = mat / 100)]
proj <- readRDS("gbd_step4_bapc.rds")$proj   # location/year/daly_proj/c_now...

## pe 反推（一次性；与 step3b 同式）
paf[, pe := paf / ((RR0 - 1) * (1 - paf))]
paf <- merge(paf, cnow, by = "location")
stopifnot(!any(is.na(paf$pe)), !any(is.na(paf$c_now)))

tpaf <- function(pe, RR, c) {
  x <- pe * (RR - 1) * (1 - ETA * c)
  x / (1 + x)
}

## 投影情景 ramp（与 step4 相同：2023→2035 线性）
proj[, c_B := pmin(c_now + (0.50 - c_now) * pmin(pmax((year - 2023) / 12, 0), 1), 1)]
proj[, c_C := pmin(c_now + (1.00 - c_now) * pmin(pmax((year - 2023) / 12, 0), 1), 1)]

eps <- 0.001
out <- rbindlist(lapply(rr_grid, function(RR) {
  ## 横断面（2023）：全部在 paf 的 [i, j] 作用域内计算
  paf[, p_untr := pe * (RR - 1) / (1 + pe * (RR - 1))]
  res <- paf[, .(
    paf_global_pct     = round(100 * sum(dalys * p_untr) / sum(dalys), 2),
    avert50_annual_2023  = sum(dalys * (tpaf(pe, RR, c_now) - tpaf(pe, RR, 0.50))),
    avert100_annual_2023 = sum(dalys * (tpaf(pe, RR, c_now) - tpaf(pe, RR, 1))),
    marginal_perpp     = sum(dalys * (tpaf(pe, RR, pmin(c_now + eps, 1)) - tpaf(pe, RR, c_now))) / eps
  )]
  ## 投影（2024-2050 累计）：全部在 proj 的作用域内计算
  res2 <- proj[, .(
    cum_attrib_2024_50 = sum(daly_proj * pe * (RR - 1) / (1 + pe * (RR - 1))),
    cum_avertB_2024_50 = sum(daly_proj * (pe * (RR - 1) / (1 + pe * (RR - 1)) - tpaf(pe, RR, c_B))),
    cum_avertC_2024_50 = sum(daly_proj * (pe * (RR - 1) / (1 + pe * (RR - 1)) - tpaf(pe, RR, c_C)))
  )]
  cbind(data.table(RR = RR), res, res2)
}))

fwrite(out, "rr_sensitivity.csv")
print(out[, lapply(.SD, function(x) if (is.numeric(x)) round(x) else x)])

## ---- 校验：RR=1.31 行必须复现锁版数（50% 年避免 35.7 万；累计 B 4.0M / C 7.9M）----
chk <- out[abs(RR - 1.31) < 1e-9]
stopifnot(abs(chk$avert50_annual_2023 / 1e4 - 35.7) < 1,
          abs(chk$cum_avertB_2024_50 / 1e6 - 4.0) < 0.15)
cat(">> 校验通过：RR=1.31 行复现锁版头条数\n")
