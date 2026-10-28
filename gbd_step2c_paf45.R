## gbd_step2c_paf45.R —— 45+ 年龄别 Levin PAF 口径修正（对齐 Wu 2026）
## 输入：mdd_age_1.csv（MDD prevalence, 分年龄带）+ ad_age_1.csv（痴呆DALY, 分年龄带）+ gbd_paf.rds
## 逻辑：PAF_a = Pe_a(RR-1)/(1+Pe_a(RR-1))；45+ PAF = Σ(DALY_a × PAF_a)/Σ DALY_a
library(data.table)

wk <- "/mnt/f/HRS/GBD/work"; setwd(wk)
RR <- 1.31
bands45 <- c("50-69 years", "70+ years")   ## 45+ 的可用带（GBD带边界，50-69近似45-69）

ing <- function(f, cause) {
  d <- fread(f)
  setnames(d, tolower(names(d)))
  d <- d[sex_name == "Both" & age_name != "All ages"]
  d[, .(location = location_name, measure = tolower(measure_name),
        metric = tolower(metric_name), age = age_name, year, val)]
}

mdd <- ing("/mnt/f/HRS/GBD/mdd_age_1.csv")
ad  <- ing("/mnt/f/HRS/GBD/ad_age_1.csv")

## MDD 患病率 Pe_a（每 1 万人 → 比例）
pe <- mdd[measure == "prevalence" & metric == "rate" & age %in% bands45,
          .(location, year, age, pe = val / 1e5)]
stopifnot("缺 prevalence rate" = nrow(pe) > 0)
cat(">> pe rows:", nrow(pe), " 范围:", round(range(pe$pe, na.rm = TRUE), 4), "\n")

## 各带痴呆 DALY 数（权重）
daly <- ad[grepl("^dalys", measure) & metric == "number" & age %in% bands45,
           .(location, year, age, daly_n = val)]

## 分带 Levin PAF
pe[, paf_a := pe * (RR - 1) / (1 + pe * (RR - 1))]
dt <- merge(pe, daly, by = c("location", "year", "age"))

## 45+ 合成（按 DALY 加权）
p45 <- dt[, .(paf45_raw = sum(daly_n * paf_a) / sum(daly_n),
              daly45    = sum(daly_n)), by = .(location, year)]

## 治疗调整（与主口径对齐）：RR_eff = 1+(RR-1)(1-ηc)，η=0.3，c=各国当前覆盖率
ETA <- 0.3
mm <- fread("/mnt/f/HRS/GBD/work/mdd_mat_204_santomauro2024.csv")
cnow <- mm[tolower(sex) == "both" & year == max(year),
           .(location = as.character(location_gbd), c = mat / 100)]
p45 <- merge(p45, cnow, by = "location")

## 逐带治疗调整后再加权（正确做法）
dt2 <- merge(dt, cnow, by = "location")
dt2[, paf_a_adj := pe * (RR - 1) * (1 - ETA * c) /
      (1 + pe * (RR - 1) * (1 - ETA * c))]
p45b <- dt2[, .(paf45_adj = sum(daly_n * paf_a_adj) / sum(daly_n)),
            by = .(location, year)]
p45 <- merge(p45[, .(location, year, paf45_raw, daly45)],
             p45b, by = c("location", "year"))

## 现有全人群 t-PAF 对比
paf <- as.data.table(readRDS("gbd_paf.rds"))
cmp <- merge(p45, paf[, .(location, year, paf_all = paf)],
             by = c("location", "year"))

## 全球 DALY 加权汇总
g <- cmp[, .(paf45_raw = sum(daly45 * paf45_raw) / sum(daly45),
             paf45_adj = sum(daly45 * paf45_adj) / sum(daly45),
             paf_all   = sum(daly45 * paf_all)   / sum(daly45)),
         by = year]

cat("\n== 全球口径对比（45+ Levin vs 现全人群 t-PAF）==\n")
print(round(g * 100, 2))  ## 百分比
cat("\n>> 分年 range paf45_adj:", round(range(cmp$paf45_adj), 4), "\n")

fwrite(cmp, "gbd_step2c_paf45_cmp.csv")
saveRDS(cmp, "gbd_step2c_paf45.rds")
cat("\n已保存 gbd_step2c_paf45.rds / _cmp.csv\n")
