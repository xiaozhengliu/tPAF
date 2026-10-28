## gbd_step4_bapc.R —— 2024-2050 痴呆 DALY 投影 + 覆盖率干预累计可规避量
## 逻辑：分国 log(DALY)~year 线性投影（数据不足15年的国用全球率收缩）
##       归因DALY = 投影DALY × PAF；可规避 = 归因 × (paf(c_now) - paf(c_target))
library(data.table)

wk <- "/mnt/f/HRS/GBD/work"; setwd(wk)

## ---- 1. 读主时间序列（All ages, Both, DALYs, number）----
files <- list.files("/mnt/f/HRS/GBD", pattern = "^ad.*\\.csv$", full.names = TRUE)
ing <- function(f) {
  d <- fread(f)
  setnames(d, tolower(names(d)))
  loccol <- grep("^location_name", names(d), value = TRUE)[1]
  d <- d[sex_name == "Both" & age_name == "All ages" &
         tolower(measure_name) %like% "daly" & tolower(metric_name) == "number"]
  d[, .(location = get(loccol), year, daly = val)]
}
ad <- rbindlist(lapply(files, ing))
ad <- ad[, .(daly = sum(daly)), by = .(location, year)]
yrs <- sort(unique(ad$year))
cat(">> 年份覆盖:", min(yrs), "-", max(yrs), " 国家数:", uniqueN(ad$location), "\n")

## ---- 2. 分国投影 log(DALY) = a + b*year ----
proj_years <- 2024:2050
fit_country <- function(d) {
  d <- d[order(year)]
  if (nrow(d) >= 15) {
    f <- lm(log(daly) ~ year, data = d)
    b <- coef(f)[2]
  } else b <- NA_real_
  list(b = b, last = d$year[nrow(d)], lastval = d$daly[nrow(d)])
}
ests <- ad[, fit_country(.SD), by = location]
b_glob <- ests[!is.na(b), sum(b * lastval) / sum(lastval)]
ests[is.na(b), b := b_glob]
cat(">> 全球平均年增长率:", round((exp(b_glob) - 1) * 100, 2), "%\n")

## 投影：从各国最后一年外推
proj <- ests[, {
  y0 <- proj_years[proj_years > last]
  .(year = y0, daly_proj = lastval * exp(b * (y0 - last)))
}, by = location]

## ---- 3. PAF 与覆盖率情景 ----
paf <- as.data.table(readRDS("gbd_paf.rds"))
pafcol <- grep("paf", names(paf), value = TRUE)[1]
p21 <- paf[year == max(year),
           .(location, paf = get(pafcol))]          ## 用最新年 PAF 代表趋势外推
proj[, location := as.character(location)]
p21[,  location := as.character(location)]
        
proj <- merge(proj, p21, by = "location")
cat(">> nrow proj after p21 merge:", nrow(proj), " paf range:", range(proj$paf, na.rm = TRUE), "\n")
tpaf <- function(pe, c, eta) pe * (RR - 1) * (1 - eta * c) /
        (1 + pe * (RR - 1) * (1 - eta * c))
RR <- 1.31; ETA <- 0.3

## 当前覆盖率（step3c 输出，若无则给默认 0.05）
matf <- file.path(wk, "gbd_step3c_doseresp.csv")
mm <- fread("/mnt/f/HRS/GBD/work/mdd_mat_204_santomauro2024.csv")
cnow <- mm[tolower(sex) == "both" & year == max(year),
           .(location = location_gbd, c_now = mat / 100)]
cnow[, location := as.character(location)]
cnow[, location := as.character(location)]
proj <- merge(proj, cnow, by = "location")
cat(">> nrow proj after cnow merge:", nrow(proj), " c_now range:", range(proj$c_now, na.rm = TRUE), "\n")
## 情景：A 现状不变；B 2035 年线性提至 50%（S1）；C 2035 至 2050 提至全民覆盖
proj[, `:=`(
  paf_now  = paf,
  attrib   = daly_proj * paf
)]
proj[, c_B := pmin(c_now + (0.50 - c_now) * pmin(pmax((year - 2023) / 12, 0), 1), 1)]
proj[, c_C := pmin(c_now + (1.00 - c_now) * pmin(pmax((year - 2023) / 12, 0), 1), 1)]

## 上面 avoid 的 tpaf 用法太绕，改为直接数值算：
pe_from_paf <- function(pafv) pafv / ((RR - 1) * (1 - pafv))
proj[, pe := pe_from_paf(paf_now)]
proj[, avoidB := attrib * (paf_now - tpaf(pe, c_B, ETA)) / paf_now]
proj[, avoidC := attrib * (paf_now - tpaf(pe, c_C, ETA)) / paf_now]
cat(">> attrib sum:", sum(proj$attrib, na.rm = TRUE), " daly_proj sum:", sum(proj$daly_proj, na.rm = TRUE), "\n")
ans <- proj[, .(
  归因DALY累计_万   = round(sum(attrib) / 1e4),
  情景B可规避_万    = round(sum(avoidB) / 1e4),
  情景C可规避_万    = round(sum(avoidC) / 1e4)
)]
cat("\n== 2024-2050 累计（全球，万 DALY）==\n"); print(ans)

## 分年轨迹（画图用）
traj <- proj[, .(attrib = sum(attrib), avoidB = sum(avoidB), avoidC = sum(avoidC)),
             by = year]
fwrite(traj, "gbd_step4_bapc_traj.csv")
saveRDS(list(proj = proj, traj = traj), "gbd_step4_bapc.rds")
cat("\n已保存 gbd_step4_bapc.rds / gbd_step4_bapc_traj.csv\n")
