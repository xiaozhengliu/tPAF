##cat > /mnt/f/HRS/GBD/work/gbd_step5_sdi.R << 'REOF'
## gbd_step5_sdi.R —— SDI 五分位分层 v1.4（真修复版）
## v1.4：值列改名 sdiv（列名 "sdi" 会在 sdi[, j] 内遮蔽数据表 sdi 本体——上轮报错根因）；
##        删除冗余 qsrc 行；其余同 v1.3（官方标签优先 + 切点近似补齐 + 全图抽查 CSV）。
wk <- "/mnt/f/HRS/GBD/work"; setwd(wk)
library(data.table)

SDI_FILE <- "IHME_GBD_2023_SDI_1950_2023_SDI_QUINTILES_Y2025M10D12.csv"
RR <- 1.31; ETA <- 0.3
tpaf <- function(pe, c) { x <- pe * (RR - 1) * (1 - ETA * c); x / (1 + x) }

det <- as.data.table(readRDS("gbd_step3b_tpaf.rds")$detail)

sdi <- tryCatch(fread(SDI_FILE, encoding = "UTF-8"),
                error = function(e) fread(SDI_FILE, encoding = "Latin-1"))
setnames(sdi, c("loc_id", "location", "sdiv", "quin"))
cat("总行数:", nrow(sdi), "；quintile 非空行:", sum(sdi$quin != ""), "\n")

sdi[location %like% "rkiye", location := "Türkiye"]
sdi[location %like% "Ivoire", location := "Côte d'Ivoire"]
det_locs <- unique(det$location)
miss <- setdiff(det_locs, sdi$location)
for (um in miss) {
  cand <- sdi$location[agrep(um, sdi$location, max.distance = 0.15, ignore.case = TRUE)]
  cand <- setdiff(cand, det_locs)
  if (length(cand) == 1) {
    sdi[location == cand, location := um]
    cat("模糊修复:", cand, "->", um, "\n")
  } else {
    cat("!! 未能匹配:", um, "候选:", paste(cand, collapse = ", "), "\n")
  }
}
stopifnot(all(det_locs %in% sdi$location))

sdi[quin == "", quin := NA_character_]
qord <- c("Low SDI", "Low-middle SDI", "Middle SDI", "High-middle SDI", "High SDI")
lab <- sdi[!is.na(quin)]
lab[, qf := factor(quin, levels = qord, ordered = TRUE)]
bounds <- lab[, .(lo = min(sdiv), hi = max(sdiv)), by = qf][order(qf)]
cuts <- (bounds$hi[1:4] + bounds$lo[2:5]) / 2
cat("由已标注行学到的近似切点:", paste(round(cuts, 4), collapse = " | "), "\n")
cat("（官方标签为参考切点制， Georgia/Kazakhstan 等即例外；以下补齐仅为近似）\n")
sdi[is.na(quin), quin := as.character(cut(sdiv, breaks = c(-Inf, cuts, Inf), labels = qord))]

m <- merge(det, sdi[, .(location, quin)], by = "location")
m[, quin := factor(quin, levels = qord, ordered = TRUE)]

res <- m[, .(
  countries         = .N,
  coverage_wmean    = round(100 * weighted.mean(c, dalys), 1),
  attrib_dalys_2023 = round(sum(attrib)),
  avoid50_annual    = round(sum(dalys * (tpaf(pe, c) - tpaf(pe, 0.50)))),
  avoid100_annual   = round(sum(dalys * (tpaf(pe, c) - tpaf(pe, 1))))
), by = quin][order(quin)]

age <- as.data.table(readRDS("gbd_age.rds"))
bands <- c("15-49 years", "50-69 years", "70+ years")
pop <- age[grepl("^dalys", measure) & age %in% bands & !is.na(pop),
           .(pop15 = sum(pop)), by = "location"]
m2 <- merge(m, pop, by = "location")
res <- merge(res,
  m2[, .(avoid50_per100k_15plus = round(1e5 * sum(dalys * (tpaf(pe, c) - tpaf(pe, 0.50))) /
                                        sum(pop15))), by = quin],
  by = "quin")[order(quin)]

fwrite(res, "Table_S6_SDI_quintiles.csv")
fwrite(sdi[location %in% det_locs][order(sdiv)], "SDI_mapping_204_countries.csv")
cat("\n==== Table S6 ====\n"); print(res)
cat("\n204 国映射已存 SDI_mapping_204_countries.csv —— 请重点核对：\n")
print(sdi[location %in% c("China","India","Brazil","Indonesia","Iran (Islamic Republic of)",
      "Egypt","South Africa","Mexico","Philippines","Kenya","Nigeria","Pakistan","Ethiopia","Sudan"),
      .(location, sdiv, quin)][order(sdiv)])
REOF
