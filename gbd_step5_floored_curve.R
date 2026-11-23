##cat > /mnt/f/HRS/GBD/work/gbd_step5_floored_curve.R << 'REOF'
## gbd_step5_floored_curve.R —— 情景地板化后的全球剂量反应曲线 v1.0
## 背景：统一目标情景（"所有国家达到 13.65%"）把高覆盖国压向低目标，
##       全球聚合在 13.65% 处出现负避免（-368/年）——凹函数 artifact。
## 修正：情景目标按各国当前覆盖取 max（已在目标以上的国家不动）。
## 影响：S1(50%)/universal 不变（最高覆盖 36.2% < 50%）；仅 13.65% 参照点
##       与 Panel A 负值区改变。负值区删除后图形更干净。
wk <- "/mnt/f/HRS/GBD/work"; setwd(wk)
library(data.table)
RR <- 1.31
tpaf <- function(pe, c, eta) { x <- pe * (RR - 1) * (1 - eta * c); x / (1 + x) }

paf <- as.data.table(readRDS("gbd_paf.rds"))[year == 2023]
mat <- fread("mdd_mat_204_santomauro2024.csv")
cnow <- mat[tolower(sex) == "both" & year == max(year),
            .(location = location_gbd, c_now = mat / 100)]
paf[, pe := paf / ((RR - 1) * (1 - paf))]
paf <- merge(paf, cnow, by = "location")
stopifnot(!any(is.na(paf$pe)), !any(is.na(paf$c_now)))

grid <- seq(0, 1, 0.02)
out <- rbindlist(lapply(c(0.2, 0.3, 0.5), function(eta_v) {
  data.table(eta = eta_v, panel = "Global", c_target = grid,
             avoid_daly = sapply(grid, function(ct)
               sum(paf$dalys * (tpaf(paf$pe, paf$c_now, eta_v) -
                                tpaf(paf$pe, pmax(paf$c_now, ct), eta_v)))))
}))
fwrite(out, "gbd_step3c_doseresp_floored.csv")

v13 <- sapply(c(0.2, 0.3, 0.5), function(eta_v)
  sum(paf$dalys * (tpaf(paf$pe, paf$c_now, eta_v) -
                   tpaf(paf$pe, pmax(paf$c_now, 0.1365), eta_v))))
cat("地板化后 1.5xWHO (13.65%) 年避免 DALYs (eta 0.2/0.3/0.5):", round(v13), "\n")

chk50 <- sum(paf$dalys * (tpaf(paf$pe, paf$c_now, 0.3) - tpaf(paf$pe, pmax(paf$c_now, 0.50), 0.3)))
chk100 <- sum(paf$dalys * (tpaf(paf$pe, paf$c_now, 0.3) - tpaf(paf$pe, 1, 0.3)))
stopifnot(abs(chk50 - 35690) < 500, abs(chk100 - 85406) < 500)
cat("校验通过：地板化后 S1 =", round(chk50), "（锁版 35,690）；universal =", round(chk100), "（锁版 85,406）\n")
##REOF
