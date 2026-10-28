## gbd_step3b_tpaf.R —— treatment-adjusted PAF V2（RR折减框架, 替代V1线性版）
## 框架: Santomauro 2023 (Lancet Psychiatry, 焦虑averted burden) 适配抑郁→痴呆
##   RR_eff(c) = 1 + (RR_untr-1)×(1-η·c);  PAF_tr = Pe(RR_eff-1)/(1+Pe(RR_eff-1))
## 基线注意: 本脚本基线=PAF_tr(c_now)（现行覆盖下的归因）, 与V1的33.3万(未处理暴露)并存报告
suppressPackageStartupMessages(library(dplyr))
wk <- "/mnt/f/HRS/GBD/work"
paf <- readRDS(file.path(wk, "gbd_paf.rds"))
mat <- read.csv(file.path(wk, "mdd_mat_204_santomauro2024.csv"), stringsAsFactors = FALSE)

RR_untr <- 1.31   ## CHARLS 轨迹锁版 HR（保守; Wu 2026 meta=1.90 作敏感性）
YR <- 2023
c_now <- mat %>% filter(sex=="both", year==2021) %>% select(location=location_gbd, c=mat) %>% distinct()

d <- paf %>% filter(year==YR) %>% left_join(c_now, by="location") %>%
  filter(!is.na(c)) %>% mutate(c = c/100,
    ## 从锁版 PAF 反推 Pe（保持与主分析逐国一致）
    Pe = paf / ((RR_untr - 1) * (1 - paf)),
    Pe = pmax(Pe, 0))

t_paf <- function(Pe, c, eta) {
  rreff <- 1 + (RR_untr - 1) * (1 - eta * c)
  Pe * (rreff - 1) / (Pe * (rreff - 1) + 1)
}
etav <- c(0.2, 0.3, 0.5)
scen <- list(S1_统一到50 = function(c) pmin(pmax(c, 0.50), 1),   # ← 新增（取max: 已达50%的不倒退）
             WHO2030_x1.5 = function(c) pmin(c*1.5, 1),
             plus25pp     = function(c) pmin(c+0.25, 1),
             universal100 = function(c) rep(1, length(c)))

cat("国家数:", nrow(d), "\n\n== 现行覆盖(c_now)下的归因负担 (万 DALY) ==\n")
for (e in etav) cat(sprintf("η=%.1f: %.1f 万 (对比未处理暴露基线 %.1f 万)\n",
    e, sum(d$dalys*t_paf(d$Pe, d$c, e))/1e4, sum(d$attrib)/1e4))

res <- list()
for (sn in names(scen)) for (e in etav) {
  z <- d %>% mutate(c_new = scen[[sn]](c),
                    avoid = dalys * (t_paf(Pe, c, e) - t_paf(Pe, c_new, e)))
  res[[paste(sn, e)]] <- z %>% summarise(avoid_g = sum(avoid)) %>% mutate(scenario=sn, eta=e)
}
R <- bind_rows(res)
cat("\n== 相对现行覆盖 可规避归因DALY (万) ==\n")
print(R %>% mutate(avoid_wan = round(avoid_g/1e4, 1)) %>% select(scenario, eta, avoid_wan))

top <- d %>% mutate(c_new = scen[["WHO2030_x1.5"]](c),
                    avoid = dalys * (t_paf(Pe, c, 0.3) - t_paf(Pe, c_new, 0.3))) %>%
  select(location, avoid) %>% arrange(desc(avoid)) %>% head(12) %>%
  mutate(avoid = round(avoid))
cat("\n== Top12 (WHO2030×1.5, η=0.3) ==\n"); print(top, row.names = FALSE)

saveRDS(list(detail = d, RR = R, top = top, RR_untr = RR_untr), file.path(wk, "gbd_step3b_tpaf.rds"))
cat("\n已保存 work/gbd_step3b_tpaf.rds\n")
