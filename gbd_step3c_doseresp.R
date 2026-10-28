## gbd_step3c_doseresp.R —— v3 English labels (2026-10-21, locked)
## Continuous treatment-coverage -> avoidable dementia DALY dose-response
## Presentation follows BMJ 2026;393 (bmj-2026-100561): continuous dose axis,
## policy target zone shading, per-country facet panels, threshold lines
## Inputs : gbd_paf.rds (year==2023, 204 countries), mdd_mat_204_santomauro2024.csv
## Outputs: gbd_step3c_doseresp.csv, gbd_step3c_doseresp.pdf
library(tidyverse)
library(gridExtra)

rr_untr <- 1.31
eta_grid <- c(0.2, 0.3, 0.5)

paf <- readRDS("gbd_paf.rds") %>% filter(year == 2023)
mat <- read.csv("mdd_mat_204_santomauro2024.csv", stringsAsFactors = FALSE)

c_now <- mat %>% filter(year == 2021, sex == "both") %>%
  transmute(location = location_gbd, c_now = mat / 100)

df <- paf %>% left_join(c_now, by = "location") %>% filter(!is.na(c_now))
stopifnot("must match 204/204" = nrow(df) == 204)

dalycol <- grep("daly|attrib", names(df), value = TRUE, ignore.case = TRUE)[1]
stopifnot("no attributable-DALY column found" = !is.na(dalycol))
df <- df %>% rename(attrib_daly = all_of(dalycol))

df <- df %>% mutate(Pe = paf / ((rr_untr - 1) * (1 - paf)))

tpaf <- function(Pe, c, eta) {
  rreff <- 1 + (rr_untr - 1) * (1 - eta * c)
  x <- Pe * (rreff - 1)
  x / (1 + x)
}

c_grid <- seq(0, 1, by = 0.02)

res <- expand_grid(i = seq_len(nrow(df)), eta = eta_grid, c_target = c_grid) %>%
  mutate(
    location = df$location[i],
    Pe = df$Pe[i], c_now = df$c_now[i], attrib_daly = df$attrib_daly[i],
    paf_t = tpaf(Pe, c_target, eta),
    paf_now = tpaf(Pe, c_now, eta),
    avoid_daly = attrib_daly * (paf_now - paf_t)   # positive = avoidable
  )

glob <- res %>% group_by(eta, c_target) %>%
  summarise(avoid_daly = sum(avoid_daly), .groups = "drop") %>% mutate(panel = "Global")

top4 <- res %>% filter(location %in% c("China", "United States of America", "India", "Japan")) %>%
  select(location, eta, c_target, avoid_daly) %>% mutate(panel = location)

write.csv(bind_rows(glob, top4), "gbd_step3c_doseresp.csv", row.names = FALSE)

c_now_global <- 9.1 / 100
who2030 <- c_now_global * 1.5
s1 <- 0.50

eps <- 0.001
m <- df %>% mutate(
  a1 = attrib_daly * (tpaf(Pe, c_now, 0.3) - tpaf(Pe, pmin(c_now + eps, 1), 0.3))
) %>% summarise(v = sum(a1) / eps) %>% pull(v)

lbl_eta <- c(`0.2` = "efficacy eta = 0.2", `0.3` = "eta = 0.3", `0.5` = "eta = 0.5")

## ---- Panel A: Global ----
pA <- ggplot(glob, aes(c_target * 100, avoid_daly / 1e4, color = factor(eta), linetype = factor(eta))) +
  annotate("rect", xmin = who2030 * 100, xmax = s1 * 100,
           ymin = -Inf, ymax = Inf, alpha = 0.12, fill = "grey50") +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = c_now_global * 100, linetype = 3, color = "grey40") +
  geom_vline(xintercept = who2030 * 100, linetype = 4, color = "grey40") +
  geom_vline(xintercept = s1 * 100, linetype = 2, color = "grey40") +
  annotate("text", x = who2030 * 100, y = Inf, label = "Policy target zone",
           vjust = 1.6, hjust = -0.05, size = 3, fontface = "italic", color = "grey40") +
  annotate("text", x = c_now_global * 100, y = Inf, label = "Current global coverage 9.1%",
           vjust = 1.5, hjust = -0.05, size = 3, color = "grey40") +
  annotate("text", x = who2030 * 100, y = Inf, label = "WHO 2030 x1.5 = 13.7%",
           vjust = 3.1, hjust = -0.05, size = 3, color = "grey40") +
  annotate("text", x = s1 * 100, y = Inf, label = "Scenario 1: universal 50%",
           vjust = 4.7, hjust = -0.05, size = 3, color = "grey40") +
  annotate("text", x = 0.98, y = 0,
           label = sprintf("Marginal gain (eta=0.3): +1pp coverage ~= %0.0f DALYs averted/yr globally", m),
           hjust = 1, vjust = -0.5, size = 3.2) +
  scale_color_manual(values = c("0.2" = "#92c5de", "0.3" = "#f4a582", "0.5" = "#ca0020"), labels = lbl_eta) +
  scale_linetype_manual(values = c("0.2" = 3, "0.3" = 1, "0.5" = 1), labels = lbl_eta) +
  labs(title = "Global: depression treatment coverage vs avoidable dementia DALYs per year",
       subtitle = "Treatment-adjusted PAF counterfactual (relative to current coverage), GBD 2023",
       x = "Target treatment coverage (%)", y = "Avoidable DALYs per year (10,000)",
       color = NULL, linetype = NULL) +
  theme_classic(base_size = 11)

## ---- Panel B: top-4 country facets, eta = 0.3 ----
top4c <- top4 %>% filter(eta == 0.3) %>%
  left_join(df %>% select(location, c_now), by = "location") %>%
  mutate(location = recode(location,
                           "United States of America" = "USA",
                           "China" = "China", "India" = "India", "Japan" = "Japan"))

pB <- ggplot(top4c, aes(c_target * 100, avoid_daly / 1e4)) +
  facet_wrap(~ location, ncol = 2, scales = "free_y") +
  annotate("rect", xmin = who2030 * 100, xmax = s1 * 100,
           ymin = -Inf, ymax = Inf, alpha = 0.12, fill = "grey50") +
  geom_line(linewidth = 1, color = "#f4a582") +
  geom_vline(aes(xintercept = c_now * 100), linetype = 3, color = "grey40") +
  geom_vline(xintercept = who2030 * 100, linetype = 4, color = "grey40") +
  geom_vline(xintercept = s1 * 100, linetype = 2, color = "grey40") +
  geom_text(data = top4c %>% group_by(location, c_now) %>% summarise(.groups = "drop"),
            aes(x = c_now * 100, y = Inf, label = sprintf("Current %.1f%%", c_now * 100)),
            vjust = 1.4, hjust = -0.05, size = 2.8, color = "grey40", inherit.aes = FALSE) +
  labs(title = "Country panels (eta = 0.3)",
       x = "Target treatment coverage (%)", y = "Avoidable DALYs per year (10,000)") +
  theme_classic(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA))

pdf("gbd_step3c_doseresp.pdf", width = 12, height = 6)
grid.arrange(pA, pB, ncol = 2, widths = c(1, 1.15))
dev.off()

cat("== dose-response v3 (locked) ==\n")
cat("Global current coverage:", c_now_global * 100, "%  WHO2030x1.5:", who2030 * 100, "%\n")
cat("Marginal gain (eta=0.3): +1pp =", round(m), "DALYs/yr\n")
chk <- glob %>% filter(eta == 0.3, c_target == 0.50) %>% mutate(wan = round(avoid_daly / 1e4, 2))
cat("Check S1(50%, eta=0.3):", chk$wan, "wan (expect 3.57)\n")
print(glob %>% filter(eta == 0.3, c_target %in% c(round(who2030, 2), 0.50, 1.00)) %>%
        mutate(wan = round(avoid_daly / 1e4, 2)))
cat("Written: gbd_step3c_doseresp.csv / gbd_step3c_doseresp.pdf\n")
