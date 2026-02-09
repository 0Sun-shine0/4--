# =========================================
# Bootstrap 内部重抽样（B=1000）乐观校正 - 多模型版
# 作者：你
# 日期：Sys.Date()
# 依赖：survival, dplyr
# =========================================











#####重新


# =========================================================
# 只保留 Apparent C-index (95% CI) 与 Optimism-corrected C-index (95% CI)
# 并比较 M0 vs M1, M0 vs M2, M1 vs M2, M2 vs M3 的差异（ΔC）
# 依赖：survival, dplyr
# =========================================================

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
})

# ---- 参数 ----
set.seed(2025)
B <- 1000                        # bootstrap 次数
data_path <- "D:/zuxue/LGE/dataz.csv"   # <-- 按需修改

# ---- 读入与基础清洗 ----
dataz <- read.csv(data_path, fileEncoding = "GBK", header = TRUE, check.names = FALSE)
names(dataz) <- trimws(names(dataz))
dataz <- dataz %>% mutate(across(where(is.character), ~ trimws(.x)))
stopifnot(all(c("Time","Events") %in% names(dataz)))

# 统一 Events 为 0/1
if (!all(unique(dataz$Events) %in% c(0,1))) {
  dataz$Events <- as.integer(as.factor(dataz$Events)) - 1L
}
dataz <- na.omit(dataz)

# 若存在 Model：低类别数时设为因子
if ("Model" %in% names(dataz)) {
  uM <- unique(dataz$Model)
  if (is.character(dataz$Model) || length(uM) <= 12) dataz$Model <- as.factor(dataz$Model)
}

# ---- 定义四个 Cox 模型（与你现用一致）----
f0 <- as.formula(Surv(Time, Events) ~ LVEF1)
f1 <- as.formula(Surv(Time, Events) ~ GLS2)
f2 <- as.formula(Surv(Time, Events) ~ Model)
f3 <- as.formula(Surv(Time, Events) ~ Model + age + NYHA)

forms <- list(
  M0_LVEF1          = f0,
  M1_GLS2           = f1,
  M2_ModelOnly      = f2,
  M3_Model_age_NYHA = f3
)

# ---- 工具：固定分层 bootstrap 索引（可复现）----
make_boot_idx <- function(events, B, seed = 2025){
  set.seed(seed)
  idx_e <- which(events == 1)
  idx_n <- which(events == 0)
  replicate(B,
            c(sample(idx_e, length(idx_e), replace = TRUE),
              sample(idx_n, length(idx_n), replace = TRUE)),
            simplify = FALSE
  )
}

# ---- 工具：收集每个模型的两类 C 及其 CI（并保留逐轮向量用于比较）----
gather_C_for_model <- function(formula, data, boot_idx) {
  # Apparent C 在原始数据上的点估计及 95%CI
  fit_app <- coxph(formula, data = data, ties = "efron", x = TRUE, y = TRUE)
  lp_app  <- predict(fit_app, type = "lp")
  cobj    <- survival::concordance(Surv(data$Time, data$Events) ~ lp_app, reverse = TRUE)
  C_app   <- as.numeric(cobj$concordance)
  se_app  <- sqrt(as.numeric(cobj$var))
  App_CI  <- if (is.finite(se_app)) C_app + c(-1,1)*1.96*se_app else c(NA_real_, NA_real_)
  
  # Bootstrap：逐轮 C_in（自助样本内）与 C_out（原始数据上）
  B <- length(boot_idx)
  C_in  <- numeric(B)
  C_out <- numeric(B)
  
  for (b in seq_len(B)) {
    ix <- boot_idx[[b]]
    dboot <- data[ix, , drop = FALSE]
    fit_b <- tryCatch(coxph(formula, data = dboot, ties="efron", x=TRUE, y=TRUE),
                      error = function(e) NULL)
    if (is.null(fit_b)) { C_in[b] <- NA; C_out[b] <- NA; next }
    lp_in  <- predict(fit_b, newdata = dboot, type = "lp")
    lp_out <- predict(fit_b, newdata = data,  type = "lp")
    C_in[b]  <- survival::concordance(Surv(dboot$Time, dboot$Events) ~ lp_in,  reverse=TRUE)$concordance
    C_out[b] <- survival::concordance(Surv(data$Time,  data$Events)  ~ lp_out, reverse=TRUE)$concordance
  }
  
  ok <- is.finite(C_in) & is.finite(C_out)
  C_in <- C_in[ok]; C_out <- C_out[ok]
  
  # 乐观校正后 C 的逐轮分布与区间（percentile）
  C_corr_vec <- C_app - (C_in - C_out)         # 每轮的“校正后 C”
  C_corr     <- mean(C_corr_vec)               # 点估计
  C_corr_CI  <- quantile(C_corr_vec, c(.025,.975), na.rm = TRUE)
  
  list(
    C_app = C_app, App_CI = as.numeric(App_CI), C_in = C_in,
    C_corr = C_corr, C_corr_CI = as.numeric(C_corr_CI), C_corr_vec = C_corr_vec
  )
}

# ---- 运行：固定索引 + 收集各模型指标 ----
boot_idx_fixed <- make_boot_idx(dataz$Events, B = B, seed = 2025)
Cstore <- lapply(forms, function(f) gather_C_for_model(f, dataz, boot_idx_fixed))
names(Cstore) <- names(forms)

# ---- 汇总：仅两类 C（点估计 + 95%CI），便于论文表格 ----
fmt_ci <- function(est, lo, hi, digits = 4){
  paste0(sprintf(paste0("%.",digits,"f"), est),
         " (", sprintf(paste0("%.",digits,"f"), lo),
         "–",  sprintf(paste0("%.",digits,"f"), hi), ")")
}
res_twoC <- do.call(rbind, lapply(names(Cstore), function(nm){
  x <- Cstore[[nm]]
  data.frame(
    Model = nm,
    Apparent_C_index_95CI  = fmt_ci(x$C_app, x$App_CI[1],   x$App_CI[2],   digits = 4),
    Corrected_C_index_95CI = fmt_ci(x$C_corr, x$C_corr_CI[1], x$C_corr_CI[2], digits = 4),
    stringsAsFactors = FALSE
  )
}))
# 按校正后 C 从高到低排序
ord <- order(sapply(Cstore, function(x) x$C_corr), decreasing = TRUE)
res_twoC <- res_twoC[ord, ]

cat("\n===== 两类 C-index（点估计 + 95%CI） =====\n")
print(res_twoC, row.names = FALSE)

# ---- 成对比较：ΔC 的均值、95%CI（percentile）、双侧 p 值 + Holm 校正 ----
pairwise_delta <- function(Cstore, pairs, which = c("apparent","corrected")){
  which <- match.arg(which)
  out <- lapply(pairs, function(pr){
    a <- pr[1]; b <- pr[2]
    if (which == "apparent") { va <- Cstore[[a]]$C_in;       vb <- Cstore[[b]]$C_in }
    else                     { va <- Cstore[[a]]$C_corr_vec;  vb <- Cstore[[b]]$C_corr_vec }
    k <- min(length(va), length(vb))
    d <- (va[seq_len(k)] - vb[seq_len(k)])
    d <- d[is.finite(d)]
    est <- mean(d)
    ci  <- quantile(d, c(.025,.975), na.rm = TRUE)
    p   <- 2 * min(mean(d <= 0), mean(d >= 0))   # 成对“符号检验”近似
    data.frame(
      contrast = paste0(a, " - ", b),
      mean_delta = round(est, 4),
      CI_L = round(ci[1], 4),
      CI_U = round(ci[2], 4),
      p_value = signif(p, 3),
      n_effective = length(d),
      prop_delta_pos = round(mean(d > 0), 3),
      stringsAsFactors = FALSE
    )
  })
  res <- do.call(rbind, out)
  res$p_adj_holm <- p.adjust(res$p_value, method = "holm")
  res
}

pairs_to_test <- list(
  c("M0_LVEF1", "M1_GLS2"),
  c("M0_LVEF1", "M2_ModelOnly"),
  c("M1_GLS2",  "M2_ModelOnly"),
  c("M2_ModelOnly", "M3_Model_age_NYHA")
)

# 1) ΔC_app（Apparent C 的差异）
delta_app <- pairwise_delta(Cstore, pairs_to_test, which = "apparent")
cat("\n===== Paired Bootstrap ΔC_app（Apparent C 的差异） =====\n")
print(delta_app, row.names = FALSE)

# 2) ΔC_corr（Optimism-corrected C 的差异）
delta_corr <- pairwise_delta(Cstore, pairs_to_test, which = "corrected")
cat("\n===== Paired Bootstrap ΔC_corr（Optimism-corrected C 的差异） =====\n")
print(delta_corr, row.names = FALSE)

# ---- 可选导出 ----
# write.csv(res_twoC,   "Cindex_two_metrics.csv", row.names = FALSE)
# write.csv(delta_app,  "DeltaC_app_pairwise.csv", row.names = FALSE)
# write.csv(delta_corr, "DeltaC_corrected_pairwise.csv", row.names = FALSE)











# =========================================================
# C-index 与 Calibration slope（含95%CI）的 Apparent 与 Optimism-corrected
# 依赖：survival, dplyr
# 假设你已有 dataz, f0,f1,f2,f3（与之前一致）
# =========================================================
suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
})

# ---- 工具 ----
harrell_c <- function(time, event, lp) {
  survival::concordance(Surv(time, event) ~ lp, reverse = TRUE)$concordance
}

make_boot_idx <- function(events, B, seed = 2025){
  set.seed(seed)
  idx_e <- which(events == 1)
  idx_n <- which(events == 0)
  replicate(B,
    c(sample(idx_e, length(idx_e), replace = TRUE),
      sample(idx_n, length(idx_n), replace = TRUE)),
    simplify = FALSE
  )
}

# 核心函数：同时计算 C 与 slope（含显著性区间）
boot_perf_with_slopes <- function(formula, data, boot_idx, show_progress = TRUE){
  # 1) 表观拟合
  fit_app <- coxph(formula, data = data, ties = "efron", x = TRUE, y = TRUE)
  lp_app  <- predict(fit_app, type = "lp")

  # 1a) Apparent C 与 95%CI（基于 concordance 方差）
  cobj    <- survival::concordance(Surv(data$Time, data$Events) ~ lp_app, reverse = TRUE)
  C_app   <- as.numeric(cobj$concordance)
  se_appC <- sqrt(as.numeric(cobj$var))
  AppC_CI <- if (is.finite(se_appC)) C_app + c(-1,1)*1.96*se_appC else c(NA_real_, NA_real_)

  # 1b) Apparent slope 与 95%CI（直接回归 Surv ~ lp_app）
  slope_app_fit <- coxph(Surv(Time, Events) ~ lp_app, data = data, ties = "efron")
  slope_app     <- coef(slope_app_fit)[1]
  se_appS       <- sqrt(vcov(slope_app_fit)[1,1])
  AppSlope_CI   <- slope_app + c(-1,1)*1.96*se_appS

  # 2) Bootstrap：逐轮统计
  B <- length(boot_idx)
  C_in <- C_out <- numeric(B)
  slope_in <- slope_out <- numeric(B)

  if (show_progress) pb <- txtProgressBar(min = 0, max = B, style = 3)
  for (b in seq_len(B)) {
    ix <- boot_idx[[b]]
    dboot <- data[ix, , drop = FALSE]

    fit_b <- tryCatch(coxph(formula, data = dboot, ties="efron", x=TRUE, y=TRUE),
                      error = function(e) NULL)
    if (is.null(fit_b)) {
      C_in[b] <- NA; C_out[b] <- NA; slope_in[b] <- NA; slope_out[b] <- NA
      if (show_progress) setTxtProgressBar(pb, b)
      next
    }

    # 线性预测子
    lp_in  <- predict(fit_b, newdata = dboot, type = "lp")
    lp_out <- predict(fit_b, newdata = data,  type = "lp")

    # C_in / C_out
    C_in[b]  <- harrell_c(dboot$Time, dboot$Events, lp_in)
    C_out[b] <- harrell_c(data$Time,  data$Events,  lp_out)

    # slope_in：在自助样本内以 lp_in 作为唯一自变量（理论上 ~1）
    slope_in_fit <- tryCatch(coxph(Surv(Time, Events) ~ lp_in, data = dboot, ties="efron"),
                             error=function(e) NULL)
    slope_in[b]  <- if (!is.null(slope_in_fit)) coef(slope_in_fit)[1] else NA_real_

    # slope_out：将自助拟合模型的 LP 应用到原始数据
    slope_out_fit <- tryCatch(coxph(Surv(Time, Events) ~ lp_out, data = data, ties="efron"),
                              error=function(e) NULL)
    slope_out[b] <- if (!is.null(slope_out_fit)) coef(slope_out_fit)[1] else NA_real_

    if (show_progress) setTxtProgressBar(pb, b)
  }
  if (show_progress) close(pb)

  ok <- is.finite(C_in) & is.finite(C_out) & is.finite(slope_in) & is.finite(slope_out)
  C_in <- C_in[ok]; C_out <- C_out[ok]
  slope_in <- slope_in[ok]; slope_out <- slope_out[ok]

  # 3) 乐观校正（C 与 slope 的“逐轮校正值”与区间）
  C_corr_vec     <- C_app - (C_in - C_out)
  C_corr         <- mean(C_corr_vec)
  C_corr_CI      <- quantile(C_corr_vec, c(.025,.975), na.rm = TRUE)

  slope_corr_vec <- slope_app - (slope_in - slope_out)  # Apparent slope - optimism
  slope_corr     <- mean(slope_corr_vec)
  slope_corr_CI  <- quantile(slope_corr_vec, c(.025,.975), na.rm = TRUE)

  # 4) 同时给出 Validation slope（= slope_out 的均值与CI，常作为“统一收缩”参考）
  val_slope      <- mean(slope_out)
  ValSlope_CI    <- quantile(slope_out, c(.025,.975), na.rm = TRUE)

  # 5) 也给出“Bootstrap-validated C”的分布（可选）
  BootC          <- mean(C_out)
  BootC_CI       <- quantile(C_out, c(.025,.975), na.rm = TRUE)

  list(
    # C
    Apparent_C = C_app, Apparent_CI = as.numeric(AppC_CI),
    Corrected_C = C_corr, Corrected_CI = as.numeric(C_corr_CI),
    BootValidated_C = BootC, BootValidated_CI = as.numeric(BootC_CI),
    # Slopes
    Apparent_Slope = slope_app, Apparent_Slope_CI = as.numeric(AppSlope_CI),
    Corrected_Slope = slope_corr, Corrected_Slope_CI = as.numeric(slope_corr_CI),
    Val_Slope = val_slope, Val_Slope_CI = as.numeric(ValSlope_CI),
    # 便于后续需要：逐轮向量
    C_in = C_in, C_out = C_out, slope_in = slope_in, slope_out = slope_out
  )
}

# ---- 运行：你的四个模型 ----
# 若尚未定义四个公式，按你之前设定再定义一次
if (!exists("f0")) f0 <- as.formula(Surv(Time, Events) ~ LVEF1)
if (!exists("f1")) f1 <- as.formula(Surv(Time, Events) ~ GLS2)
if (!exists("f2")) f2 <- as.formula(Surv(Time, Events) ~ Model)
if (!exists("f3")) f3 <- as.formula(Surv(Time, Events) ~ Model + age + NYHA)

forms <- list(
  `Model 0 (LVEF)`                 = f0,
  `Model 1 (GLS)`                  = f1,
  `Model 2 (GLS + LVEF)`           = f2,
  `Model 3 (GLS + LVEF + Age + NYHA)` = f3
)

B <- if (exists("B")) B else 1000
boot_idx_fixed <- make_boot_idx(dataz$Events, B = B, seed = 2025)

res_list2 <- lapply(forms, function(f) boot_perf_with_slopes(f, dataz, boot_idx_fixed, show_progress = TRUE))
names(res_list2) <- names(forms)

# ---- 汇总成表（含 C 与 slope 的 apparent & corrected 两种口径）----
fmt_ci <- function(est, ci, digits=3){
  paste0(sprintf(paste0("%.",digits,"f"), est),
         " (", sprintf(paste0("%.",digits,"f"), ci[1]),
         "–",  sprintf(paste0("%.",digits,"f"), ci[2]), ")")
}

tab_out <- do.call(rbind, lapply(names(res_list2), function(nm){
  x <- res_list2[[nm]]
  data.frame(
    Model = nm,
    Apparent_C_index_95CI          = fmt_ci(x$Apparent_C,      x$Apparent_CI,      3),
    Optimism_corrected_C_index_95CI= fmt_ci(x$Corrected_C,     x$Corrected_CI,     3),
    Apparent_Calibration_Slope_95CI= fmt_ci(x$Apparent_Slope,  x$Apparent_Slope_CI,3),
    Optimism_corrected_Slope_95CI  = fmt_ci(x$Corrected_Slope, x$Corrected_Slope_CI,3),
    Validation_Slope_95CI          = fmt_ci(x$Val_Slope,       x$Val_Slope_CI,     3),
    stringsAsFactors = FALSE
  )
}))
print(tab_out, row.names = FALSE)

# 可选导出
# write.csv(tab_out, "Cindex_and_Slope_Apparent_vs_Corrected.csv", row.names = FALSE)


















###比较C指数


# =========================================
# Pairwise Bootstrap Comparison of C-index Differences
# 修正版：同时输出 Apparent & Optimism-corrected，顺序正确
# 依赖：survival, dplyr
# =========================================
suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
})

# ---------- 配置 ----------
set.seed(2025)
B <- 1000
SHOW_ABS_DIFF <- TRUE   # 表中 ΔC 是否显示为绝对值；若想看带符号差，改为 FALSE

# ---------- 读入数据（若已存在 dataz 则跳过） ----------
if (!exists("dataz")) {
  data_path <- "D:/zuxue/LGE/dataz.csv"   # <-- 按需修改
  dataz <- read.csv(data_path, fileEncoding = "GBK", header = TRUE, check.names = FALSE)
  names(dataz) <- trimws(names(dataz))
  dataz <- dataz %>% mutate(across(where(is.character), ~ trimws(.x)))
  stopifnot(all(c("Time","Events") %in% names(dataz)))
  if (!all(unique(dataz$Events) %in% c(0,1))) {
    dataz$Events <- as.integer(as.factor(dataz$Events)) - 1L
  }
  dataz <- na.omit(dataz)
  if ("Model" %in% names(dataz)) {
    uM <- unique(dataz$Model)
    if (is.character(dataz$Model) || length(uM) <= 12) dataz$Model <- as.factor(dataz$Model)
  }
}

# ---------- 定义四个 Cox 模型（若已存在 f0~f3 则使用已定义者） ----------
if (!exists("f0")) f0 <- as.formula(Surv(Time, Events) ~ LVEF1)                    # Model 0
if (!exists("f1")) f1 <- as.formula(Surv(Time, Events) ~ GLS2)                     # Model 1
if (!exists("f2")) f2 <- as.formula(Surv(Time, Events) ~ Model)                    # Model 2
if (!exists("f3")) f3 <- as.formula(Surv(Time, Events) ~ Model + age + NYHA)       # Model 3

forms <- list(
  `Model 0` = f0,
  `Model 1` = f1,
  `Model 2` = f2,
  `Model 3` = f3
)

# ---------- 小工具 ----------
harrell_c <- function(time, event, lp) {
  survival::concordance(Surv(time, event) ~ lp, reverse = TRUE)$concordance
}

# 固定分层 bootstrap 索引（保持事件比例；可复现）
make_boot_idx <- function(events, B = 1000, seed = 2025){
  set.seed(seed)
  idx_e <- which(events == 1)
  idx_n <- which(events == 0)
  replicate(B,
            c(sample(idx_e, length(idx_e), replace = TRUE),
              sample(idx_n, length(idx_n), replace = TRUE)),
            simplify = FALSE)
}

# 收集各模型：C_in、C_out、Apparent C
gather_C_for_model <- function(formula, data, boot_idx) {
  fit_app <- coxph(formula, data = data, ties="efron", x=TRUE, y=TRUE)
  lp_app  <- predict(fit_app, type="lp")
  C_app   <- harrell_c(data$Time, data$Events, lp_app)
  
  B  <- length(boot_idx)
  C_in  <- numeric(B)
  C_out <- numeric(B)
  
  for (b in seq_len(B)) {
    ix <- boot_idx[[b]]
    dboot <- data[ix, , drop = FALSE]
    fit_b <- tryCatch(coxph(formula, data = dboot, ties="efron", x=TRUE, y=TRUE),
                      error=function(e) NULL)
    if (is.null(fit_b)) { C_in[b] <- NA; C_out[b] <- NA; next }
    lp_in  <- predict(fit_b, newdata=dboot, type="lp")
    lp_out <- predict(fit_b, newdata=data,  type="lp")
    C_in[b]  <- harrell_c(dboot$Time, dboot$Events, lp_in)
    C_out[b] <- harrell_c(data$Time,  data$Events,  lp_out)
  }
  ok <- is.finite(C_in) & is.finite(C_out)
  C_in <- C_in[ok]; C_out <- C_out[ok]
  
  # 每轮“乐观校正后 C”：C_app - (C_in - C_out)
  C_corr_vec <- C_app - (C_in - C_out)
  
  list(C_in = C_in, C_corr_vec = C_corr_vec)
}

# 成对比较：输出 mean ΔC、95%CI、p、Holm p（表观 or 校正）
pairwise_delta <- function(Cstore, pairs, metric = c("apparent","corrected"),
                           show_abs = TRUE) {
  metric <- match.arg(metric)
  rows <- lapply(pairs, function(pr){
    A <- pr[1]; B <- pr[2]
    va <- if (metric=="apparent") Cstore[[A]]$C_in else Cstore[[A]]$C_corr_vec
    vb <- if (metric=="apparent") Cstore[[B]]$C_in else Cstore[[B]]$C_corr_vec
    k  <- min(length(va), length(vb))
    d  <- va[seq_len(k)] - vb[seq_len(k)]            # 带符号差（用于 p 值）
    d  <- d[is.finite(d)]
    
    # 展示用：绝对差或带符号差
    d_show <- if (show_abs) abs(d) else d
    
    est <- mean(d_show)
    ci  <- quantile(d_show, c(.025,.975), na.rm = TRUE)
    # p 值（双侧符号检验）基于“带符号”的 d
    p   <- 2 * min(mean(d <= 0), mean(d >= 0))
    
    data.frame(
      Contrast = paste0(A, " vs. ", B),
      Metric   = if (metric=="apparent") "Apparent C-index" else "Optimism-corrected C-index",
      mean_delta = est, CI_L = ci[1], CI_U = ci[2],
      p_value = p, n_effective = length(d),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$Holm_adjusted_p <- p.adjust(out$p_value, method = "holm")
  out
}

# ---------- 运行 ----------
boot_idx <- make_boot_idx(dataz$Events, B = B, seed = 2025)
Cstore <- lapply(forms, function(f) gather_C_for_model(f, dataz, boot_idx))
names(Cstore) <- names(forms)

# 需要的对比顺序（与你表一致）
pairs_to_test <- list(
  c("Model 1","Model 0"),
  c("Model 0","Model 2"),
  c("Model 1","Model 2"),
  c("Model 2","Model 3")
)

# 两类指标分别做成对比较
tab_app  <- pairwise_delta(Cstore, pairs_to_test, metric = "apparent",  show_abs = SHOW_ABS_DIFF)
tab_corr <- pairwise_delta(Cstore, pairs_to_test, metric = "corrected", show_abs = SHOW_ABS_DIFF)

# 合并为最终长表，并格式化 + 稳健排序（每个对比下：Apparent 在前，Corrected 在后）
contrast_levels <- c(
  "Model 1 vs. Model 0",
  "Model 0 vs. Model 2",
  "Model 1 vs. Model 2",
  "Model 2 vs. Model 3"
)
metric_levels <- c("Apparent C-index", "Optimism-corrected C-index")

tab_all <- bind_rows(tab_app, tab_corr) %>%
  mutate(
    Contrast = factor(Contrast, levels = contrast_levels),
    Metric   = factor(Metric,   levels = metric_levels)
  ) %>%
  arrange(Contrast, Metric) %>%
  mutate(
    `mean ΔC` = sprintf("%.3f", mean_delta),
    `95% CI`  = paste0("(", sprintf("%.3f", CI_L), " to ", sprintf("%.3f", CI_U), ")"),
    `p-value` = ifelse(p_value < 0.001, "< 0.001", sprintf("%.3f", p_value)),
    `Holm-adjusted p` = ifelse(Holm_adjusted_p < 0.001, "< 0.001", sprintf("%.3f", Holm_adjusted_p))
  ) %>%
  select(Contrast, Metric, `mean ΔC`, `95% CI`, `p-value`, `Holm-adjusted p`)

cat("\n===== Table 2. Pairwise Bootstrap Comparison of Model Performance Differences =====\n")
print(tab_all, row.names = FALSE)

# 可选导出：
# write.csv(tab_all, "Table2_Pairwise_Cindex_Differences.csv", row.names = FALSE)

















