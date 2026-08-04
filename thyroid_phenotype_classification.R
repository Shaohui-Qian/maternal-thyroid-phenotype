# ==============================================================================
# 核心流程：
#   1. 多源数据表关联合并
#   2. TSH/FT4/TPOAb 按孕期分层判定
#   3. 基础单病编码 + 单人多条记录聚合 + 0编码优先级过滤
#   4. 共病标准化重编码（结节/癌优先、两病专属码、三病归8、四病归9）
#   5. 结果导出 + 异常病例质控核查
library(tidyverse)
library(readxl)
library(writexl)
# FT4 不同孕期参考区间
ft4_cutoff <- list(
  early      = list(lower = 11.8, upper = 21.0),   # 孕早
  mid        = list(lower = 10.6, upper = 17.6),   # 孕中
  late       = list(lower = 9.2,  upper = 16.7),   # 孕晚
  postpartum = list(lower = 10,    upper = 31)     # 产后
)

# TSH 不同孕期参考区间
tsh_cutoff <- list(
  early      = list(lower = 0.03, upper = 4.51),
  mid        = list(lower = 0.05, upper = 4.50),
  late       = list(lower = 0.47, upper = 4.54),
  postpartum = list(lower = 0.3,  upper = 4.2)
)

# TPOAb 阳性阈值
tpoab_threshold <- 60


#' FT4 按孕期分层判定
classify_ft4 <- function(gestation_stage, ft4_value) {
  case_when(
    gestation_stage == "孕早" & ft4_value < ft4_cutoff$early$lower      ~ "降低",
    gestation_stage == "孕早" & between(ft4_value, ft4_cutoff$early$lower, ft4_cutoff$early$upper) ~ "正常",
    gestation_stage == "孕早" & ft4_value > ft4_cutoff$early$upper      ~ "偏高",
    
    gestation_stage == "孕中" & ft4_value < ft4_cutoff$mid$lower        ~ "降低",
    gestation_stage == "孕中" & between(ft4_value, ft4_cutoff$mid$lower, ft4_cutoff$mid$upper) ~ "正常",
    gestation_stage == "孕中" & ft4_value > ft4_cutoff$mid$upper        ~ "偏高",
    
    gestation_stage == "孕晚" & ft4_value < ft4_cutoff$late$lower       ~ "降低",
    gestation_stage == "孕晚" & between(ft4_value, ft4_cutoff$late$lower, ft4_cutoff$late$upper) ~ "正常",
    gestation_stage == "孕晚" & ft4_value > ft4_cutoff$late$upper       ~ "偏高",
    
    gestation_stage == "产后" & ft4_value < ft4_cutoff$postpartum$lower ~ "降低",
    gestation_stage == "产后" & between(ft4_value, ft4_cutoff$postpartum$lower, ft4_cutoff$postpartum$upper) ~ "正常",
    gestation_stage == "产后" & ft4_value > ft4_cutoff$postpartum$upper ~ "偏高",
    
    is.na(ft4_value) ~ "缺失",
    TRUE ~ NA_character_
  )
}

#' TSH 按孕期分层判定
classify_tsh <- function(gestation_stage, tsh_value) {
  case_when(
    gestation_stage == "孕早" & tsh_value < tsh_cutoff$early$lower      ~ "降低",
    gestation_stage == "孕早" & between(tsh_value, tsh_cutoff$early$lower, tsh_cutoff$early$upper) ~ "正常",
    gestation_stage == "孕早" & tsh_value > tsh_cutoff$early$upper      ~ "偏高",
    gestation_stage == "孕早" & is.na(tsh_value)                        ~ "缺失",
    
    gestation_stage == "孕中" & tsh_value < tsh_cutoff$mid$lower        ~ "降低",
    gestation_stage == "孕中" & between(tsh_value, tsh_cutoff$mid$lower, tsh_cutoff$mid$upper) ~ "正常",
    gestation_stage == "孕中" & tsh_value > tsh_cutoff$mid$upper        ~ "偏高",
    gestation_stage == "孕中" & is.na(tsh_value)                        ~ "缺失",
    
    gestation_stage == "孕晚" & tsh_value < tsh_cutoff$late$lower       ~ "降低",
    gestation_stage == "孕晚" & between(tsh_value, tsh_cutoff$late$lower, tsh_cutoff$late$upper) ~ "正常",
    gestation_stage == "孕晚" & tsh_value > tsh_cutoff$late$upper       ~ "偏高",
    gestation_stage == "孕晚" & is.na(tsh_value)                        ~ "缺失",
    
    gestation_stage == "产后" & tsh_value < tsh_cutoff$postpartum$lower ~ "降低",
    gestation_stage == "产后" & between(tsh_value, tsh_cutoff$postpartum$lower, tsh_cutoff$postpartum$upper) ~ "正常",
    gestation_stage == "产后" & tsh_value > tsh_cutoff$postpartum$upper ~ "偏高",
    gestation_stage == "产后" & is.na(tsh_value)                        ~ "缺失",
    
    TRUE ~ NA_character_
  )
}

#' 0 编码优先级过滤（良性编码不参与疾病共病统计）
#'
#' 规则：
#'   - 全部为 0X 良性编码 → 有 00 取 00，否则取第一个
#'   - 混杂 0X 和真实疾病编码 → 剔除全部 0 开头编码，只保留疾病编码
filter_zero_code <- function(code_str) {
  if (is.na(code_str) || code_str == "") return(code_str)
  
  codes <- unlist(str_split(code_str, "&"))
  all_zero <- all(str_detect(codes, "^0[0-9]"))
  
  if (!all_zero) {
    severe <- codes[!str_detect(codes, "^0")]
    return(paste(severe, collapse = "&"))
  }
  
  if ("00" %in% codes) return("00")
  return(codes[1])
}

#' 共病标准化重编码（表三核心规则）
#' 规则优先级从高到低：
#'   1. 含结节(6)/癌(7) → 优先归入 6 组或 7 组，分配专属共病码
#'   2. 两病组合 → 分配固定专属编码
#'   3. 三种疾病 → 统一归 8 组
#'   4. 四种及以上 → 统一归 9 组
#'   5. 单一疾病 → 直接返回原编码
map_comorbidity_code <- function(code_str) {
  if (is.na(code_str) || code_str == "分类不明") return(code_str)
  
  codes <- unlist(str_split(code_str, "&"))
  categories <- sort(unique(substr(codes, 1, 1)))
  n_disease <- length(categories)
  cat_key <- paste(categories, collapse = "&")
  
  # ---------- 优先级1：含甲状腺结节(6) → 优先归 6 组 ----------
  if ("6" %in% categories) {
    nodule_map <- list(
      "1&6"       = "61",   # 桥本 + 结节
      "2&6"       = "62",   # 甲减 + 结节
      "3&6"       = "63",   # 甲状腺毒症 + 结节
      "4&6"       = "64",   # 甲亢 + 结节
      "5&6"       = "65",   # 低甲状腺素血症 + 结节
      "1&2&6"     = "67",   # 桥本 + 甲减 + 结节
      "1&3&6"     = "68",   # 桥本 + 毒症 + 结节
      "1&4&6"     = "69",   # 桥本 + 甲亢 + 结节
      "1&5&6"     = "691",  # 桥本 + 低甲素 + 结节
      "2&3&6"     = "692",  # 甲减 + 毒症 + 结节
      "2&5&6"     = "693",  # 甲减 + 低甲素 + 结节
      "3&4&6"     = "694",  # 毒症 + 甲亢 + 结节
      "3&5&6"     = "692"   # 毒症 + 低甲素 + 结节
    )
    if (!is.null(nodule_map[[cat_key]])) return(nodule_map[[cat_key]])
    return("60")
  }
  
  # ---------- 优先级2：含甲状腺癌(7) → 优先归 7 组 ----------
  if ("7" %in% categories) {
    cancer_map <- list(
      "1&7"       = "71",   # 桥本 + 癌
      "2&7"       = "72",   # 甲减 + 癌
      "3&7"       = "73",   # 毒症 + 癌
      "5&7"       = "74",   # 低甲素 + 癌
      "6&7"       = "75",   # 结节 + 癌
      "1&3&7"     = "76",   # 桥本 + 毒症 + 癌
      "1&5&7"     = "77",   # 桥本 + 低甲素 + 癌
      "1&6&7"     = "78",   # 桥本 + 结节 + 癌
      "2&3&7"     = "79",   # 甲减 + 毒症 + 癌
      "2&5&7"     = "791",  # 甲减 + 低甲素 + 癌
      "3&4&7"     = "792",  # 毒症 + 甲亢 + 癌
      "3&5&7"     = "793"   # 毒症 + 低甲素 + 癌
    )
    if (!is.null(cancer_map[[cat_key]])) return(cancer_map[[cat_key]])
    return("70")
  }
  
  # ---------- 优先级3：两病组合 → 专属固定编码 ----------
  pair_map <- list(
    "1&2" = "13",   # 桥本 + 甲减
    "1&3" = "14",   # 桥本 + 甲状腺毒症
    "1&4" = "15",   # 桥本 + 甲亢
    "1&5" = "16",   # 桥本 + 低甲状腺素血症
    "2&3" = "26",   # 甲减 + 甲状腺毒症
    "2&5" = "27",   # 甲减 + 低甲状腺素血症
    "3&4" = "33",   # 甲状腺毒症 + 甲亢
    "3&5" = "34",   # 甲状腺毒症 + 低甲状腺素血症
    "4&5" = "43"    # 甲亢 + 低甲状腺素血症
  )
  if (n_disease == 2 && !is.null(pair_map[[cat_key]])) {
    return(pair_map[[cat_key]])
  }
  
  # ---------- 优先级4：三种疾病 → 统一归 8 组 ----------
  if (n_disease == 3) {
    triple_map <- list(
      "1&2&3" = "80",   # 桥本 + 甲减 + 毒症
      "1&2&5" = "81",   # 桥本 + 甲减 + 低甲素
      "1&3&4" = "82",   # 桥本 + 毒症 + 甲亢
      "1&3&5" = "83",   # 桥本 + 毒症 + 低甲素
      "1&4&5" = "84",   # 桥本 + 甲亢 + 低甲素
      "2&3&5" = "85",   # 甲减 + 毒症 + 低甲素
      "3&4&5" = "86"    # 毒症 + 甲亢 + 低甲素
    )
    if (!is.null(triple_map[[cat_key]])) return(triple_map[[cat_key]])
    return("80")
  }
  
  # ---------- 优先级5：四种及以上 → 统一归 9 组 ----------
  if (n_disease >= 4) {
    quadruple_map <- list(
      "1&2&3&5"   = "90",   # 桥本+甲减+毒症+低甲素
      "1&2&5&6"   = "91",   # 桥本+甲减+低甲素+结节
      "1&3&4&5"   = "92",   # 桥本+毒症+甲亢+低甲素
      "1&3&4&6"   = "93",   # 桥本+毒症+甲亢+结节
      "1&3&5&6"   = "94",   # 桥本+毒症+低甲素+结节
      "3&4&5&6"   = "95",   # 毒症+甲亢+低甲素+结节
      "3&4&6&7"   = "96",   # 毒症+甲亢+结节+癌
      "1&3&4&5&6" = "97"    # 五种病
    )
    if (!is.null(quadruple_map[[cat_key]])) return(quadruple_map[[cat_key]])
    return("90")
  }
  
  # ---------- 单一疾病 → 直接返回 ----------
  return(codes[1])
}

#' 根据精细编码提取疾病大类（首数字）
extract_disease_category <- function(code_str) {
  if (is.na(code_str) || code_str == "") return(code_str)
  if (code_str == "分类不明") return("分类不明")
  
  codes <- unlist(str_split(code_str, "&"))
  categories <- unique(substr(codes, 1, 1))
  paste(categories, collapse = "&")
}
# 步骤1：数据读取与合并
# 读取原始数据表
df_sample_id <- read_xlsx("0_标准样本ID.xlsx")
df_diagnosis <- read_xlsx("6_诊断记录.xlsx")
df_lab_other <- read_xlsx("4_检验-其他.xlsx")
# 诊断表 + 检验表关联合并（同一人允许多条检测记录）
df_raw <- left_join(
  df_diagnosis,
  df_lab_other,
  by = "NIPTID",
  relationship = "many-to-many"
)
# 步骤2：实验室指标按孕期分层判定
df_stratified <- df_raw %>%
  mutate(
    TPOAb_group = case_when(
      `甲状腺过氧化物酶抗体TPOA` <= tpoab_threshold ~ "正常",
      `甲状腺过氧化物酶抗体TPOA` >  tpoab_threshold ~ "偏高",
      is.na(`甲状腺过氧化物酶抗体TPOA`)             ~ "缺失",
      TRUE ~ NA_character_
    ),
    FT4_group = classify_ft4(`分期.x`, `血清游离甲状腺素FT4`),
    TSH_group = classify_tsh(`分期.x`, `促甲状腺激素TSH`)
  )

# 步骤3：基础单病编码
df_basic_coded <- df_stratified %>%
  mutate(
    disease_code = case_when(
      # ---- 0 组：非甲功异常 ----
      TSH_group == "正常" & FT4_group == "正常" & TPOAb_group == "正常" &
        (is.na(记录内容) | !str_detect(记录内容, "甲|桥本|甲状腺|Graves")) ~ "00",
      
      TSH_group == "正常" & FT4_group == "正常" & TPOAb_group == "缺失" &
        (is.na(记录内容) | !str_detect(记录内容, "甲|桥本|甲状腺|Graves")) ~ "01",
      
      TSH_group == "正常" & FT4_group == "正常" & TPOAb_group == "正常" &
        str_detect(记录内容, "甲") ~ "08",
      
      ((TSH_group == "缺失" | FT4_group == "缺失") & TPOAb_group == "缺失") &
        (is.na(记录内容) | !str_detect(记录内容, "甲|Graves")) ~ "09",
      
      # ---- 1 组：桥本氏甲状腺炎 ----
      TPOAb_group == "偏高" & str_detect(记录内容, "桥本") ~ "10",
      TPOAb_group == "偏高" & (is.na(记录内容) | !str_detect(记录内容, "桥本")) ~ "11",
      TPOAb_group %in% c("正常", "缺失") & str_detect(记录内容, "桥本") ~ "12",
      
      # ---- 2 组：甲状腺功能减退 ----
      TSH_group == "偏高" & FT4_group == "降低" & str_detect(记录内容, "甲状腺功能低下|甲减") ~ "20",
      TSH_group == "偏高" & FT4_group == "降低" & (is.na(记录内容) | !str_detect(记录内容, "甲状腺功能低下|甲减")) ~ "21",
      (TSH_group == "缺失" | FT4_group == "缺失") & str_detect(记录内容, "甲状腺功能低下|甲减") ~ "22",
      
      # ---- 2 组：亚临床甲减 ----
      TSH_group == "偏高" & FT4_group == "正常" & str_detect(记录内容, "亚临床甲减|亚临床甲状腺功能减退") ~ "23",
      TSH_group == "偏高" & FT4_group == "正常" & (is.na(记录内容) | !str_detect(记录内容, "亚临床甲减|亚临床甲状腺功能减退")) ~ "24",
      (TSH_group == "缺失" | FT4_group == "缺失") & str_detect(记录内容, "亚临床甲减|亚临床甲状腺功能减退") ~ "25",
      
      # ---- 3 组：妊娠期甲状腺毒症 ----
      TSH_group == "降低" & FT4_group == "正常" & str_detect(记录内容, "甲状腺毒症") ~ "30",
      TSH_group == "降低" & FT4_group == "正常" & (is.na(记录内容) | !str_detect(记录内容, "甲状腺毒症")) ~ "31",
      (TSH_group == "缺失" | FT4_group == "缺失") & str_detect(记录内容, "甲状腺毒症") ~ "32",
      
      # ---- 4 组：甲状腺功能亢进 ----
      TSH_group == "降低" & FT4_group == "偏高" & str_detect(记录内容, "甲状腺功能亢进|甲亢|Graves病") ~ "40",
      TSH_group == "降低" & FT4_group == "偏高" & (is.na(记录内容) | !str_detect(记录内容, "甲状腺功能亢进|甲亢|Graves病")) ~ "41",
      (TSH_group == "缺失" | FT4_group == "缺失") & str_detect(记录内容, "甲状腺功能亢进|甲亢|Graves病") ~ "42",
      
      # ---- 5 组：低甲状腺素血症 ----
      FT4_group == "降低" & TSH_group == "正常" & str_detect(记录内容, "低甲状腺素血症") ~ "50",
      FT4_group == "降低" & TSH_group == "正常" & (is.na(记录内容) | !str_detect(记录内容, "低甲状腺素血症")) ~ "51",
      (FT4_group %in% c("正常", "缺失")) & (TSH_group %in% c("正常", "缺失")) & str_detect(记录内容, "低甲状腺素血症") ~ "52",
      
      # ---- 6 组：甲状腺结节 ----
      str_detect(记录内容, "甲状腺结节") ~ "60",
      
      # ---- 7 组：甲状腺癌 ----
      str_detect(记录内容, "甲状腺癌") ~ "70",
      
      # ---- 其他 ----
      TRUE ~ "分类不明"
    )
  )

# 步骤4：按 NIPTID 聚合单人多条记录
df_subject <- df_basic_coded %>%
  group_by(NIPTID) %>%
  summarise(
    # 多条记录编码合并：相同保留一个，不同用 & 连接
    code_concat = ifelse(
      n_distinct(disease_code) == 1,
      first(disease_code),
      paste(sort(unique(disease_code)), collapse = "&")
    ),
    .groups = "drop"
  ) %>%
  mutate(
    # 清理拼接中的「分类不明」无效标签
    code_concat = map_chr(code_concat, function(x) {
      if (x == "分类不明") return(x)
      if (!str_detect(x, "&")) return(x)
      parts <- str_split(x, "&")[[1]]
      parts <- parts[parts != "分类不明"]
      if (length(parts) == 0) return("分类不明")
      paste(parts, collapse = "&")
    }),
    
    # 过滤 0 开头良性编码
    code_filtered = sapply(code_concat, filter_zero_code, USE.NAMES = FALSE),
    
    # 共病标准化重编码（表三核心）
    disease_code = sapply(code_filtered, map_comorbidity_code, USE.NAMES = FALSE),
    
    # 提取疾病大类
    disease_category = sapply(disease_code, extract_disease_category, USE.NAMES = FALSE)
  )

# 步骤5：关联样本信息 & 结果导出
df_final <- left_join(df_subject, df_sample_id, by = "NIPTID")
df_export <- df_final %>%
  select(NIPTID, disease_code, disease_category, everything())
write_xlsx(df_export, "甲状腺疾病分类_3.xlsx")

# 步骤6：质控核查
# 1. 分类不明病例核查
case_unknown <- df_export %>% filter(disease_category == "分类不明")
case_unknown_detail <- inner_join(case_unknown, df_basic_coded, by = "NIPTID")
# 2. 诊断文本同时含「甲亢+甲减」冲突病例核查
case_conflict <- df_basic_coded %>%
  filter(str_detect(记录内容, "亢") & str_detect(记录内容, "减"))
# 3. 重点亚组提取
case_hashimoto_hypo  <- df_export %>% filter(disease_code == "13")  # 桥本合并甲减
case_hashimoto_hyper <- df_export %>% filter(disease_code == "15")  # 桥本合并甲亢
case_nodule          <- df_export %>% filter(disease_category == "6")  # 结节组
# 4. 编码分布统计
cat("===== 精细编码分布 =====\n")
print(table(df_export$disease_code))
cat("\n===== 疾病大类分布 =====\n")
print(table(df_export$disease_category))
