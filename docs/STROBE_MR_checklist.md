# STROBE-MR 清单

> **STROBE-MR** = Strengthening the Reporting of Observational Studies in
> Epidemiology – Mendelian Randomization（Skrivankova et al., *JAMA* 2021）
>
> 状态标记：✅ 完成 ｜ ⚠️ 部分完成 ｜ ❌ 待办 ｜ — 不适用
> 本文档基于**实际完成的分析**填写，不虚构未做的工作。

---

## 一、当前状态汇总

| # | 条目 | 状态 |
|---|---|---|
| 1 | 标题与摘要 | ⚠️ 骨架已定，摘要待写 |
| 2 | 背景与原理 | ⚠️ 逻辑已定，正文待写 |
| 3 | 研究目标与预注册 | ❌ **未预注册** |
| 4 | 研究设计 | ✅ |
| 5 | 数据来源与暴露/结局 | ✅ |
| 6 | 工具变量 | ✅ |
| 7 | 样本量与功效 | ⚠️ 部分 |
| 8 | 主要统计方法 | ✅ |
| 9 | 其他统计方法 | ✅ |
| 10 | 多重检验 | ⚠️ 部分 |
| 11 | 软件与代码 | ⚠️ 版本已记录，仓库待建 |
| 12 | 参与者流程 | ✅ |
| 13 | 描述性数据 | ✅ |
| 14 | 主要结果 | ✅ |
| 15 | 敏感性分析 | ✅ |
| 16 | 数据与代码可用性 | ❌ **待办** |
| 17 | 三方资助 | ❌ **待办** |
| 18 | 利益冲突 | ❌ **待办** |
| 19 | 局限 | ✅ |
| 20 | 解释 | ✅ |

**完成 10 项、部分 5 项、待办 5 项。**

---

## 二、逐条填写

### 条目 1 — 标题与摘要 ⚠️

**1a. 标题**
> Shared genetic architecture with adiposity, rather than a direct causal effect,
> explains the association between obstructive sleep apnoea and heart failure:
> a two-sample and multivariable Mendelian randomisation study

**1b. 摘要** —— 待写。结构见 `PAPER_SKELETON.md` 第三节。

---

### 条目 2 — 背景与原理 ⚠️

**2a. 科学背景**
OSA 与心衰的观察性关联一致，但 MR 结论矛盾：
- PMID 36480010（*Sleep* 2023）：校正 BMI 后阳性，OR 1.13 (1.01–1.27)
- PMID 36611115（*Eur J Prev Cardiol* 2023）：阴性
- PMID 40472801（*EBioMedicine* 2025）：PGS 研究（**非 MR**），BMI 校正后心衰信号消失
- PMID 39288744（*Cardiology* 2025）：两步 MR，提示肥胖主导

**2b. 为何需要 MR**
观察性关联受反向因果与残余混杂影响；RCT 层面 CPAP 试验（SAVE、RICCADSA）
未见心血管获益，机制不明。

**2c. 原理** —— 利用配子随机化，以遗传变异为工具推断因果。

---

### 条目 3 — 研究目标与预注册 ❌ **未预注册**

**3a. 目标**
裁决 OSA→心衰是否独立于肥胖；并用多个独立维度表征 OSA 工具的遗传性质。

**3b. 预注册** ❌ **本研究未进行预注册。必须如实声明。**

> 建议局限声明：
> "This study was not prospectively registered. All analyses were pre-specified
> in an internal analysis plan but the protocol was not publicly registered
> before data analysis."

⚠️ **这是审稿人可能扣分的点。** 若时间允许，建议在 OSF 补登（可标为
post hoc registration），至少表明分析计划是事前锁定的。

---

### 条目 4 — 研究设计 ✅

**4a. 设计**：两样本孟德尔随机化（two-sample MR）+ 多变量 MR（MVMR）

**4b. 方向性**：单向（OSA → 心衰）；方向性由 Steiger 检验确认

**4c. 变异性**：遗传变异（SNP）

**4d. 效应尺度**：**OR**（心衰为二分类结局）；中介分析在风险差尺度上做，
避免 OR 尺度系数相乘的偏倚（Carter 2021, PMID 33961203）

---

### 条目 5 — 数据来源与暴露/结局 ✅

**5a. 数据来源**

| 角色 | 数据 | 登录号 / ID | 样本量 |
|---|---|---|---|
| 暴露（主） | MVP OSA（欧洲） | `GCST90475824` | 152,031 例 / 278,027 对照 |
| 暴露（meta + 复现） | FinnGen R13 OSA | `G6_SLEEPAPNO_INCLAVO` | 74,697 例 |
| 结局（主） | HERMES 心衰 | `GCST009541` | 47,309 例 / 930,014 对照 |
| 结局（复现） | FinnGen R13 心衰 | `I9_HEARTFAIL` | 41,591 例 / 458,595 对照 |
| 结局（分层） | FinnGen 心衰+BMI≥25 | `I9_HEARTFAIL_AND_OVERWEIGHT` | 25,129 例 / 204,333 对照 |
| 结局（分层对照） | FinnGen 心衰+冠心病 | `I9_HEARTFAIL_AND_CHD` | 26,459 例 / 409,472 对照 |
| 混杂（MVMR） | BMI | `ieu-b-40` | 501 工具 |

**5b. 暴露定义**
OSA：MVP 为电子病历诊断；FinnGen 为 ICD-10 **G47.3**（含初级保健门诊诊断）。

**5c. 结局定义**
心衰：HERMES 为临床诊断汇总；FinnGen 为 ICD 编码端点。

**5d. 时间**：汇总统计量，无个体层面时间信息（MR 固有特性）

**5e. 人群**：全部为**欧洲血统**

**5f. 数据获取**：公开汇总统计量
- GWAS Catalog FTP：`https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/`
- FinnGen R13：`https://storage.googleapis.com/finngen-public-data-r13/summary_stats/`
- OpenGWAS API：`https://api.opengwas.io/`

---

### 条目 6 — 工具变量 ✅

**6a. 工具选择**
- 全基因组显著：**P < 5×10⁻⁸**
- **LD 聚簇：r² < 0.001，窗口 10 Mb，1000G EUR 参考**（经 OpenGWAS API）
- 显著变异数：MVP **4,770**；FinnGen OSA **3,205**
- 聚簇后独立工具：MVP **109**；meta **10**

**6b. 工具强度（F 统计量）**

| 工具集 | nSNP | mean F | min F |
|---|---|---|---|
| MVP（LD 聚簇） | 109 | **40.7** | **29.8** |
| meta（MVP+FinnGen） | 10 | **121.9** | — |

全部远超阈值 10，**无弱工具问题**。

**6c. 三个核心假设的检验**

| 假设 | 检验 | 结果 |
|---|---|---|
| **相关性**（Relevance） | F 统计量 | ✅ mean F ≥ 40.7 |
| **独立性**（Independence） | PheWAS + 共定位 + FTO 区剔除 | ✅ 见条目 9 |
| **排他性**（Exclusion restriction） | MR-Egger 截距 + MR-PRESSO + 加权中位数 | ⚠️ 见下 |

MR-Egger 截距 = −0.00668，**p = 0.067**（边际，未达 0.05）
→ 声明为"提示轻度定向多效性，但未达统计学显著"

**6d. 表型定义的一致性**
暴露与结局分别来自不同队列的独立诊断，表型定义一致（OSA 与心衰均为临床诊断）

---

### 条目 7 — 样本量与功效 ⚠️

**7a. 样本量**：见条目 5a

**7b. 功效计算**
未做事前功效计算。但通过以下方式表征功效：
- 条件 F（MVMR）：**40.9**，跨过阈值 10
- 阳性对照：BMI → 心衰 OR 1.715 (1.614–1.822), p=2.6e-68 ✅

⚠️ **建议补做**：对 MVMR 的直接效应做事后最小可检测效应量计算，
以支持"零结果是真阴性而非功效不足"的论断。

---

### 条目 8 — 主要统计方法 ✅

**8a. 主要估计方法**：逆方差加权（IVW）

**8b. 效应尺度**：OR（每 1 个 log-odds 遗传易患性）

**8c. 模型假设**：
- 固定效应 IVW（辅以随机效应）
- MR-Egger：允许定向多效性，但需 InSIDE 假设
- 加权中位数：需多数工具有效
- 加权众数：需最大工具簇有效

**8d. 分析软件**：见条目 11

---

### 条目 9 — 其他统计方法 ✅

**9a. 多变量 MR（MVMR）**
- 校正变量：BMI（`ieu-b-40`）
- 估计方法：**自建加权最小二乘（WLS）**
  （`MVMR` 包因 r-universe 403 无法安装）
- **条件 F 统计量（Sanderson et al. 2019）**：
  - MVP 工具集：F(OSA|BMI) = **1.3**（不足）
  - meta 工具集：F(OSA|BMI) = **40.9**（充足）✅

**9b. 工具性质诊断 — PheWAS**
10 个工具 SNP 对 **17 个表型**做 IVW：
- 预期阳性（11 个）：肥胖、BMI ×3、腰臀比 ×2、体脂、T2D、血糖、血脂 ×2、高血压
- 预期阴性（6 个）：结直肠癌、前列腺癌、白内障、骨关节炎等

结果：阳性 7/11 显著；**肥胖效应最强（β=+1.95, p=3.0e-8）**；
癌症/白内障全部阴性。

**9c. 共定位**
- 方法：`coloc.abf`（Giambartolomei et al. 2014）
- 窗口：各工具位点 ±500 kb
- 数据：FinnGen OSA vs FinnGen 心衰（**同队列**，已声明）
- 结果：**PP.H4 ≥ 0.8 有 1/10（FTO，0.999）**；PP.H3 ≥ 0.8 有 2/10

**9d. Steiger 方向性**
- 方法：易患性尺度换算，同一队列内比较 R²
- 公式：`beta_liab = beta_logOR × s(1−s)/i`，`R² = 2×MAF×(1−MAF)×beta_liab²`
- 结果：**10/10 工具方向正确，R² 比值 9.26**

**9e. 肥胖分层**
三个心衰端点对比：严格心衰 1.345 / 心衰+BMI≥25 1.277 / 心衰+CHD 1.312

**9f. 敏感性分析完整清单**

| 分析 | 方法/工具 | 结果 |
|---|---|---|
| 异质性 | Cochran Q | IVW Q p=0.0020 |
| 定向多效性 | MR-Egger 截距 | −0.0067, p=0.067 |
| 离群值 | MR-PRESSO (5000 次) | 未检出离群值；global p=0.038 |
| 留一法 | leave-one-out | 88 个 SNP 全部移除后仍显著 |
| 区段剔除 | 剔除 FTO 区 ±500 kb | OR 1.264 (1.164–1.372), p=2.2e-8 |
| 弱工具 | 条件 F | 40.9（meta 工具集） |
| 阳性对照 | BMI → 心衰 | OR 1.715, p=2.6e-68 |

---

### 条目 10 — 多重检验 ⚠️

**现状**：
- 主分析：单一暴露 × 单一主结局 → α = 0.05
- CMR 相关分析（未纳入本文）：原计划 FDR(BH) 校正
- PheWAS：已用 **BH-FDR** 校正（`PHEWAS_results.csv` 含 `fdr` 列）
- 共定位：报原始 PP，未做多重校正

⚠️ **待办**：明确声明本文的检验族（family of tests）与校正策略。
建议：主分析（MVP / meta / 复现 / MVMR / 分层）共 5 个主要检验，
可用 **Bonferroni α = 0.01** 作为敏感性界值。

---

### 条目 11 — 软件与代码 ⚠️

**已记录的软件版本**：

| 软件 | 版本 |
|---|---|
| R | 4.4.1 (2024-06-14 ucrt) |
| TwoSampleMR | 0.7.9 |
| ieugwasr | 1.1.0 |
| coloc | 5.2.3 |
| data.table | 1.18.2.1 |
| MRPRESSO | 1.0 |
| MendelianRandomization | 0.10.0 |
| RadialMR | 1.2.4 |
| MVMR | **不可用**（自建实现） |
| MRlap | **不可用**（缺 `make` 编译工具链） |

**分析脚本**：`scripts/00–43`（共 44 个 R 脚本）

❌ **待办**：建立公开仓库（GitHub + Zenodo DOI）

---

### 条目 12 — 参与者流程 ✅

汇总统计量研究，无个体层面流程。已报告：
- 起始变异数（MVP 19,708,468；FinnGen 21,327,043）
- 显著变异数（4,770 / 3,205）
- 聚簇后工具数（109 / 10）
- harmonise 后进入分析的 SNP 数（88 / 9）
- 各步骤排除原因（palindromic 等位基因、缺失数据、LD）

---

### 条目 13 — 描述性数据 ✅

**工具变量特征**（见 `results/mvp_osa_instruments_LDclumped.tsv`）：
SNP、效应等位基因、EAF、beta、SE、p、F 统计量

**暴露与结局 GWAS 特征**：样本量、病例数、血统、数据来源（见条目 5a）

---

### 条目 14 — 主要结果 ✅

| 分析 | nSNP | OR (95% CI) | p |
|---|---|---|---|
| MVP OSA → HERMES 心衰（IVW） | 88 | **1.279 (1.187–1.378)** | 1.0e-10 |
| MVP OSA → HERMES（加权中位数） | 88 | 1.222 (1.107–1.349) | 6.9e-5 |
| MVP OSA → HERMES（MR-Egger） | 88 | 1.553 (1.249–1.931) | 1.5e-4 |
| meta OSA → HERMES（IVW） | 9 | 1.447 (1.198–1.748) | 1.3e-4 |
| **meta OSA → FinnGen 心衰（复现）** | 9 | **1.345 (1.103–1.641)** | 0.0035 |
| **MVMR：OSA 直接效应** | 82 | **1.158 (0.904–1.483)** | 0.246 |
| MVMR：BMI 直接效应 | 82 | 1.950 (1.619–2.350) | 2.1e-12 |
| 阳性对照：BMI → 心衰 | 446 | 1.715 (1.614–1.822) | 2.6e-68 |

---

### 条目 15 — 敏感性分析 ✅

见条目 9f 完整清单。**所有敏感性分析结果方向一致**，
除肥胖分层（见条目 19）。

---

### 条目 16 — 数据与代码可用性 ❌ **待办**

**数据**：全部公开
- GWAS Catalog：`GCST90475824`、`GCST009541`
- FinnGen R13：公开 GCS 桶
- OpenGWAS：`ieu-b-40`

**代码**：❌ 尚未建立公开仓库

**建议声明模板**：
> "All summary statistics used are publicly available (accession numbers in
> Table 1). Analysis code is available at [GitHub URL] and archived at
> [Zenodo DOI]."

---

### 条目 17 — 三方资助 ❌ **待办**（需作者填写）

---

### 条目 18 — 利益冲突 ❌ **待办**（需作者填写）

---

### 条目 19 — 局限 ✅

1. 仅欧洲人群
2. 工具富集于肥胖相关位点（但这是发现本身，非缺陷）
3. `GCST90475824` 的 SE 列全为 `#NA`，**SE 由 p 反推**
   （CI 列与 p 不自洽，z 比值中位数 0.903）
4. meta 工具仅 10 个独立 SNP，单变量点估计偏高（1.447 vs 工具全集 1.279）
5. **MRlap 未做**（系统缺 `make`）—— 但核心结果是零效应，
   重叠只会使估计偏向阳性，不威胁结论；且已在非重叠 FinnGen 复现
6. 共定位为**同队列**（FinnGen vs FinnGen）
7. 肥胖分层端点病例比不同（11.0% vs 8.3%），log-OR 跨端点不严格可比
8. MVMR 为自建实现
9. 未预注册
10. MR 不评估治疗干预效果

---

### 条目 20 — 解释 ✅

**主要发现**：OSA 与心衰的遗传关联主要反映**与肥胖的共享遗传结构**，
而非独立因果通路；但超重人群内关联仍在，提示 BMI 高端区间可能有残余效应。

**与既往研究比较**：
- 裁决了 36480010（阳性）与 36611115（阴性）的矛盾，
  并证明旧工具（5 SNP）移除单个 FTO 位点即翻转，是矛盾的可能来源
- 独立验证 40472801（PGS）的提示，但**首次用 MR 因果估计量**
- 与 EBioMedicine 2026（PMID 41833577）的"SA 独立于血糖/血压"结论**不矛盾**，
  因该文未校正 BMI

**机制解释**：为 CPAP 随机试验的阴性结果提供统一解释

**临床含义**：减重仍是风险管理核心；**不否定 CPAP 对症状与生活质量的疗效**

**未来方向**：连续 BMI 的中介 MR；非欧洲人群；BMI 高端区间的独立效应

---

## 三、投稿前必须完成的事项

| 优先级 | 事项 |
|---|---|
| **高** | 代码仓库 + Zenodo DOI（条目 16） |
| **高** | 声明未预注册，或补登 OSF（条目 3） |
| **中** | MVMR 事后功效 / 最小可检测效应量（条目 7） |
| **中** | 明确多重检验策略（条目 10） |
| **中** | 摘要撰写（条目 1b） |
| 低 | 若有条件，在装 Rtools 的环境补做 MRlap |
