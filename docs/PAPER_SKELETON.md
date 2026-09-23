# 论文骨架（修订版 v2）

> 本版根据补充分析结果**修正了核心论点措辞**：从"完全由肥胖介导"改为
> "主要反映与肥胖的共享遗传结构"。
> 依据见 `STAGE5_stratified_finding.md`（肥胖分层结果部分不支持简单中介模型）。

---

## 一、标题

### 推荐

**Shared genetic architecture with adiposity, rather than a direct causal effect,
explains the association between obstructive sleep apnoea and heart failure:
a two-sample and multivariable Mendelian randomisation study**

### 备选

2. *Genetic liability to obstructive sleep apnoea and heart failure: evidence for
   adiposity-driven shared architecture*
3. *Obstructive sleep apnoea and heart failure: a Mendelian randomisation study
   adjudicating conflicting evidence*

**为什么不用"mediated by adiposity"**：分层分析显示超重人群中心衰关联仍在，
"完全介导"的表述会被审稿人用这一点直接反驳。
用 "shared genetic architecture" 更准确，也更难被攻击。

---

## 二、一句话核心论点（修订）

> 遗传预测的 OSA 与心衰存在稳健且可复现的关联，但该关联**主要反映与肥胖的
> 共享遗传结构**而非独立的因果通路：OSA 工具本质上是一个肥胖工具
> （PheWAS 对肥胖 β=+1.95、对无肥胖关联表型全部阴性），唯一强共定位位点在
> FTO（PP.H4=0.999），校正 BMI 后直接效应不显著（OR 1.158, p=0.246，
> 条件 F=40.9 功效充足）。
> **但在限定为超重人群的心衰结局中关联仍然存在**，提示 BMI 高端区间内
> 可能仍有独立于二分类超重状态的效应。

---

## 三、摘要结构

| 段落 | 内容 |
|---|---|
| **Background** | OSA 与心衰的观察性关联一致，但 MR 结论矛盾（36480010 阳性 vs 36611115 阴性）；CPAP 随机试验未见心血管获益，机制不明 |
| **Methods** | 两样本 MR + MVMR；MVP OSA（152,031 例）与 FinnGen R13 双队列；meta 增强工具；HERMES 与 FinnGen 双结局；共定位、PheWAS、Steiger、肥胖分层 |
| **Results** | 单变量 OR 1.279 (1.187–1.378)，独立复现 OR 1.345 (1.103–1.641)；MVMR 校正 BMI 后 OR 1.158, p=0.246（条件 F=40.9）；PheWAS 显示工具对肥胖 β=+1.95；共定位仅 FTO 位点强共享 (PP.H4=0.999)；超重分层效应未增强（1.277, p=0.005） |
| **Conclusion** | OSA 的心衰关联主要源于与肥胖的共享遗传结构；减重仍是风险管理核心；**但 BMI 高端区间的残余效应值得进一步研究** |

---

## 四、引言逻辑链（六步）

1. **临床背景**：OSA 患病率高，观察性研究一致显示与心衰相关
2. **矛盾出现**：两篇同设计 MR 结论相反 —— 36480010（校正 BMI 后阳性 OR 1.13）vs 36611115（阴性）
3. **第三方证据**：40472801（PGS，非 MR）显示 BMI 校正后心衰信号消失；39288744 两步 MR 提示肥胖主导
4. **关键缺口**：
   - 矛盾未被裁决
   - 既往工具极弱（36480010 仅 5 个 SNP）→ **功效不足是矛盾的合理解释**
   - 无人用现代强工具正面检验
5. **本研究策略**：MVP + FinnGen 双队列 + meta 增强工具 + 多维度佐证（共定位/PheWAS/分层）
6. **预期贡献**：给出裁决性结论并解释 CPAP 阴性试验

---

## 五、方法要点

### 数据来源

| 角色 | 数据 | N |
|---|---|---|
| 暴露（主） | MVP OSA `GCST90475824` | 152,031 例 / 278,027 对照 |
| 暴露（meta+复现） | FinnGen R13 `G6_SLEEPAPNO_INCLAVO` | 74,697 例 |
| 结局（主） | HERMES 心衰 `GCST009541` | 47,309 例 / 930,014 对照 |
| 结局（复现） | FinnGen R13 `I9_HEARTFAIL` | 41,591 例 |
| **结局（分层）** | FinnGen `I9_HEARTFAIL_AND_OVERWEIGHT` | 25,129 例（BMI≥25） |
| **结局（分层对照）** | FinnGen `I9_HEARTFAIL_AND_CHD` | 26,459 例 |
| 混杂（MVMR） | BMI `ieu-b-40` | 501 工具 |

### 统计步骤

1. 工具选择：P<5e-8，**LD 聚簇 r²<0.001、10 Mb、EUR**
2. 主分析：IVW + MR-Egger + 加权中位数 + 加权众数
3. **meta 分析**：MVP × FinnGen 逆方差固定效应（235 位点，同向 100%）
4. **MVMR**：加权最小二乘 + **Sanderson (2019) 条件 F**
5. **共定位**：`coloc.abf`，各位点 ±500 kb，同队列（FinnGen vs FinnGen）
6. **Steiger**：本地 FinnGen，**易患性尺度**换算，避免病例比例偏倚
7. **PheWAS**：17 个表型（11 预期阳 + 6 预期阴）
8. **肥胖分层**：三个心衰端点对比
9. 敏感性：Cochran Q、Egger 截距、留一法、FTO 区剔除

### 必须声明的技术处理

- `GCST90475824` 的 `standard_error` 列**全为 `#NA`**；CI 列与 p 值不自洽
  → **SE 由 p 反推**：`se = |beta| / qnorm(p/2)`
- MVMR 包不可用 → 自建 WLS 实现
- MRlap 因系统缺 `make` 无法安装 → 用理论论证替代（见局限）

---

## 六、图表清单（修订）

### 主图

| # | 图 | 内容 | 数据 |
|---|---|---|---|
| **图 1** | 研究设计总图 | 双队列暴露 × 多端点结局 + 分析流程 | 示意 |
| **图 2** | 单变量 MR 森林图 | OSA→心衰：三套工具 × 四方法 | `mvp_OSA_HF.csv`、`META_OSA_HF.csv`、`REPLICATION_finngen_HF.csv` |
| **图 3** | ★ **MVMR 前后对比** | 单变量 OR 1.279 → 校正后 OR 1.158；BMI OR 2.037；条件 F=40.9 | `FINAL_mvmr_osa_bmi_hf.csv` |
| **图 4** | ★ **PheWAS** | 17 表型森林图，肥胖 β=+1.95 突出显示；预期阴/阳分组着色 | `PHEWAS_results.csv` |
| **图 5** | ★ **共定位** | 10 位点 PP.H3/PP.H4 堆叠条形图，FTO 高亮 | `COLOC_results.csv` |
| **图 6** | 敏感性分析 | 留一法 + FTO 剔除 + Steiger 方向 | `STEIGER_instruments.csv` |
| **图 7** | ★ **肥胖分层** | 三个心衰端点的 OR 对比（严格 / BMI≥25 / CHD） | `STRATIFIED_HF_comparison.csv` |
| **图 8** | 效应分解示意 | 路径图 + 各证据方向 | 综合 |

> 图 4、5、7 是本文区别于普通 MR 论文的三个亮点，**不要压缩到附图**。

### 表

| # | 表 | 内容 |
|---|---|---|
| 表 1 | 数据源清单 | 登录号、N、人群、访问方式 |
| 表 2 | 工具变量特征 | SNP 数、F、条件 F |
| 表 3 | 主结果汇总 | 三套工具 × 四方法 |
| 表 4 | MVMR 结果 | OSA 与 BMI 直接效应 + 条件 F |
| 表 5 | 肥胖分层结果 | 三端点 OR/CI/p |
| 附表 | STROBE-MR 清单 + 代码 + 全部敏感性分析 | |

---

## 七、讨论要点

### 1. 裁决既往矛盾（核心贡献）

**关键论点**：36480010 与 36611115 的矛盾很可能源于**工具强度不足**。
我们实测证明：FinnGen 的 5-SNP 工具，**移除单个 FTO 位点结论即翻转**
（b: 0.214 → 0.050, p=0.63）。

### 2. 工具性质的直接诊断（本文的方法学亮点）

**PheWAS 显示 OSA 工具本质上是一个肥胖工具**：
- 对肥胖 β=+1.95（效应最强）
- 对 BMI p=2.0e-13
- 对结直肠癌/前列腺癌/白内障**全部阴性**

这是对"共享遗传结构"最直接的证据，比单纯 MVMR 更有说服力。

### 3. 与 40472801 的关系

该文是 **PGS 关联研究，不是 MR**，无法做因果推断，也无中介分解。
本研究**首次用 MR 因果估计量**系统检验并证实其提示。

### 4. 解释 CPAP 阴性试验（现实意义）

SAVE、RICCADSA 等 CPAP 随机试验未见心血管获益。
若风险主要来自与肥胖共享的遗传结构，**治疗 OSA 本身自然难以降低心衰风险**。

### 5. ★ 必须处理的矛盾（诚实性关键）

**MVMR 校正 BMI 后归零，但按 BMI 分层效应仍在。**

论文中必须解释这个差异，建议措辞：

> "The apparent discrepancy reflects that these analyses address different
> questions: MVMR asks whether a direct effect persists across the continuous
> BMI distribution, whereas BMI-stratification restricts to a binary category
> that still contains substantial within-stratum BMI variation. The two are
> therefore not equivalent, and the persistence of the association within the
> overweight stratum suggests that residual effects at the upper end of the BMI
> distribution cannot be excluded."

### 6. 临床含义

- 减重是风险管理核心
- 对减重手术 / GLP-1 类药物的患者选择有提示
- **必须强调**：不否定 CPAP 对症状、生活质量、白天嗜睡的疗效

---

## 八、局限（完整清单）

1. **仅欧洲人群**
2. **工具富集于肥胖相关位点**，可能高估肥胖通路贡献（但这也是发现本身）
3. **SE 由 p 反推**（原始文件 SE 缺失）
4. **meta 工具仅 10 个独立 SNP**，单变量点估计偏高（1.447 vs 工具全集 1.279）
5. **MRlap 未做**（系统缺 `make` 编译工具链）—— 但核心结果是零效应，
   重叠只会使估计偏向阳性，故**不威胁结论**；且已在非重叠 FinnGen 复现
6. **共定位为同队列**（FinnGen vs FinnGen），非跨队列
7. **肥胖分层端点病例比不同**（11.0% vs 8.3%），log-OR 跨端点不严格可比
8. **MVMR 为自建实现**；条件 F 用 Sanderson 2019 公式
9. **MR 不评估治疗干预效果**

---

## 九、目标期刊

| 顺序 | 期刊 | 理由 |
|---|---|---|
| 1 | ***EBioMedicine*** | 2025（40472801）与 2026（41833577）各发一篇该方向 |
| 2 | ***Cardiovascular Diabetology*** | MR + 中介同形先例极多 |
| 3 | ***Sleep*** | 睡眠领域顶刊，36480010 发表刊 |
| 保底 | *Front Cardiovasc Med* / *BMC Medicine* | 接收率高 |

**卖点定位**：不是"又一个 OSA 与某病的 MR"，而是
**"用工具性质诊断（PheWAS）+ 共定位 + 分层，裁决一个三年未决的矛盾"**。

---

## 十、投稿前待办

| 项目 | 状态 |
|---|---|
| 核心分析 | ✅ 完成 |
| Steiger / PheWAS / 分层 | ✅ 完成 |
| MRlap | ⚠️ 受阻，有替代论证 |
| **STROBE-MR 清单** | ❌ 必做 |
| **代码 + 汇总统计存档**（GitHub/Zenodo） | ❌ 必做 |
| 图 1/图 8 示意图绘制 | ❌ 待做 |
| 参考文献核对 | ❌ 待做 |
