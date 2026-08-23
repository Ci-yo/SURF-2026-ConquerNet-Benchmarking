# Week9 A：直接四分类 posterior 与特征重要性分析

## 1. 任务目标

使用更新后的四类技术评估方法，直接计算每项技术的 TN、FN、FP、TP 联合 posterior；对四类分别应用 global FDR `<0.05` 得到高置信位点，再用这些标签重复随机森林特征重要性分析。

## 2. 为什么使用 no-cali 数据

IVT 将用于独立验证模型产生的 TP/FP 标签。若 full ConquerBN 训练阶段已经使用 IVT calibration，再用 IVT 验证会形成信息循环。因此正式分析使用未经过 IVT calibration、但与 IVT 位点坐标匹配的 no-cali benchmark 和对应 full model。

## 3. 数据

```text
候选位点：2,788,705
技术：14
初始 OmicsFeatures：224
RF screening features：185
```

185 个特征由 224 个初始特征删除 15 个 RBP 和 24 个 histone modification 特征得到。Chromosome 只用于 blocked folds；原始 chromosome/coordinate/strand 不作为 predictor。

## 4. 四分类 posterior

设 `R=TrueM6A root state`、`S=strategy state`、`T=technology state`、`e=全部证据`。利用 `R→S→T`：

```text
P(R=r,T=t | e)
= Σs P(R=r,S=s | e) × P(T=t,S=s | e) / P(S=s | e)
```

| R | T | Class |
|---:|---:|---|
| 0 | 0 | TN |
| 1 | 0 | FN |
| 0 | 1 | FP |
| 1 | 1 | TP |

R/T 的 0 或 1 是联合潜在状态的取值，不是把边际 posterior 以 0.5 硬切分。每个位点首先保留四个完整概率。

数值审计：

```text
max |TN+FN+FP+TP−1|         = 7.2×10^-15
max |TP+FN−P(R=1)|          = 1.4×10^-14
max |TP+FP−P(T=1)|          = 2.2×10^-16
```

## 5. Global-FDR calling 与 coverage

对每项技术和类别 `c`：

```text
local error_i(c) = 1 − P_i(c)
global FDR_k(c)  = 前 k 个最小 local error 的累计平均
```

保留 global FDR `<0.05` 的 class call。只有恰好一个类别入选的位点获得最终标签；本次所有技术 overlap 均为 0。随后对对应技术独立要求 `Total>10`。

```text
0=TN, 1=FN, 2=FP, 3=TP, NA=unclassified
```

本次不使用 MAP/Viterbi fallback。老师最新确认的 `>50` 是是否运行 ML 的门槛，不是 global-FDR calling 的组成部分，也不改变任何位点标签。

## 6. 三个 RF 任务

| Task | Positive | Negative | 目的 |
|---|---|---|---|
| TP vs FP | TP | FP | 技术阳性是否得到 root 状态支持 |
| FP vs TN | FP | TN | 技术误报机制 |
| TP vs TN | TP | TN | 总体阳性/阴性结构差异 |

运行标准：

```text
纳入模型：当前任务的两类都 >50（代码 min_reliable=51）
低样本：两类均非零但至少一类 ≤50，不纳入正式运行
不可训练：任一类为0
```

这里“类似 leave-one-out”通过 chromosome-blocked 5-fold 实现：每轮留下一个未见 chromosome block 作验证，不是逐位点 leave-one-out。

## 7. RF 设计与 QC

```text
chromosome-blocked 5-fold
4 folds training + 1 fold validation
training 每类最多100,000
validation保留自然prevalence
300 trees
permutation importance
每个位点恰好获得一次OOF prediction
```

26 个完成模型均有 5 folds；OOF row index 无重复，预测无缺失，概率均在 `[0,1]`。26 个模型均导出 Top20，共 520 行，完整性审计全部通过。

## 8. 结果

| Task | 完成 | 低样本未跑 | 单类别 | Macro AUROC | Macro AUPRC |
|---|---:|---:|---:|---:|---:|
| TP vs FP | 8 | 1 | 5 | 0.8938 | 0.9529 |
| FP vs TN | 7 | 2 | 5 | 0.7988 | 0.3163 |
| TP vs TN | 11 | 1 | 2 | 0.9319 | 0.7676 |

新加入的 m6A-SAC 与 DeepRM 结果：

```text
m6A-SAC TP-vs-FP  AUROC=0.7853；FP-vs-TN AUROC=0.5663
DeepRM  TP-vs-FP  AUROC=0.9114；FP-vs-TN AUROC=0.7539
```

m6A-SAC 只有 62 个 FP，DeepRM 只有 236 个 FP。它们符合老师的运行标准，但估计方差和 fold-to-fold 稳定性不可能与上万 FP 的技术相同。TP-vs-FP AUPRC 接近 1 主要受 TP prevalence 接近 1 影响；FP-vs-TN AUPRC 很低则受 FP prevalence 极低影响。

## 9. 稳定重要特征

候选稳定规则为至少 4/5 folds importance `>0` 且至少 4/5 folds 进入 Top20。

TP-vs-FP（8项技术）：

```text
8/8 meta_tx_topology
7/8 log2_length_exons
7/8 log2_dist3prime_exons
6/8 log2_dist5prime_exons
```

FP-vs-TN（7项技术）：

```text
7/7 meta_tx_topology
6/7 relativePOS_exonicThreePrimeUTR
5/7 log2_dist5prime_exonicThreePrimeUTR
5/7 gc_content_exonicCDS
5/7 relativePOS_exonicCDS
```

这些仍是 185-feature screening 中的稳定候选，不是相关性去冗余后的最终 reduced feature set。Permutation importance 给重要程度，不给作用方向；方向需要 SHAP 或单变量效应补充。

## 10. 结论

直接四分类 posterior 标签下，OmicsFeatures 能够区分部分技术的 TP/FP 和 FP/TN，且转录本结构、meta-transcript 位置和局部序列组成反复出现。新增小 FP 技术的结果说明模型可以运行，但也展示了小类别下性能与重要性估计的不稳定性。所有类别仍是 ConquerBN posterior-derived 内部标签，不是独立湿实验真值。
