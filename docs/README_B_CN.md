# Week9 B：FP 与 TP 的标准 metagene 分布（每类 >50）

## 1. 输入

- 四分类矩阵：`/home/xiyao/15data/W8/results/01_four_class_global_fdr/four_class_matrix_Total_gt_10.rds`
- no-cali、IVT-matched benchmark：`/home/xichen/ConquerNet/ConquerBN_basic/m6AConquer_hg38_BenchmarkDataSet_Combined_SE_no_cali_IVT_matched.rds`
- OmicsFeatures：`/home/xiyao/15data/W8/input_data/features/m6Aconquer_omicsFeatures_hg38.rds`
- 三者均为2,788,705个位点，坐标和行顺序完全一致。

FP取class code 2，TP取class code 3；每项技术再次核对`Total>10`，不使用fallback。

## 2. 标准坐标

对具有互斥5′UTR、CDS或3′UTR注释的位点，将`meta_tx_topology`映射为：

```text
5′UTR [0,1/3]
CDS    [1/3,2/3]
3′UTR  [2/3,1]
```

编码审计显示三类注释位点100%落在各自预期区间。TP/FP分别计算密度，因此曲线表示类别内部位置分布。

## 3. 纳入标准

老师确认不必要求每类1,000，只需两类都大于50。因此主图包含8项技术。代码使用`min_class=51`，与整数计数`>50`严格等价。

| Technology | FP | TP | 说明 |
|---|---:|---:|---|
| eTAM | 30,424 | 220,701 | 大样本 |
| m6A-REF | 29,628 | 20,827 | 大样本 |
| m6A-SAC | 62 | 81,098 | 小FP，谨慎解释 |
| m6ACE | 15,400 | 217,236 | 大样本 |
| MAZTER | 7,071 | 2,880 | 可比较 |
| MeRIP | 15,229 | 255,812 | 大样本 |
| scDART | 3,745 | 4,984 | 可比较 |
| DeepRM | 236 | 72,881 | 小FP，谨慎解释 |

主图标注的n是通过对应转录本区域annotation mask后的assigned n，因此会略低于这张标签总数表。

## 4. 效应量

方向统一为`TP−FP`；KS D只表示分布差异大小，Cliff's delta为正表示TP位置整体更高。

| Technology | Meta-tx mean diff | KS D | Cliff δ | Full-gene mean diff | KS D | Cliff δ |
|---|---:|---:|---:|---:|---:|---:|
| eTAM | +0.008 | 0.086 | +0.021 | +0.018 | 0.056 | +0.052 |
| m6A-REF | −0.010 | 0.087 | −0.040 | +0.018 | 0.067 | +0.052 |
| m6A-SAC | −0.099 | 0.298 | −0.330 | −0.087 | 0.204 | −0.225 |
| m6ACE | −0.022 | 0.056 | −0.029 | −0.052 | 0.067 | −0.079 |
| MAZTER | +0.028 | 0.139 | +0.092 | +0.082 | 0.170 | +0.196 |
| MeRIP | −0.014 | 0.047 | −0.012 | −0.051 | 0.065 | −0.075 |
| scDART | −0.092 | 0.287 | −0.373 | −0.073 | 0.120 | −0.177 |
| DeepRM | +0.073 | 0.269 | +0.234 | +0.070 | 0.272 | +0.242 |

## 5. 解释

- 技术之间不存在统一方向。
- scDART、m6A-SAC和DeepRM显示较大的meta-transcript差异。
- 但m6A-SAC和DeepRM的FP很少，密度曲线和效应量对少数位点敏感，只能作为初步结果。
- eTAM、m6A-REF、m6ACE和MeRIP的差异较小；MAZTER为中等差异。
- assignment rate在FP/TP之间可能不同，因此图和效应量必须与`metagene_sample_audit.csv`一起审查。

## 6. 结果路径

```text
/home/xiyao/15data/W8/results/02_fp_tp_meta_distributions/
Total_gt_10_no_cali_IVT_matched/standard_metagene_min51/

/home/xiyao/15data/W8/results/02_fp_tp_meta_distributions/
Total_gt_10_no_cali_IVT_matched/effect_sizes_min51/
```

不要使用`Total_gt_10_INVALID_wrong_dataset`。
