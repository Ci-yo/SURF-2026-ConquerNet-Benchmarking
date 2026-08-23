# Week9 A+B 组会汇报稿（每类 >50）

## A：直接四分类 posterior 与特征重要性

本周使用更新后的方法，直接计算每项技术 TN、FN、FP、TP 四个联合 posterior。对每个类别，以 `1−posterior` 为 local error，排序后计算累计平均错误率，并以 global FDR 小于 0.05 定义高置信类别。只保留唯一类别，同时要求对应技术 `Total>10`；不使用 MAP/Viterbi fallback。

老师进一步明确，机器学习不需要等到两类都达到1,000，只要当前比较的两类都大于50即可运行。因此我们保留原有22个模型，只补跑 m6A-SAC 和 DeepRM 的 TP-vs-FP、FP-vs-TN，共新增4个模型。最终完成 TP-vs-FP 8项、FP-vs-TN 7项、TP-vs-TN 11项，共26个模型；全部通过 chromosome-blocked 5-fold OOF 完整性检查。

三个任务的 macro AUROC 分别为0.894、0.799和0.932。m6A-SAC只有62个FP，DeepRM只有236个FP，因此二者虽达到运行门槛，仍必须结合类别数解释，尤其不能把TP-vs-FP接近1的AUPRC当成接近完美的性能。

稳定重要特征仍以转录本结构为主。TP-vs-FP中 `meta_tx_topology` 在8/8技术稳定；FP-vs-TN中它在7/7技术稳定。exon长度、5′/3′距离、3′UTR相对位置和CDS GC含量也反复出现。

## B：8项技术的 FP/TP 标准 metagene

使用A的同一套 no-cali posterior-derived 标签，把FP与TP位点标准化映射到5′UTR、CDS和3′UTR。按两类都大于50的标准，主图纳入8项技术：eTAM、m6A-REF、m6A-SAC、m6ACE、MAZTER、MeRIP、scDART和DeepRM。

结果仍显示明显的technology-specific模式，没有所有技术一致的统一方向。scDART的meta-transcript差异最大；m6A-SAC与DeepRM也显示较明显曲线差异，但二者FP分别只有62和236个，因此只能谨慎解释，不能与大样本技术等量比较。

需要强调：这里的FP/TP是ConquerBN内部posterior-derived标签。当前结果描述的是模型标签与转录本位置及OmicsFeatures的关联，不是实验确认的真实假阳性机制。
