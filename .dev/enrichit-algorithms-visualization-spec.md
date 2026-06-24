# enrichplot × enrichit 新算法可视化衔接 Spec

## 1. 目标

这份文档讨论的是：

**`enrichplot` 应该如何消费 `enrichit` 新增算法所产生的结果对象与 explanation-ready 数据结构。**

关注范围包括：

- 多组学早期融合后的常规富集结果
- 多组学晚期融合返回的 `enrichResult`
- 单层网络富集 `nseaResult`
- 多层网络富集 `mnseaResult`
- `mnsea` explanation helpers 所产生的中间数据

## 2. 设计前提

### 2.1 职责边界

`enrichplot` 的职责不是做算法计算，而是做：

- 可视化
- plot-facing helper
- fortify / tidy conversion
- explanation-ready 数据的图形映射

换句话说：

- **计算在 `enrichit`**
- **知识库语义入口在 `clusterProfiler`**
- **图形表达在 `enrichplot`**

### 2.2 本 spec 的关键判断

对于 `enrichplot` 来说，新算法接入的关键不是“支持多少新函数名”，而是“支持哪些新数据结构”。

现阶段最重要的对象有两类：

1. 已兼容现有绘图生态的对象
   - `enrichResult`
   - `gseaResult`
   - `nseaResult`
2. 需要专门扩展的对象
   - `mnseaResult`
   - `get_mnsea_contribution()` 的返回表
   - `extract_mnsea_subnetwork()` 的返回列表

## 3. 总体策略

### 3.1 基本原则

`enrichplot` 对新算法功能的接入分两层进行：

#### 第一层：零成本复用

凡是已经返回：

- `enrichResult`
- `gseaResult`
- `nseaResult`

且现有 `fortify()` / `dotplot()` / `cnetplot()` / `heatplot()` 语义仍然成立的对象，尽量直接复用现有图形方法。

例如：

- `aggregate_enrichment()` 返回的 `enrichResult`
- early fusion 后接 `gsea()` / `ora()` / `nsea()`

这些本来就应该尽量“不额外发明图”。

#### 第二层：为 `mnseaResult` 增加专门图形语义

`mnseaResult` 有额外结构：

- layer_scores
- collapsed_scores
- pathway_contribution
- feature_contribution
- coupling_table

这意味着它不能只被当成“另一种 `gseaResult`”。

所以 `enrichplot` 需要围绕 `mnseaResult` 设计新的 plot-facing helper 和图形方法。

## 4. 不需要新增图形支持的部分

### 4.1 晚期融合

`aggregate_enrichment()` 返回的是 `enrichResult`。

因此它天然可以复用：

- `dotplot()`
- `barplot()`
- `emapplot()`
- `cnetplot()`
- `heatplot()`

不建议为了“late fusion”单独发明新 plot。

真正需要补的只是：

- 文档说明：late fusion 结果可以直接进现有 enrichplot 工作流
- 若有必要，在 vignette 中补一小段示例

### 4.2 早期融合

早期融合下游还是普通：

- `enrichResult`
- `gseaResult`
- `nseaResult`

所以也不需要在 `enrichplot` 新增专门图形类型。

## 5. 需要专门扩展的部分：mnsea

## 5.1 为什么 `mnseaResult` 需要单独对待

`mnseaResult` 不是简单的 GSEA 结果加几个注释列。

它包含的是：

- 多层扩散后的 pathway 结果
- pathway-level contribution table
- feature-level contribution table
- explanation-ready subnetwork extraction

用户真正想问的是：

- 哪一层驱动了这个 term
- 哪些 feature 在驱动
- 层间耦合关系如何参与了信号传播

这些问题超出了现有 `dotplot(gseaResult)` 的语义。

## 5.2 设计原则

### 5.2.1 helper 放在 enrichplot，而不是 enrichit

如果 helper 只是为了可视化服务，就应该放在 `enrichplot`。

例如：

- 把 pathway contribution 转成 term × layer matrix
- 把 subnetwork list 转成 ggraph / tidygraph 可吃的 node / edge tibble
- 为分层热图、分层网络图整理 aesthetics 所需字段

而这些 helper 不应该反向塞回 `enrichit`。

### 5.2.2 不把 plotting concern 带回计算层

不建议在 `enrichit` 中加入：

- node color mapping
- shape mapping
- panel splitting
- legend-related metadata

这些都属于 `enrichplot` 的职责。

## 6. enrichplot 中拟新增的 helper

### 6.1 `fortify.mnseaResult()`

建议新增：

```r
fortify.mnseaResult <- function(
    model,
    data,
    showCategory = 10,
    by = c("p.adjust", "NES", "contribution"),
    level = c("result", "pathway"),
    layer = NULL,
    ...
)
```

#### 语义

- `level = "result"`：
  - 返回常规 term-level 结果表
  - 类似 `fortify.gseaResult()`
- `level = "pathway"`：
  - 返回 pathway contribution table
  - 适用于 layer-aware dotplot / heatplot

#### 设计目的

让 `mnseaResult` 先进入 `fortify()` 体系，再决定上层绘图函数如何消费。

### 6.2 `fortify_mnsea_contribution()`

如果不想把所有语义都压进 `fortify.mnseaResult()`，建议再加一个明确 helper：

```r
fortify_mnsea_contribution <- function(
    x,
    level = c("pathway", "feature"),
    pathway_id = NULL,
    ...
)
```

内部直接调用：

- `get_mnsea_contribution()`

再补：

- 排序
- 因子顺序
- 画图所需列名标准化

### 6.3 `fortify_mnsea_subnetwork()`

建议新增：

```r
fortify_mnsea_subnetwork <- function(
    x,
    pathway_id = NULL,
    include_couplings = TRUE,
    include_isolated = FALSE,
    ...
)
```

内部调用：

- `extract_mnsea_subnetwork()`

返回标准化列表：

- `nodes`
- `edges`

并补足可视化常用列，例如：

- node_type
- layer
- abs_score
- sign
- edge_type

## 7. enrichplot 中拟新增图形接口

这部分是 spec，不代表要一次全做完。

建议分成三个层次。

### 7.1 第一优先级：mnsea dotplot

建议新增：

```r
dotplot.mnseaResult <- function(
    object,
    x = "share",
    color = "contribution",
    showCategory = 10,
    layer = NULL,
    split = NULL,
    ...
)
```

#### 图形语义

- y 轴：term
- x 轴：layer share 或 contribution
- color：contribution / NES / p.adjust
- facet 或 split：layer

#### 为什么优先做它

因为这是最轻量、最容易落地、也最容易解释的 `mnsea` 可视化入口。

### 7.2 第二优先级：mnsea heatplot

建议新增：

```r
heatplot.mnseaResult <- function(
    x,
    pathway_id = NULL,
    showCategory = 10,
    showTop = NULL,
    value = c("score", "abs_score", "share"),
    ...
)
```

#### 图形语义

两类模式：

1. term × layer 的 contribution heatmap
2. feature × layer 的 pathway-specific heatmap

这里的 `heatplot` 更像对现有 `heatplot.enrichResult()` 的扩展，而不是完全新物种。

### 7.3 第三优先级：mnsea cnetplot / layered cnet

建议新增：

```r
cnetplot.mnseaResult <- function(
    x,
    pathway_id = NULL,
    include_couplings = TRUE,
    node_label = "all",
    ...
)
```

但这里要明确：

- 它不是现有 `cnetplot.enrichResult()` 的简单复用
- 它的输入应该来自 `fortify_mnsea_subnetwork()`

#### 视觉语义

- term 节点
- feature 节点
- layer 信息
- coupling edge / intra edge 区分

这部分复杂度最高，建议放到后期。

## 8. 与现有函数的衔接策略

### 8.1 现有 `dotplot()`

现有 `dotplot.enrichResult()` / `dotplot.gseaResult()` 本质上依赖：

- `fortify()`
- 标准列：
  - `Description`
  - `Count`
  - `GeneRatio`
  - `p.adjust`
  - `NES`

所以如果要支持 `mnseaResult`，最稳的做法是：

- 先把 `mnseaResult` 映射到 `fortify()` 体系
- 再在 `dotplot.mnseaResult()` 中定义自己的 aesthetics

### 8.2 现有 `heatplot()`

现有 `heatplot.enrichResult()` 假设输入是：

- gene set list
- optional foldChange / pvalue

而 `mnsea` 的 pathway-specific explanation 更像：

- feature × layer 矩阵

所以这里更适合做**专门方法**，不要硬复用旧逻辑。

### 8.3 现有 `cnetplot()`

现有 `cnetplot.enrichResult()` 核心输入是：

- term -> gene set mapping

而 `mnsea` 子网已经是更底层的 node/edge 数据。

因此：

- `mnsea` 的网络图语义应基于 `extract_mnsea_subnetwork()` 及其 fortify helper
- 不建议强行塞进旧的 `extract_geneSets()` 逻辑

## 9. 文档策略

### 9.1 README

README 不需要立即加一堆 `mnsea` 图。

建议只补一句定位：

- enrichplot can visualize explanation-ready outputs generated by enrichit topology-aware workflows

### 9.2 vignette

建议在 `vignettes/enrichplot.qmd` 中增加一个章节：

`Visualizing topology-aware enrichment results`

分三小节：

1. late fusion results as ordinary `enrichResult`
2. pathway-level contribution plot for `mnseaResult`
3. subnetwork-ready data for layered network visualization

### 9.3 man / examples

新方法一旦实现，示例必须用小 toy data，保持轻量。

## 10. 测试策略

### 10.1 helper tests

建议新增：

- `test-fortify-mnsea.R`
- `test-fortify-mnsea-subnetwork.R`

测试内容：

- 返回列名是否稳定
- 排序和筛选是否正确
- `pathway_id = NULL` 时默认行为是否明确

### 10.2 plotting tests

建议新增：

- `test-dotplot-mnsea.R`
- `test-heatplot-mnsea.R`

测试目标：

- 返回 `ggplot` 对象
- 数据映射列存在
- 不对图形细节做过度脆弱断言

### 10.3 不在 enrichplot 重复测试 enrichit 算法

这里的测试应只关心：

- helper 是否正确整理输入
- 图形方法是否能消费这些输入

不应在 `enrichplot` 重复验证 `mnsea` 统计计算本身。

## 11. 分阶段实施建议

### Phase 1: helper-first

- `fortify.mnseaResult()`
- `fortify_mnsea_contribution()`
- `fortify_mnsea_subnetwork()`
- vignette 草稿

### Phase 2: 最小图形支持

- `dotplot.mnseaResult()`
- 简单的 pathway contribution heatmap

### Phase 3: 复杂网络图

- `cnetplot.mnseaResult()`
- layered network visualization

## 12. 一句话总结

`enrichplot` 对新算法功能的最佳衔接方式，不是把计算搬进来，而是：

**围绕 `mnseaResult` 和 explanation-ready helpers 建立新的 fortify 与 plot-facing helper 体系，让 `enrichit` 负责算，`enrichplot` 负责把这些数据讲清楚。**
