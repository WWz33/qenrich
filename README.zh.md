[English](README.md) | 简体中文

# qenrich

面向非模式生物的快速基因富集工具（GO / KEGG / Pfam / InterPro），统计基于 [decoupler](https://github.com/scverse/decoupler) 的 ORA（Fisher 精确检验，BH 校正）。

注释文件自动识别、解析一次并缓存；基因列表一列一个基因集，一次跑完所有列。

## Install

```bash
git clone https://github.com/WWz33/qenrich.git
cd qenrich
pip install -e .            # 自动安装 decoupler 等依赖
```

## Usage

```bash
# 解析 + 缓存 + 富集 go
qenrich -i emapper.annotations.tsv --genelist gene_list.txt

# 同一文件的 KEGG
qenrich -i emapper.annotations.tsv -f kegg --genelist gene_list.txt

# 复用已解析的对象
qenrich -i go --genelist gene_list.txt --db qenrich_db/

# 解析成对象库以便复用
qenrich parse emapper.annotations.tsv -o qenrich_db/

# GO 传播 + 术语名（GO 分析建议加）
qenrich -i emapper.annotations.tsv --genelist gene_list.txt --obo go-basic.obo

# 折叠有显著子节点的父节点；给 KEGG/Pfam 补名；去 .N 后缀
qenrich -i go --genelist gene_list.txt --obo go-basic.obo --drop-parents \
        --desc ko_ids.txt --strip-suffix

# eggNOG：只保留一个注释层级
qenrich -i emapper.annotations.tsv --eggnog-lvl 2759 --genelist gene_list.txt

# 自定义背景 + 出图
qenrich -i go --genelist gene_list.txt --bg universe.txt -o out/ --plot

# 只跑指定列（列名或 1-based 序号）
qenrich -i emapper.annotations.tsv --genelist gene_list.txt -c salt_stress_up,2
```

`gene_list.txt` 空白或制表符分隔，一列一个基因集。表头自动识别，列名即集合名（无表头记为 `set1..setN`）。`gene,weight` 列（第二字段为数值，如 log2FC）走 GSEA，纯 ID 列走 ORA。`.gz` 直接读取。

纯 ID 列：

```
salt_stress_up  control_down
Glyma.01G000100 Glyma.01G000400
Glyma.01G000200 Glyma.01G000500
Glyma.01G000300 Glyma.01G000600
```

带权列（GSEA）：

```
deg_up
Glyma.01G000100,3.2
Glyma.01G000200,1.8
Glyma.01G000300,-0.5
```

每个集合写入 `<set>_enrichment.tsv`（ORA）或 `<set>_gsea.tsv`（GSEA），合并为 `summary.tsv`。加 `--desc go_zh.tsv` 后输出带英文和中文两列 name（`--no-name-zh` 去掉中文列）：

```
term          name                          name_zh          term_size  overlap  genes           pvalue    log_or  padj
GO:0006950    response to stress            对胁迫的响应      18         1        glyma…Gm…0028… 4.02e-07  -3.53   1.75e-05
GO:0043565    sequence-specific DNA binding  序列特异性DNA结合 14         1        glyma…Gm…0032… 4.52e-05  -2.97   0.000196
GO:0048519    negative regulation of bio...  生物过程的负调控   17         3        glyma…Gm…0031… 4.41e-04  -2.2    0.00174
```

`--style` 切换绘图风格：

<img src="data/png/matplotlib_en.png" width="500" alt="matplotlib 风格">
<img src="data/png/enrichplot_en.png" width="500" alt="enrichplot 风格">
<img src="data/png/heatmap_en.png" width="500" alt="汇总热图">

标签由 `--labels {name,id}` 控制（默认 `name`）。

## 中文标签

`go_zh.tsv`（38,092 行）是 go-basic.obo 全部术语名的 LLM 翻译，未经人工校对，引用前请核对。该文件在仓库根目录，用 `--desc` 指定：

```bash
qenrich -i emapper.annotations.tsv --genelist gene_list.txt \
        --desc go_zh.tsv --plot --style enrichplot
```

<img src="data/png/enrichplot_zh.png" width="500" alt="中文点图">
<img src="data/png/heatmap_zh.png" width="500" alt="中文热图">

## 支持的注释格式

| 格式 | 识别 | 产出对象 |
|---|---|---|
| eggNOG-mapper `*.emapper.annotations.tsv` | `#query` 表头 + `GOs` 列 | go, kegg, pathway, pfam, cog |
| InterProScan TSV（`--goterms`） | 无表头，第 2 列 = 32 位 md5 | go, interpro, pfam |
| Blast2GO `.annot` | 4 列，第 3 列 = `InterPro` | go, interpro |
| Blast2GO / OmicsBox tabular | `Sequence Name` 表头 | go |
| PANNZER2 `.out` | `qpid`/`qseqid` 表头 | go |
| Trinotate report | `gene_id` + `transcript_id` header | go |
| GAF 2.x（TAIR / UniProt / Ensembl Plants） | `!` 注释行，15 列，GO 第 5 列 | go |
| GFF3 `Ontology_term=` | `##gff-version`；mRNA 向上归到 gene | go, interpro |
| KofamKOALA / GhostKOALA / KAAS | 第 2 列 = `K\d{5}` | kegg |
| 通用（基因列 + ID 列） | 兜底 | 按正则分拣 go/kegg/pfam/interpro |
| 标准长表 TSV | `source`,`target` 表头 | net |

`data/format/` 每种格式示例：

```bash
qenrich -i data/format/emapper.annotations.tsv \
        --genelist data/format/gene_list.txt --tmin 3
```

## 注意

`--tmin`（默认 5）丢弃目标基因过少的 term；小注释集调低。
