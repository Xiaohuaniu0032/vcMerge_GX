# vcMerge_GX

### vcMerge plugin on Genexus

### 插件目的
如果一张芯片上了多个新冠样本，为了检查哪些突变是共检出，哪些突变是只有少数样本检出，需要肉眼比较不同样本的VCF文件，比较耗时。

本插件的目的是合并一个报告内的多个样本的新冠变异检测结果，方便比较变异检测结果。

### 插件过程
1) get `results_dir` and `analysis_dir` info from `startplugin.json`
2) get all `Snvindel.tsv` from `results_dir`
3) check if `Snvindel.tsv` is SARS sample. if not, skip this sample
4) record all SARS samples' variants info
5) output the final file `*.TSVC_variants.merged.vcf.xls`

### 关于Snvindel.tsv文件
1) skip `/ABSENT/` line
2) skip top2 line (blank line and header line)
3) skip line not belong to `2019-nCoV`

### 输出文件格式
```
Chrom	Position	Ref	Variant	training-test-flf-1_LibXXX	training-test-flf-2_LibXXX
2019-nCoV	29416572	T	C	0.999	0.999
2019-nCoV	140476936	G	A	0.996	NA

```
> 对于一个变异位点，如果该样本未检测到，则频率为NA。



### 结果下载

> 点击蓝色链接下载即可。xls格式。

![variantCallerMerge](https://github.com/Xiaohuaniu0032/vcMerge/blob/main/variantCallerMerge.png)


### 插件安装
> 假设当前目录为`/data/fulongfei/git_repo/vcMerge_GX`，你已经做完了全部的修改，下面可以打包插件，步骤如下：
1. `cd ..`
2. `zip -r vcMerge_GX.zip vcMerge_GX -x \*.git\*`
3. `mkdir vcMerge_GX_v1.1`
4. `cp vcMerge_GX.zip ./vcMerge_GX_v1.1`
5. add `asset-metadata.json` file
6. `zip -r vcMerge_GX_v1.1.zip vcMerge_GX_v1.1`
7. cp `vcMerge_GX_v1.1.zip` into GX's `/media/usbinstall/`
8. on TS, Configure -> Software Updates -> Software Updates -> USB (click APP STORE) -> install
