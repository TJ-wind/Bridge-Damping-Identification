# Bridge Damping Identification from Exciter Tests

本项目提供用于激振试验的主梁阻尼比识别程序。激振器将目标模态激起后停止驱动，振子仍与主梁连接。程序用主梁和振子在这一阶段的同步自由衰减记录，识别桥梁–激振器系统的耦合模态，并重构主梁自身的阻尼比。

This MATLAB code estimates bridge girder damping from forced-vibration tests in which the exciter remains attached during free decay. It fits the complex Fourier spectra of the bridge and exciter responses, then reconstructs the uncoupled bridge parameters while accounting for non-classical damping.

**使用本程序或发表由本程序获得的结果时，必须引用下方论文。 If you use this code, you must cite the accompanying paper.**

## 程序与运行环境

| 文件 | 用途 |
| --- | --- |
| [`main.m`](main.m) | 导入双通道数据、设置激振器参数，运行辨识并输出诊断图。 |
| [`main_function.m`](main_function.m) | 交互选择自由衰减窗口，拟合双通道复数频谱，调用物理参数重构并保存结果。 |
| [`refit_record5_unequal_free_k2.m`](refit_record5_unequal_free_k2.m)、[`+physical_refit_support/`](+physical_refit_support/) | 重构主梁和激振器的物理参数。请保持文件和包目录的相对位置。 |

使用支持图形窗口的 MATLAB。完整识别主梁阻尼比需要 Optimization Toolbox（`fmincon`）。Global Optimization Toolbox 可选。安装并授权后，程序使用 `particleswarm` 进行搜索，否则使用随代码附带的差分进化实现。

## 使用方法

1. 准备主梁、激振器振子两路同步响应。两路均以固定参考系为准，数据为有限实数向量，使用相同的均匀采样频率。首个样点对应同一实际时刻，记录长度可以不同。
2. 在 MATLAB 中将当前文件夹切换到仓库根目录。打开 `main.m`，将数据输入区中的示例 `load(...)` 路径改为自己的文件，并设置 `structureSignal`（主梁）、`exciterSignal`（振子）和 `fs`（Hz）。示例读取 MAT 文件中的 `beam` 和 `damper`，实测记录需自行提供。
3. 按实测工况更新 `physicalOptions`。`m2` 为已知振子质量（kg），`hasInerter` 表示是否接入惯容器。若接入惯容器并使用已知调谐条件（`useKnownAttachment = true`），填写裸振子固定底座频率和调谐后的频率（均为 Hz）。`responseOrders` 按主梁、振子顺序填写，`0` 为位移、`1` 为速度、`2` 为加速度。`inputScaleToSI` 同样按通道顺序填写，倍率包含单位换算与测量方向的符号。脚本里的采样率和参数属于原示例工况。
4. 运行 `main.m`。在图中依次点击停振后两路自由衰减的共同起点、主梁通道终点、振子通道终点。共同起点用于保留两路响应之间的相位关系。每路选区至少保留 32 个采样点，并且响应须有变化。选择终点时，主梁通道应保留可辨认的慢衰减，振子通道应覆盖早期可辨认的相对运动。

已在工作区导入两路信号并设置好试验参数时，也可以直接调用 `result = main_function(structureSignal, exciterSignal, fs, physicalOptions)`。

## 查看主梁阻尼比

`result.physicalFit.model.uncoupledBridgeDamping` 为无量纲主梁阻尼比。乘以 100 即为百分数。

```matlab
bridgeDampingPercent = 100 * result.physicalFit.model.uncoupledBridgeDamping;
```

| 结果文件 | 内容 |
| --- | --- |
| `physical_search_parameters.csv` | `zetaB` 为主梁自身阻尼比，`fB_Hz` 为主梁自身固有频率。 |
| `physical_parameters.csv` | 重构得到的质量、刚度、阻尼系数和惯容参数。 |
| `physical_fit_quality.csv` | 两通道时域和所选频带的拟合误差。 |
| `modal_parameters.csv` | 桥梁–激振器系统两阶耦合模态的固有频率和阻尼比。 |

`main_function.m` 将 MAT 结果与 CSV 表格保存到当前目录下的 `main_function_results/`。`main.m` 将诊断图和结果另存到 `main_plots/`。解读结果时，结合两路时间历程和频谱的拟合误差，以及搜索边界诊断。若振子与主梁几乎同动，振子相对运动偏弱，同一次试验可能无法可靠识别主梁模态质量。

## 引用 · Citation

Tan, J., Cui, W., and Zhao, L. *Bayesian spectral identification of structural damping ratio by forced-vibration exciter with non-classical damping*. Manuscript.
