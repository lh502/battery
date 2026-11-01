多电芯独立复现实验（2026-10-12）

先打开 Multi_Cell_Results.html 查看结果、条件和图表。
主运行入口：Run_Multi_Cell_Replication.ipynb（也可直接运行 code/run_multicell.py）。
代码保留原核心算法，适配 Arbin Excel 格式，原始输入已经包含，无须重新下载。
新增电芯：CS2-33/34（0.5C）、CS2-35/36/37/38（1C）。原CS2-21/8作为数值复核。
固定前101个候选放电循环，先编号后质控，容量61-101测试，电压81-101确认。

依次运行：
python code/run_multicell.py
python code/plot_results.py
python code/build_report.py

使用现有相容Python环境，或先运行 python -m pip install -r requirements.txt。
notebook需要另有Jupyter环境；源代码依赖只列实际分析库。

重要：本包容量结果来自电压复现链的容量诊断，不能替代正文主实验11.36→5.79。
CS2-8最终固定记录的电压曲线迁移MAE实算为97.654 mV，旧259.6 mV不适用于该记录。
完整预测、失败/排除记录、源数据哈希和绘图代码均保留。
这不是全寿命、跨化学体系或端到端零样本迁移验证。

实现修正记录：首轮检查发现Arbin放电容量列跨循环累计；现取每循环末值减初值，
并与电流积分核对。该修正未改变事先固定的0.88Ah阈值，最终全流程复跑后指标不变。
该核对的2mAh容差是读取实现核验项；其引入未新增任何循环排除。
数据计数和结论以最终保存结果为准。
