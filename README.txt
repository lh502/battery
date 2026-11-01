本次四代码重跑结果
优先打开运行结果.html，并查看本次运行说明.txt。
代码与原始四个文件逐字节相同。所有必要原始数据已保留，复跑无需联网。
protocol_before_run.json是本次从代码提取的配置记录，不是原实验协议。
安装requirements.txt依赖后，在本目录依次运行：
python code/locked_core.py --raw-file data/CS2_21_6_8_10.txt --output-dir core_results/CS2-21 --cell-label CS2-21
python code/locked_core.py --raw-file data/CS2_8_6_8_10.txt --output-dir core_results/CS2-8 --cell-label CS2-8
python code/run_multicell.py
python code/plot_results.py
python code/build_report.py
最后一句会重新生成原报告和原README，请以本次运行说明.txt作为本次实算记录。
results目录为完整数值，figures目录为四幅图，core_results为两次独立核心运行结果，logs为运行日志。
