#!/usr/bin/env python3
"""Build the self-contained Chinese report and compact tables from saved results."""
from pathlib import Path
import base64
import html
import json
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]

def table(frame):
    return frame.to_html(index=False,escape=False,border=0,classes='data',float_format=lambda x:f'{x:.3f}')

def main():
    s=pd.read_csv(ROOT/'results/multicell_summary.csv')
    new=s[~s.cell.isin(['CS2-21','CS2-8'])].copy()
    capacity=new[['cell','capacity_test_cycles','recent_state_capacity_MAE_mAh','linear_capacity_MAE_mAh','saturating_capacity_MAE_mAh','capacity_MAE_reduction_percent']].copy()
    capacity.columns=['电芯','测试循环数','近期基线 MAE / mAh','线性模型 MAE / mAh','饱和模型 MAE / mAh','相对近期基线降低 / %']
    voltage=s[['cell','confirmation_cycles','pooled_voltage_MAE_mV','pooled_voltage_RMSE_mV','pooled_voltage_bias_mV','CS2_21_transfer_MAE_mV']].copy()
    voltage.columns=['电芯','确认循环数','电压 MAE / mV','电压 RMSE / mV','偏差 / mV','CS2-21 曲线迁移 MAE / mV']
    capacity.to_csv(ROOT/'Capacity_Results.csv',index=False,encoding='utf-8-sig')
    voltage.to_csv(ROOT/'Voltage_Results.csv',index=False,encoding='utf-8-sig')
    ci=new[['cell','cycle_mean_voltage_MAE_mV','cycle_mean_MAE_CI95_low_mV','cycle_mean_MAE_CI95_high_mV','capacity_paired_gain_CI95_low_mAh','capacity_paired_gain_CI95_high_mAh']].copy()
    ci.columns=['电芯','循环均值 MAE / mV','循环均值 CI 下限 / mV','循环均值 CI 上限 / mV','容量误差改善 CI 下限 / mAh','容量误差改善 CI 上限 / mAh']
    endpoints=new[['cell','confirmation_cycles','endpoint_evaluable_count','endpoint_MAE_min_diagnostic']].copy()
    endpoints.columns=['电芯','确认循环数','可评价终点数','终点 MAE / min（仅诊断）']
    records=[]
    for cell in new.cell:
        q=pd.read_csv(ROOT/'results'/cell/'cycle_quality_control.csv')
        records.append({'电芯':cell,'筛选前循环':101,'有效循环':int(q.accepted.sum()),
                        '排除循环':','.join(q.loc[~q.accepted,'cycle'].astype(str)) or '无',
                        '起始日期':str(q.date_start.iloc[0])[:10],
                        '第101次记录结束日期':str(q.date_end.iloc[-1])[:10]})
    figures=[]
    captions=[
      '图1　新增六个电芯的容量模型比较，以及八个电芯的循环平均电压误差与95%区间。右图是循环均值，和表内合并采样点的 MAE 统计口径不同。',
      '图2　电芯独立标定与 CS2-21 电压曲线直接迁移。两种设置均保留目标电芯自己的容量模型；不是整套模型的零样本迁移。',
      '图3　新增六个电芯的容量实测轨迹与冻结后的预测。灰色为1–60开发区间。',
      '图4　新增电芯的实测/重构端电压轨迹。展示最接近第91次的有效循环，以循环位置选择，未按误差选择。']
    for f,caption in zip(sorted((ROOT/'figures').glob('*.png')),captions):
        data=base64.b64encode(f.read_bytes()).decode('ascii')
        figures.append(f'<figure><img src="data:image/png;base64,{data}" alt="{html.escape(caption)}"><figcaption>{caption}</figcaption></figure>')
    body=f'''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>多电芯独立复现实验结果</title><style>
body{{font-family:-apple-system,BlinkMacSystemFont,"PingFang SC","Microsoft YaHei",sans-serif;color:#263640;line-height:1.75;background:#f2f5f7;margin:0}}
main{{max-width:1150px;margin:28px auto;background:white;padding:42px 48px;border-radius:12px}}h1{{font-size:29px;line-height:1.3;margin:0 0 12px}}h2{{font-size:21px;margin:34px 0 12px;color:#245F91}}p{{margin:12px 0}}.meta{{color:#677c89;font-size:14px}}.note{{border-left:4px solid #245F91;padding:12px 18px;background:#edf4f8}}.caution{{border-left:4px solid #c76b47;padding:12px 18px;background:#fbf3ed}}table{{border-collapse:collapse;width:100%;font-size:13px}}th{{background:#eaf0f4;text-align:left}}td,th{{padding:9px 10px;border-bottom:1px solid #dce5eb}}td:not(:first-child){{font-variant-numeric:tabular-nums}}.scroll{{overflow-x:auto}}figure{{margin:28px 0}}img{{max-width:100%;height:auto}}figcaption{{color:#536773;font-size:13px}}code{{background:#eef2f4;padding:2px 5px}}li{{margin:7px 0}}a{{color:#245F91}}@media(max-width:700px){{main{{padding:24px 18px;margin:0}}}}@media print{{body{{background:white}}main{{padding:0;margin:0}}figure,table{{break-inside:avoid}}}}
</style><main><h1>多电芯独立复现实验结果</h1><div class="meta">执行日期：2026-10-12　|　CALCE CS2 系列　|　6个新增电芯 + 2个基准复核</div>
<p class="note"><b>已实际运行：</b>新增 CS2-33、CS2-34（0.5C）以及 CS2-35、CS2-36、CS2-37、CS2-38（1C）。六个新增电芯的电压MAE为8.858–62.678 mV；容量MAE均优于近期基线，但CS2-37未优于线性基线。全部结果均保留。结论适用于本次同系列电芯与固定101循环记录窗口。</p>
<h2>1. 代码是否需要重写</h2><p>无需重写核心算法。保留最终包中的容量拟合、电压–SOC关系构建和评价定义，增加 Arbin Excel 读取、累计容量增量处理、局部循环号的时间排序、倍率适配、批量运行和结果汇总。原始算法原样保存在 <code>code/locked_core.py</code>。</p>
<p>原有 CADEX 文件仍走原版读取器。新增 Arbin 文件按文件日期读取、按实际放电起始时间排列，先编号再质控；不因排除循环而重新压缩编号。首批新增名单、101循环上限和评价区间在新模型结果产生前固定。收到继续增加电芯的要求后，另行固定CS2-37/38名单，沿用全部规则，没有根据前四个电芯的结果调参。</p>
<h2>2. 固定实验设计与数据来源</h2><ul>
<li>1–40：初始拟合；41–60：开发验证；1–60重新拟合后冻结。61–101：容量测试；61–80：电压稳健性检查；81–101：电压确认。</li>
<li>每个目标电芯用自己的开发数据标定容量模型和电压曲线；测试容量没有作为主预测输入，测试电压没有参与拟合。</li>
<li>统一取首次3.0 V穿越并线性插值，积分实测电流得到相应容量。新文件通常实际放电至2.7 V；这里3.0 V是统一分析终点。</li>
<li>新电芯质控：至少80个放电采样点；时长至少0.75/倍率小时（0.5C为1.5 h，1C为0.75 h）；每次放电容量增量至少0.88 Ah；必须穿越3.0 V；电流变异系数≤5%，中位电流偏离名义值≤10%。</li>
<li>Arbin容量计数可能跨循环累积，故用本循环末值减初值。增量与电流积分的允许差为2 mAh，实际最大差小于0.045 mAh。该读取修正后全流程复跑，主要结果不变。</li>
<li>饱和模型沿用原实现：C∞在700–1000 mAh内、幅度非负，τ在3–500循环的1500点对数网格上选取。近期基线为开发末10个有效循环的容量中位数；线性基线用同一1–60开发集拟合。</li></ul>
<div class="scroll">{table(pd.DataFrame(records))}</div>
<p>排除原因为不完整放电或未达到3.0 V等预设质控条件；详见每个电芯的 <code>cycle_quality_control.csv</code>。实际循环记录来自2010年；本次执行时间不代表采集时间。原 CS2-21、CS2-8 使用固定的2010-06-08记录，其局部循环起点不等于新 Arbin 数据的寿命起点。</p>
<p>官方来源：<a href="https://calce.umd.edu/battery-data">CALCE Battery Data</a>；新增原始压缩包为 <a href="https://web.calce.umd.edu/batteries/data/CS2_33.zip">CS2_33.zip</a>、<a href="https://web.calce.umd.edu/batteries/data/CS2_34.zip">CS2_34.zip</a>、<a href="https://web.calce.umd.edu/batteries/data/CS2_35.zip">CS2_35.zip</a>、<a href="https://web.calce.umd.edu/batteries/data/CS2_36.zip">CS2_36.zip</a>、<a href="https://web.calce.umd.edu/batteries/data/CS2_37.zip">CS2_37.zip</a>、<a href="https://web.calce.umd.edu/batteries/data/CS2_38.zip">CS2_38.zip</a>。包内保存本次实际读取的原始工作簿、校验值和来源清单，不需要联网即可复跑。</p>
<h2>3. 新增容量测试结果</h2><div class="scroll">{table(capacity)}</div>
<p>容量单位为mAh；相对降低=(近期基线MAE−饱和模型MAE)/近期基线MAE。这里评价的是电压复现链中的“积分至3.0 V容量”诊断；不能替换论文原主容量实验的11.36→5.79 mAh。</p>
<h2>4. 电压确认结果与直接迁移</h2><div class="scroll">{table(voltage)}</div>
<p>电压主指标为各循环共同支持区间上的采样点合并统计；偏差=预测−实测。共同支持止于实测3.0 V时刻和预测容量耗尽时刻中较早者；各循环覆盖率保存在明细中。CS2-21自身一行的“迁移”数只是同一电压曲线的恒等检查。</p>
<p class="caution"><b>旧迁移数字需更正：</b>CS2-21曲线用于最终固定 <code>CS2_8_6_8_10.txt</code>，本次实算MAE为 <b>97.654 mV</b>。259.6 mV属于旧记录，不能写成当前固定记录的迁移结果。转移的是电压关系，目标电芯的容量模型仍独立标定；不同测试仪和1C工况也可能影响迁移误差，不能把误差全部归因于电芯个体差异。</p>
<h2>5. 不确定性与终点指标</h2><div class="scroll">{table(ci)}</div>
<p>循环均值区间通过循环级循环移动块bootstrap计算，20,000次、块长5、种子20260806。容量改善区间基于同一循环上“基线绝对误差−饱和模型绝对误差”。这些是给定电芯和窗口内的不确定性，不是电芯总体泛化置信区间；表内循环均值和上表合并点MAE不得混称。</p>
<div class="scroll">{table(endpoints)}</div>
<p>电压关系强制U(0)=3.0 V，预测终点与容量耗尽耦合；终点指标仅作一致性诊断，且只包含记录支持范围内可达的终点，不能作为独立电压触发关机或真实手机TTE验证。</p>
<h2>6. 图表</h2>{''.join(figures)}
<h2>7. 可以支持的结论与仍需保留的限制</h2><ul>
<li>这次增加了6个独立电芯的流程复现，覆盖同系列下0.5C与1C放电；容量模型是否优于基线需按电芯逐一判断，不把负结果排除。</li>
<li>CS2-37的饱和模型容量MAE为20.198 mAh，线性模型为15.294 mAh；CS2-38的电压MAE为62.678 mV，明显高于其余新增电芯。该电芯替换为实测容量的仅评价诊断仍有59.183 mV误差，提示容量预测误差不能解释全部电压失配，不能声称冻结电压关系在所有电芯上始终保持同等精度。</li><li>原CS2-21与CS2-8电压结果被数值复现；优化后的批量评价与原评价函数在新增Arbin循环上也核对一致，CS2-8迁移结果另用原函数完整复核。</li>
<li>这仍是CS2系列、有限记录窗口和电芯独立重新标定条件下的证据。不能写成跨化学体系通用模型、总体统计泛化、全寿命预测或端到端零样本迁移成功。</li>
<li>6个新增电芯的C∞拟合值均靠近原1000 mAh上界，因此不宜赋予其已经验证的长期物理饱和值含义。未据测试结果放宽边界或改拟合方法。</li>
<li>本次容量诊断中CS2-21饱和模型MAE为10.667 mAh，近期基线为9.889 mAh；CS2-8线性模型MAE为4.559 mAh，略低于饱和模型的4.931 mAh。完整结果均保留，不能宣称饱和模型在所有电芯、所有基线上最优。</li></ul>
<h2>8. 复跑与交付文件</h2><p>打开 <code>Run_Multi_Cell_Replication.ipynb</code> 从上到下运行，或在本文件夹运行：</p><pre>python code/run_multicell.py
python code/plot_results.py
python code/build_report.py</pre><p>依赖版本在 <code>requirements.txt</code>。四幅图均提供PNG、PDF和SVG；完整数值在 <code>results/multicell_summary.csv</code>，逐循环预测、拟合参数和质控记录在各电芯子目录。SHA256SUMS用于核对文件完整性。</p>
</main></html>'''
    (ROOT/'Multi_Cell_Results.html').write_text(body,encoding='utf-8')
    readme='''多电芯独立复现实验（2026-10-12）

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
'''
    (ROOT/'README.txt').write_text(readme,encoding='utf-8')

if __name__=='__main__':main()
