# 电力电子化电力系统分析工具——阻抗法 模型及技术文档

本项目是我个人整理的阻抗法相关论文的模型算例及技术文档，在遵守 GPL v3.0 协议的情况下可自由使用，详见[License](https://github.com/tommyjiang/Impedance/blob/master/LICENSE).

# 项目说明

## Document
本项目的技术文档，在 Mac OS Sierra(10.12.5) + TeXLive 2017 下利用 XeLaTeX 编译即可，注意需要 Adobe 的四套字体（宋体、仿宋、黑体、楷体），详见 [GitHub Font Repo](https://github.com/dolbydu/font/tree/master/unicode).

## Model
本项目的模型文件，环境为 MATLAB 2015b。
- 各文件说明（前面均为作者姓+年份+期刊+期号）：
  - *.slx：模型
    - *_m.m：参数
	  - *_sim.m：仿真脚本

	  相关问题：
	  - 编码：由于 Mac 下 MATLAB 使用 UTF8 编码，而 Windows 下为 GB2312 编码，为避免出现乱码问题，程序和注释全部为英文。
	  - Simulink 模型：slx 文件，由于模型 slx 文件为二进制文件，目前不接受 Pull Request 修改模型，可以开 Issue 讨论。
