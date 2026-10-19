# OpenFoam_stratifiedFoam
基于Boussinesq假设的温盐分层流求解器
控制方程
连续性方程：
[图片]
动量方程：
[图片]
温度输运方程：
[图片]
盐度输运方程：
[图片]
状态方程：
[图片]

1、本求解器基于OpenFoam-v2106版本开发，最好在这个版本下编译，否则可能出现报错；
2、将求解器压缩包解压到任意位置，推荐放置路径为：OpenFoam-v2106/application/solvers/heatTransfer；
3、在stratifiedFoam文件夹中进入终端，加载OpenFoam环境，运行命令：wmake，无报错信息即为编译成功。后续在算例计算时，使用stratifiedFoam命令即可运行求解器开始计算。

进入文件夹后直接运行脚本./Allrun即可运行算例
注意事项：
1、分层环境初始场的设置使用的是swak4Foam工具包中的funkySetFields命令，要运行该命令需在OpenFoam中编译安装，详见下方链接：https://www.topcfd.cn/41293/。
2、由于该求解器为温盐双驱决定密度场，除了一般的边界条件外，还需要在0文件夹中制定温度T、盐度S、温度扩散率alphat、盐度扩散率kappaS的边界条件。
3、相比一般的算例，使用该求解器还需在fvSchemes和fvSolution文件中补充温度输运方程和盐度输运方程的离散格式和求解参数。标准算例中已配置好。
4、使用该求解器还需在constant/transportProperties文件中额外指定热扩散系数betaT，盐收缩系数betaS，参考温度TRef，参考盐度SRef，层流普朗特数Pr，层流施密特数Sc，湍流普朗特数Prt，湍流施密特数Sct。标准算例中已设置好推荐值。
5、标准算例由于网格数量较少，求解出的尾迹不明显，待后续网格优化后更新。
