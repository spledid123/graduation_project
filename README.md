这个文件的引用要[Obsidian](https://obsidian.md/)软件才能显示。
程序参考了[用 C 语言画光（一）：基础 - 知乎 (zhihu.com)](https://zhuanlan.zhihu.com/p/30745861)和[Inigo Quilez :: computer graphics, mathematics, shaders, fractals, demoscene and more (iquilezles.org)](https://iquilezles.org/articles/distfunctions2d/)

### analyse_onepartical_path.m
* 分析单粒子轨迹，输入单粒子轨迹文件
* 分为圆阵列，六边形孔喉结构两个模块，需要输入结构参数
* 功能待补充

### analyse_particals.m
* 分析多粒子轨迹表
* 分析粒子的分布，看是否稳定
* 分析速度分布

### analyse_particals_path.m
* 由路径库,赋予时间项，考虑吸附
	* 吸附时间指数分布，内置get_T_ab函数赋予
* 得到均方位移-时间，扩散系数
* 记录了吸附时间-扩散系数的数据
* 记录了粒子数-扩散系数的数据，与REV有关。

### analyse_six_pore_thermo_data.m
* 输入六边形孔喉结构的路径文件
* 分析其热力学性质(数密度/压强)与空腔的差别
* 分析差别产生原因
* 分析蒙特卡洛模拟的计算误差
* 得到孔隙度
* 依赖内置函数per_walk得到每次碰撞的数据
	* 时间项依赖内置get_T_ab函数
	* 分析位置便于区分孔依赖[scene_porethroat.m](#scene_porethroat.m函数)
* 通过times_one_particalpath得到单孔数据
* 通过内置get_th函数获得修正越孔速度方向

### arcSDF.m函数
* 圆弧的SDF
* 依赖[[#circleSDF.m函数]]

### Boltzmann.m函数
* 输出出射速度大小，满足3D/2D玻尔兹曼分布/匀速
* 3D玻尔兹曼分布就是2D泄流分布

### boxSDF.m函数
* 计算矩形的SDF
* 做了倒角，依赖[[#circleSDF.m函数]]

### change_files_name.m
* 文件名批处理

### cir.mat
* 储存开放空间内随机圆多孔介质的信息
* 圆心位置，半径
* cir_1,cir_4最后一行是孔隙度，xy范围是$-205--205$
* cir_2,cir_3是四圆

### circleSDF.m函数
* 输入：点的位置，圆心位置，圆的半径
* 输出：点到圆的距离，在圆外为正

### cube.mat
* 储存开放空间内随机矩形多孔介质的信息
* 中心点坐标x,y,半长,半宽,旋转角
* 最后一行是孔隙度，和xy范围a即$[-a,a]$

### exl_to_txt.m
 * 读取文件夹里所有xlsx文件(包括子文件夹)，将其改为txt文件
 * xlsx格式读取写入效率很慢
 * 使用readtable和writetable

### find_ReV.m
* 根据孔隙度找REV
* 找到扩散系数有统计意义需要的粒子数

### fun_exp_bol_pdf.m函数
* 画出指数分布和玻尔兹曼分布以及他们的和
* 没啥用，但算法有意思，保留

### fun_partical_one.m函数
* 单粒子运动，给定最大碰撞次数，给出初始位置与速度方向,形状参数,碰撞足够多次后停止
* 输入参数：(x0,gx0,N_pmax,filename,R,Tem),分别为粒子的初始位置，运动的初始方向,最大碰撞次数(判停条件)，输出文件名,形状参数,温度(速度参数)
* 输出参数:(filename,x,gx,s),表示文件名,末位置,末运动方向,末运动位置与原点的距离
* 每次碰撞时输出位置
* 注释参考[onegas.m](#onegas.m函数)，不考虑时间
	1. 调用[gradient.m](#gradient.m函数)找出边界法线
	2. 调用[reflect.m](#reflect.m函数)输出每次反射的速度方向
	3. 调用[scene.m](#scene.m函数)输出与边界的最近距离

### get_cubecir.m
* 从在圆形阵列中单粒子轨迹得到有关参数构建随机游走模型
* 因为圆形阵列孔喉几何结构划分不明显，所以没办法由几何得到参数（评价为这种做法没啥用）
* 待改进

### get_D_porethroat_randomwalk.m
* 参考说明文档中最后的算法即可

### get_frac_cir.m
* 得到并画出满足分形规律的圆阵列
* 调用内置函数get_fractal_cir_R得到圆半径

### get_particals_D.m
* 由多粒子的不同时刻的位置，得到扩散系数
* 计算不同方向的扩散系数
* 扩散系数与圆阵列R的关系

### get_porethroat_para.m
* 得到六边形孔喉模型的参数，孔喉的单步时间/概率；平均单步时间

###  get_poro.m函数
* 输入：图片gca,已经画好的多孔介质图，固体用一种颜色填充，'color'为RGB
* 输出孔隙度

### get_poro_six.m函数
* 输入六边形孔喉结构的参数，返回孔隙度

### get_randomwalk_D.m
* 根据孔喉双重结构随机游走模型的参数计算扩散系数
* 参数计算依赖[[#get_porethroat_para.m]]

### gradient.m函数
 * 输入边界点的位置，输出边界法线的单位向量

### make_random_v.m
*  生成随机数,满足3D/2D玻尔兹曼分布,储存在相应excel文件中
* 需要更改[Boltzmann.m](#Boltzmann.m函数)文件更改分布

### matlab.mat
* 储存开放空间内多孔介质稳定后的数密度/粒子数信息
* 归一化：$r=500,\Delta T=1/2000$(大圆半径，射入频率)
* 具体信息见[[#plot_Box.m]]

### Ngas.m函数
* 功能与[onegas.m](#onegas.m函数)类似，但是砍掉了记录轨迹的功能

### onegas.m函数
 * 输入位置与速度，以及时间间隔T0，粒子运动T0期间，每次在边界脱附的形成位置表，并输出结束时的位置与速度。考虑每次在边界的吸附，赋予停留时间。
	1. 调用[gradient.m](#gradient.m函数)找出边界法线
	2. 调用[reflect.m](#reflect.m函数)输出每次反射的速度方向
	3. 调用[Boltzmann.m](#Boltzmann.m函数)输出每次反射的速度大小
	4. 调用[T_absorb.m](#T_absorb.m函数)输出每次碰撞的吸附时间
	5. 调用[scene.m](#scene.m函数)输出与边界的最近距离

### partical_N.m
* 描述N个粒子的运动，赋予速度，固定时间间隔记录粒子位置，速度
* 粒子并行计算
* 考虑周期性边界
* 调用[Ngas.m](#Ngas.m函数)对单粒子在某个时间间隔的轨迹进行计算
* 可以直接读取纪录粒子信息的表格文件，直接开始模拟

### partical_one.m
 * 对于单个粒子，给出初始位置与速度方向,形状参数,得到初始时间T=0到sum_j$\times$N_T时间内，每次碰撞的位置表(时间,坐标x 坐标y 与边界的距离)，调用[onegas.m](#onegas.m函数)函数

### partical_one_hete.m
* 单粒子在非均匀圆阵列中路径，半径均匀分布
* 调用内置函数fun_partical_one_hete输出轨迹
* fun_partical_one_feijunyun函数
	* 给定时间间隔，让粒子在非均匀随机生成半径的圆形阵列里运动，边运动边生成半径。
	* 调用内置函数scene_hete，gradient_hete;
		* 半径的生成调用内置函数get_ran_r.m
	* 调用[reflect.m](#reflect.m函数)
* 全局变量记录已经生成的半径

### partical_one_para.m
* 并行地得到一系列粒子的轨迹，只输出碰撞时的位置
* 需要给出每个粒子的形状参数，给出粒子数，最大碰撞次数
* 对于每一个粒子,调用[fun_partical_one.m](#fun_partical_one.m函数)计算

### paths_1D.m
* 模拟多粒子在2维直管的运动
* 考虑吸附时间占主导与运动时间为主的情况
* 可以进一步修改结合二者

### pipe_entrance_effect.m
* 直管"入口效应"与单孔长度的关系
* 从单孔碰撞次数/是否回到原方向两方面去看

### plot_Box.m
* 配套[[#analyse_particals.m]]中关于开放空间内多孔介质内外粒子数+数密度的分析
* 将以下语句分别贴到命令行窗口，把数密度信息存入matlab.mat
* 画箱型图


### plot_cir.m函数
* 给定圆心半径，画圆
* 能够选择是否填充以及颜色；
* 可选参数'flag',flag == 1表示填充
* 名称数值对：'var'为plot的可选参数；'fillvar'为填充颜色,'linewidth'为线宽

### plot_cross_boundary
* 画出某次跨越边界的轨迹
* 需要补充，修改

### plot_cube_porethroat
* 画出正方体排列的孔喉阵列

### plot_distribution.m函数
* 输入：提供数组kTTT；flag=0为连续，flag=1为离散整数
* 若为连续，则区间取edges个;
* 输出图：连续的纵坐标为频率密度,离散整数时为频率
* 输出：plot(xxfit,yy)可以画出对应频率密度分布折线图

### plot_location_probability.m
* 根据某个时刻多粒子的位置表(或者单粒子固定时间间隔的位置)，画出二维的概率密度图

### plot_partical_path.m
* 画多孔介质：圆形阵列；六边形孔喉
	* 调用[plot_cir.m](#plot_cir.m函数)画圆
	* 调用内置plot_onecir_in_six和plot_pore_cir两个函数
	* plot_onecir_in_six函数画出六边形孔喉结构的单个圆及三个半喉
* 根据轨迹文件画出例子轨迹
	* 内置get_T函数根据轨迹位置赋予时间,可以调整吸附时间
	* 能够根据时间给轨迹附上颜色
* 其它一些多孔介质模型

### plot_particals.m
* 根据多粒子位置，画粒子

### plot_particals_gif.m
* 画出不同时刻粒子位置以及多孔介质,画gif

### plot_randomwalk.m
* 画示意图——六边形孔喉随机游走模型

### plot_rec.m函数
* 画长方形，可选是否填充，以及颜色；
* 输入：矩形中心点，半长半宽，旋转角
* 可选参数'flag',flag == 1表示填充
* 名称数值对：'var'为plot的可选参数；'fillvar'为填充颜色

### plot_spe.m函数
* 根据函数$f(t)=x$,t为等距采样，为等差数列,画出频谱图，找周期
* 输入为$t,x$一一对应的采样点，输出为频谱图
* 函数内做了$x$的平均化$x=x-mean(x)$,所以消去了0处峰值

### plot_v_cloud_map.m
* 根据某时刻的粒子位置速度，画速度云图，可能会有归一化
* 画需要的多孔介质

### randomwalk_1D.m
* 模拟多粒子的1维随机游走，考虑吸附
* 画出粒子在每一步的吸附概率与扩散系数的关系
* 线性关系

### randomwalk_porethroat_network.m
* 通过随机游走模拟得到单元时间概率分布，单步时间可设为相同的或不同的
* 计算概率分布，得到单元停留时间/次数pdf
* 输入参数：当前在孔，某一步在原地的概率为Pp;当前在喉，某一步留在原地的概率为Pt;在孔内的单步时间为Tp,在喉内的单步时间为Tt

### reflect.m函数
 * 输入射入速度方向与法线方向，输出出射方向

### scene.m函数
 * 输入位置与结构参数，输出点到边界的最短距离
 * 返回点到边界的最短距离
	1. 圆形阵列，R=(Rx,Ry)为xy方向的圆心间距/直径，调用[circleSDF.m](#circleSDF.m函数)
	2. 孔喉六边形结构，R(1)为喉道长度/直径；R(2)为喉道宽度/直径，找出点所在的基本单位（两圆+与其相连的喉道），并调用[twocircleSDF.m](#circleSDF.m函数)算距离
	3. 周期性空腔圆阵列交错整个周期性边界大致为九宫格：存在两种相同的矩形区块，分别有不同的形态，任意一种的上下左右都是另一种，一个区块是圆阵列为numx$\times$numy,另一个区块是空白;$R = [numx,numy]$
	4. 周期性边界,为九宫格：存在两种边长相等正方形区块，分别有不同的形态，任意一种的上下左右都是另一种；为了对称性，先暂定两种区块都是均质圆阵列且xy方向性质一样；多孔介质参数$R=[右边的圆半径 右边的圆半径+喉道宽度一半 右边的一排/列圆数量 左边的圆半径 左边的圆半径+喉道宽度一半 左边的一排/列圆数量]$
	5. 周期性边界，为九宫格：存在两种边长相等正方形区块，分别有不同的形态，任意一种的上下左右都是另一种，且之间存在空余；为了对称性，先暂定两种区块都是均质圆阵列且xy方向性质一样。多孔介质参数$R=[右边的圆半径 右边的圆半径+喉道宽度一半 右边的一排/列圆数量 左边的圆半径 左边的圆半径+喉道宽度一半 左边的一排/列圆数量 空出的长度]$
	6. bulk-porous media体系,bulk为一个圆，粒子以某个频率向中间发射，中间是多孔介质的体系，也算作方形圆阵列。$R = [r,Rp,nump]$,r为bulk圆半径,圆心在原点,Rp为阵列的孔喉比,nump为圆阵列每行个数，圆的半径为1
	7. 圆形空腔,R表示半径

### scene_porethroat.m函数
 * 输入点的位置x；R为喉道长度/直径；a为喉道宽度/直径
 * 输出sd，sd(1):返回点到六边形孔喉结构边界的最短距离,sd(2),sd(3)为点最靠近的圆心的阵列编号(x,y方向);sd(4)表示其是否为孔，是为1，否为0

### segmentSDF.m函数
 * 计算线段AB的SDF

### solve_pde_1Ddiffusion
*  求解1D的扩散偏微分方程

### T_absorb.m函数
* 产生吸附时间随机数，满足指数分布/相同时间

### triangleSDF.m函数
* 计算三角形的SDF
* 依赖[[#segmentSDF.m函数]]

### twocircleSDF.m函数
 * 输入：点的位置，其所在基本单元：两个圆+与其相连的喉道，两个圆的圆心，喉道的长度/宽度
 * 输出：sd:返回点到结构单元结构边界的最短距离,kx,ky为点最靠近的圆心的阵列编号(x,y方向);ispore表示其是否为孔，是为1，否为0

