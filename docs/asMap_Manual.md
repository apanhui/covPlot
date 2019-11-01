## asMap使用手册

#### 介绍

asMap是用来绘制基因的深度分布和splice reads分布情况的一个画图程序。
asMap的输入是bam文件和gtf文件，一次只能绘制一个基因的图。
asMap支持多样本，同时会展示一个基因的所有转录本（根据gtf文件自动识别）。
例图如下：

![asMap例图](E:/画图/可变剪切reads分布图.png)


#### 使用方法

asMap主要通过配置文件来控制输入和图形的样式，为了方便使用部分重要的参数可以通过命令行直接设置。使用方法如下(可以通过perl asMap.pl --help唤起)：

```
Usage:   perl asMap.pl [options] --conf your.conf [ --gene target_gene ]

Import Options:
         --help                      print the simple usage info
         --conf   <FILE>             the configuration file, ["/Bio/User/aipeng/bin/asMap/asMap.conf"]
         --gene   <STR>              set the gene name which you want to draw (must be in the gtf file)

Input Options:
         --gtf    <FILE>             set the gtf file

Styles Options (it's better to be set in conf file):
         --intron_fix         <INT>      fix the intron sizse as [200] nt
         --intron_scale       <FLOAT>    zoom the intron size as [auto] ratio to raw size, [scale the introns total size == exons]
         --intron_background             display the intron backgroud or not
         --width              <INT>      the width size of the figure, [1000]
         --from_top                      the connect junction lines draw from top first, otherwise will from bottom
         --unit_depth_height  <INT>      the unit depth height of each sample, [100]

```

其他的参数详见配置文件的说明。

使用注意事项：

1. 推荐大家每个项目填写一个配置文件，调试图片样式的时候可以通过部分参数来实现快速调试，但是在固定之后还是通过配置文件来固定。然后--gene这个参数可以帮助大家通过命令行来实现批量绘制基因。
2. 图片的大小主要通过width和unit_depth_heigth来设定，其中unit_depth_height是指每个样本的高度，会随着样本数量的变化而变化。
3. intron_fix和intron_scale这两个参数主要影响的是外显子和内含子的显示比例，由于默认情况下外显子的长度都远远小于内含子的长度，导致如果真实的比例来显示的话，外显子就很难看清楚，而reads的分布恰恰主要是在外显子区，所以主要通过这2个参数来控制。intron_fix是固定所有内含子的长度为某个值，intron_scale是将内含子的长度按照固定的比例来缩小或者放大，默认将内含子的长度均一化成与外显子的长度一致）。不确定的时候建议使用默认的参数。
4. 配置文件中设置颜色时，如果是"#"开头的十六进制RGB颜色，一定要记得在"#"前加上反义符“\”。
5. 关于marker和legend的功能目前还不是很完善，颜色和大小以及支持的类型都在代码中定义死了，如果有需要的话直接与我联系。（大家在跑的时候记得把legend设置成no）
6. 目前的输入文件除了包含基因的gtf文件外，还需要输入每个样本的比对结果和junction.bed文件，其中比对结果是用来判定reads深度的，junctions.bed文件是用来读取splice read位置和深度信息的。目前仅支持tophat2的结果，其他软件的结果没有测试。