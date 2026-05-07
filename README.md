#                                         CPGView (Website: http://www.1kmpg.cn/cpgview/)

The goal of the entire project is to achieve genetic visualization, The entire program needs to be executed in the Linux environment. Below are the instructions for use of  CPGView:


### 1. Installation

​	1.1  Install Pixi

Install `pixi` if necessary

```bash
# install pixi
curl -fsSL https://pixi.sh/install.sh -o install_pixi.sh
sh install_pixi.sh
```

Getting the repository including sub-modules
```bash
git clone https://github.com/michoug/CPGView
cd CPGView
```

Install the required `pixi` environments

```bash
pixi install
```

Activate the `pixi` environment

```bash
pixi shell
```


### 2. The process of implementation 

​	2.1  Execute program

   "plasdrawmap.pl" is a perl implementation file(in the  "Linux"  folder), containing the phase of three image generated.  "plasdrawmap.pl" is a simple script provided using streaming API. It takes four arguments: the plasdrawmap.pl file, the input GeneBank file(Relative or absolute path), file name/ID and the name of output folder(Relative or absolute path). Since three jobs need to be running, this process maybe take 20 -40 seconds. Below is an example:

```shell
format: perl plasdrawmap.pl example.gb example_name result_folder
example: perl plasdrawmap.pl sequence.gb Arabidopsis Arabidopsis_folder 
```

2.4  Prompt information

​			During the execution of "plasdrawmap.pl", some information will appear for users' reference, such as:

```
---  plasdrawmap completed! ---
```



### 3. The output folder

​			The output folder contains some images and some data folder. Users only need to pay attention to the files ending with "cc02.pdf", "cis.pdf" and "trans.pdf". Other files are generated in the middle for reference only.

### 4. Contact

​			Shengyu Liu - shengyuliu558@163.com 
​			Chang Liu - cliu6688@yahoo.com


