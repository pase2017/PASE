[TOC]

## git 使用
简单介绍 git 的使用方法.

### 1. 创建仓库
初始化仓库.
``` bash
$ git init

Initialized empty Git repository in /home/XXX/pase/.git/
```

连接仓库和 github 的账号.
``` bash
$ git config --global user.name "pase2017"
$ git config --global user.email pase2017@XXX.com
```

### 2. 查看版本库当前状态
``` bash
$ git status

# On branch master
#
# Initial commit
#
# Untracked files:
#   (use "git add <file>..." to include in what will be committed)
#
#	.DS_Store
#	.pase.h.swp
#	Makefile
#	backup/
#	libpase.a
#	pase.h
#	pase_cg.c
#	pase_cg.h
#	pase_cg.o
#	pase_hypre.h
#	pase_mv.c
#	pase_mv.h
#	pase_mv.o
#	pase_parpack.c
#	pase_parpack.h
#	pase_parpack.o
#	pase_pcg.c
#	pase_pcg.h
#	pase_pcg.o
#	test/
nothing added to commit but untracked files present (use "git add" to track)
```
这是因为刚创建的代码仓库还是空的, 需要将原有的文件添加到仓库.

### 3. 添加文件到暂存区
git 有个重要的概念叫 "暂存区", 用于将当前的文件改动提交到一个暂时性的区域.

刚才 ```git status``` 命令的输出中表示所有文件都还没有被跟踪 (track), 因此我们需要将所有的文件先一次性提交.

``` bash
$ git add *
```

然后再次查看仓库状态.

``` bash
$ git status

# On branch master
#
# Initial commit
#
# Changes to be committed:
#   (use "git rm --cached <file>..." to unstage)
#
#	new file:   Makefile
#	new file:   backup/eigen-arpack.c
#	new file:   backup/eigen-parpack.c
#	new file:   libpase.a
#	new file:   pase.h
#	new file:   pase_cg.c
#	new file:   pase_cg.h
#	new file:   pase_cg.o
#	new file:   pase_hypre.h
#	new file:   pase_mv.c
#	new file:   pase_mv.h
#	new file:   pase_mv.o
#	new file:   pase_parpack.c
#	new file:   pase_parpack.h
#	new file:   pase_parpack.o
#	new file:   pase_pcg.c
#	new file:   pase_pcg.h
#	new file:   pase_pcg.o
#	new file:   test/Makefile
#	new file:   test/serPASE_ver01
#	new file:   test/serPASE_ver01.c
#	new file:   test/serPASE_ver01.o
#
# Untracked files:
#   (use "git add <file>..." to include in what will be committed)
#
#	.DS_Store
#	.pase.h.swp
```

输出 ```Changes to be committed:``` 说明改动已经添加到暂存区, 下面要做的是把暂存区的文件提交到仓库.

输出表明还有文件未被跟踪
``` bash
#	.DS_Store
#	.pase.h.swp
```
这些文件其实跟我们的代码无关, 后面在讲怎么处理.

### 4. 提交暂存区
把暂存区的更改提交.
``` bash
$ git commit -m "First commit: all files are commited to the master branch."

[master (root-commit) 279e504] First commit: all files are commited to the master branch.
 22 files changed, 3014 insertions(+), 0 deletions(-)
 create mode 100644 Makefile
 create mode 100644 backup/eigen-arpack.c
 create mode 100644 backup/eigen-parpack.c
 create mode 100644 libpase.a
 create mode 100644 pase.h
 create mode 100644 pase_cg.c
 create mode 100644 pase_cg.h
 create mode 100644 pase_cg.o
 create mode 100644 pase_hypre.h
 create mode 100644 pase_mv.c
 create mode 100644 pase_mv.h
 create mode 100644 pase_mv.o
 create mode 100644 pase_parpack.c
 create mode 100644 pase_parpack.h
 create mode 100644 pase_parpack.o
 create mode 100644 pase_pcg.c
 create mode 100644 pase_pcg.h
 create mode 100644 pase_pcg.o
 create mode 100644 test/Makefile
 create mode 100644 test/serPASE_ver01
 create mode 100644 test/serPASE_ver01.c
 create mode 100644 test/serPASE_ver01.o
```

### 5. 忽略文件
有些文件可能存在于仓库中, 但不能提交. 比如缓存文件, 密码配置文件, 目标文件等. 此时可以在仓库的根目录下创建 ```.gitignore``` 文件, 然后输入需要忽略的文件名, 以后 ```git``` 就会自动忽略.

Google 提供了一些配置好的 ```.gitignore``` 文件, 地址 [https://github.com/github/gitignore](https://github.com/github/gitignore)

我们这里以学习为目的, 于是自己手动配置.
``` bash
*.DS_Store
*.swp
*.o
*.out
*.exe
```
同样需要提交.
``` bash
$ git add .gitignore
$ git commit -m "Add .gitignore file."

[master f728e7e] Add .gitignore file.
 1 files changed, 5 insertions(+), 0 deletions(-)
 create mode 100644 .gitignore
```

此时再次查看仓库状态.
``` bash
$ git status

# On branch master
nothing to commit (working directory clean)
```

这说明所有需要的文件已经提交至仓库.

### 6. 创建分支
使用 git 进行多人协作开发时, 每个人都应该有自己的分支, 这样自己的改动只有自己能看到, 不会影响到别人. 刚才的命令输出中可以多次看到 ```On branch master```, 意思是当前位于主分支 ```matser```. 对主分支进行改动一定要慎重再慎重. 比较好的做法是, 主分支 ```matser``` 只用于大版本的更新, 分支 ```dev``` 用于小版本的更新, 个人的分支用于个人的更新. 分支 ```dev``` 还不存在, 因此需要首先创建.
``` bash
$ git branch dev
```

然后查看所有分支.
``` bash
$ git branch

dev
* master
```
目前有两个分支: ```matser``` 和 ```dev```,  其中 ```matser``` 前的星号 ```*``` 说明当前位于 ```matser``` 分支. 接下来把分支切换到 ```dev``` 上并查看.
```bash
$ git checkout dev

Switched to branch 'dev'

$ git branch

* dev
  master
```

在创建分支时, 可以同时切换到创建的分支.
``` bash
$ git checkout -b XXX

Switched to a new branch 'XXX'
```
再次查看所有分支.
``` bash
$ git branch
  dev
  master
* XXX
```

### 7. 在分支上工作
创建了属于我自己的 ```XXX``` 分支后, 就可以自己 "偷偷" 工作了, 比如创建 ```README.md``` 文件, 然后添加并提交. 查看目录下的文件.
``` bash
$ ls

backup	   Makefile   pase_cg.h  pase.h        pase_mv.c  pase_mv.o	  pase_parpack.h  pase_pcg.c  pase_pcg.o  test
libpase.a  pase_cg.c  pase_cg.o  pase_hypre.h  pase_mv.h  pase_parpack.c  pase_parpack.o  pase_pcg.h  README.md
```
可以看到多出了名为 ```README.md``` 的文件.

重新切换回 ```dev``` 分支并查看目录.
``` bash
$ git checkout dev

Switched to branch 'dev'

$ ls

backup	   Makefile   pase_cg.h  pase.h        pase_mv.c  pase_mv.o	  pase_parpack.h  pase_pcg.c  pase_pcg.o
libpase.a  pase_cg.c  pase_cg.o  pase_hypre.h  pase_mv.h  pase_parpack.c  pase_parpack.o  pase_pcg.h  test
```
新创建的 ```README.md``` 文件没有了. 这就是多人协作的好处之一, 自己的改动不会影响到别的人.

将工作合并到 ```dev``` 分支, 需要切换到 ```dev``` 分支并指定从哪个分支合并.
``` bash
$ git checkout dev

Switched to branch 'dev'

$ git merge XXX

Updating f728e7e..a7f6518
Fast-forward
 README.md |    1 +
 1 files changed, 1 insertions(+), 0 deletions(-)
 create mode 100644 README.md
```
查看当前目录.
``` bash
ls

backup	   Makefile   pase_cg.h  pase.h        pase_mv.c  pase_mv.o	  pase_parpack.h  pase_pcg.c  pase_pcg.o  test
libpase.a  pase_cg.c  pase_cg.o  pase_hypre.h  pase_mv.h  pase_parpack.c  pase_parpack.o  pase_pcg.h  README.md
```
文件 ```README.md``` 出现了.

### 8. 与 GitHub 远程仓库连接
已在 ```GitHub``` 上创建名为 ```pase``` 的仓库. 将本地的 ```master``` 分支与远程库 ```origin``` 连接起来.

首先添加远程仓库.
``` bash
$ git remote add origin git@github.com:pase2017/pase.git
```
然后推送 ```master``` 分支.
``` bash
$ git push -u origin master

Warning: Permanently added the RSA host key for IP address '192.30.255.112' to the list of known hosts.
Counting objects: 29, done.
Delta compression using up to 2 threads.
Compressing objects: 100% (28/28), done.
Writing objects: 100% (29/29), 1.82 MiB | 360 KiB/s, done.
Total 29 (delta 8), reused 0 (delta 0)
remote: Resolving deltas: 100% (8/8), done.
To git@github.com:pase2017/pase.git
 * [new branch]      master -> master
Branch master set up to track remote branch master from origin.
```

以后提交作业就可以简单地进行了.
``` bash
$ git push origin master
```

由于大家都要在 ```dev``` 分支的基础上工作, 因此该分支也需要与远程仓库连接起来. 但是个人的仓库一般不需要同步到远程, 视具体情况而定.

### 9. 从远程仓库克隆

从远程仓库克隆代码到本地.
``` bash
$ git clone git@github.com:pase2017/pase.git

Cloning into 'pase'...
remote: Counting objects: 27, done.
remote: Compressing objects: 100% (18/18), done.
remote: Total 27 (delta 6), reused 27 (delta 6), pack-reused 0
Receiving objects: 100% (27/27), 1.77 MiB | 358.00 KiB/s, done.
Resolving deltas: 100% (6/6), done.
```

将远程的 ```dev``` 分支 ```checkout``` 并在本地起名为 ```dev```.
``` bash
$ git checkout -b dev origin/dev

Branch dev set up to track remote branch dev from origin.
Switched to a new branch 'dev'
```

查看当前目录的工作区状态和分支情况.
``` bash
$ git status

On branch dev
Your branch is up-to-date with 'origin/dev'.
nothing to commit, working tree clean

$ git branch

* dev
  master
```

查看当前目录的所有文件.
``` bash
$ ls

Makefile       backup         pase_cg.c      pase_hypre.h   pase_mv.h      pase_parpack.h pase_pcg.h
README.md      pase.h         pase_cg.h      pase_mv.c      pase_parpack.c pase_pcg.c     test
```
所有文件已经到本地了.

### 10. 本地提交修改至远程仓库.

将自己的 ```SSH``` 公钥添加至 ```GitHub``` 的 ```SSH and GPG keys``` 里即可. 需要自己登录 ```GitHub```, 具体方法自行搜索即可.

提交到远程仓库很简单 (再次提醒: **慎重**对 ```master``` 分支做修改, 尽量提交至 ```dev``` 分支).
``` bash
$ git push

Counting objects: 3, done.
Delta compression using up to 4 threads.
Compressing objects: 100% (3/3), done.
Writing objects: 100% (3/3), 3.48 KiB | 0 bytes/s, done.
Total 3 (delta 1), reused 0 (delta 0)
remote: Resolving deltas: 100% (1/1), completed with 1 local object.
To pase.github.com:pase2017/pase.git
   5844e91..a9b5e84  dev -> dev
```

### 11. 使用多个 Github 账户

如果有多个 ```GitHub``` 账户 (每个账户对应一个邮箱), 不同账户是无法添加相同的 ```SSH``` 公钥的. 为了将对不同账号下仓库的修改做相应地提交, 就需要进一步配置.

首先需要生成多个公钥 (将邮箱和生成公钥的名字对应起来).
``` bash
$ ssh-keygen -t rsa -C [mail] -f [filename]
```
回车两次即可. 此时 ```~/.ssh/``` 文件夹下会出现不同的 ```.pub``` 结尾的公钥文件, 比如说 ```pase.pub``` 和 ```myamg.pub``` 对应不同 ```GitHub``` 账号下的两个仓库 ```pase``` 和 ```myamg```.

编辑文件 ```~/.ssh/config```.
``` bash
Host pase.github.com
  HostName github.com
  PreferredAuthentications publickey
  IdentityFile ~/.ssh/pase

Host myamg.github.com
  HostName github.com
  PreferredAuthentications publickey
  IdentityFile ~/.ssh/myamg
```
这里 ```Host``` 是自己取的名字, ```HostName``` 才是真正的域名.

克隆仓库时针对不同的仓库进行 ( 注意 ```@``` 之后的部分 ).
``` bash
$ git clone git@pase.github.com:pase2017/pase.git

Cloning into 'pase'...
remote: Counting objects: 27, done.
remote: Compressing objects: 100% (18/18), done.
remote: Total 27 (delta 6), reused 27 (delta 6), pack-reused 0
Receiving objects: 100% (27/27), 1.77 MiB | 267.00 KiB/s, done.
Resolving deltas: 100% (6/6), done.
```

之后就可以按照前面讲的进行开发和提交了.


### 12. 日常命令

git clone git@github.com:pase2017/pase-1.0.git
git checkout -b dev origin/dev
git pull
git add <files>
git commit -m "<commits>"
git push
git log
git reset --hard <commit id>
