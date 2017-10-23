# Git常用命令与知识记录

### 设置使用Git的姓名和邮箱地址

```shell
$ git config --global user.name "Firstname Lastname"
$ git config --global user.email "your_email@example.com"
```

注意，如果只想对当前git目录生效而不是全局生效，去掉`--global`选项。

该命令会在`~/.gitconfig`中输出设置文件

```
[user]
name = Firstname Lastname
email = your_email@example.com
```

可以通过直接编辑该文档进行信息的修改。

### 提高命令输出的可读性

```shell
$ git config --global color.ui auto
```

“~/.gitconfig”中会增加下面一行。

```
[color]
ui = auto
```

这样一来，各种命令的输出就会变得更容易分辨。 



### 设置SSH Key

GitHub 上连接已有仓库时的认证，是通过使用了 SSH 的公开密钥认证方式进行的。现在让我们来创建公开密钥认证所需的 SSH Key，并将其添加至 GitHub。运行下面的命令创建 SSH Key。

```shell
$ ssh-keygen -t rsa -C "your_email@example.com"
Generating public/private rsa key pair.
Enter file in which to save the key
(/Users/your_user_directory/.ssh/id_rsa): 按回车键
Enter passphrase (empty for no passphrase): 输入密码
Enter same passphrase again: 再次输入密码 
```

“your_email@example.com”的部分请改成您在创建账户时用的邮
箱地址。密码需要在认证时输入，请选择复杂度高并且容易记忆的组合。
输入密码后会出现以下结果。

```shell
Your identification has been saved in /Users/your_user_directory/.ssh/id_rsa.
Your public key has been saved in /Users/your_user_directory/.ssh/id_rsa.pub.
The key fingerprint is:
fingerprint值 your_email@example.com
The key's randomart image is: 
+--[ RSA 2048]----+
| .+ + |
| = o O . |
略
```

id_rsa 文件是私有密钥， id_rsa.pub 是公开密钥。

在 GitHub 中添加公开密钥，今后就可以用私有密钥进行认证了。 



### 提交

通过`git add`命令将文件加入暂存区，然后`git commit`对文档进行提交操作。

添加成功后，可以通过`git log`命令查看提交日志。



### 进行push

只要执行push操作，GitHub上的仓库就会被更新。

```shell
$ git push
```



### git基本操作

**初始化仓库**

```shell
$ git init
```

如果初始化成功会在init的目录下生成`.git`目录，里面存放管理当前目录内容所需的仓库数据。

**查看仓库的状态**

```shell
$ git status
```

**向暂存区中添加文件** 

```shell
$　git add filename
```

**保存仓库的历史记录**

```shell
$ git commit -m "message you want to record for this commit"
```

`-m`选项后接字符串文字记述详细提交信息。

**查看提交日志**

```shell
$ git log
```

如果只想让程序显示第一行简述信息，可以在 git log命令后加上 `--pretty=short`。 只要在 `git log`命令后加上目录名，便会只显示该目录下的日志。如果加的是文件名，就会只显示与该文件相关的日志。 如果想查看提交所带来的改动，可以加上` -p`参数，文件的前后差别就会显示在提交信息之后。

**查看更改前后的差别**

`git diff`命令可以查看工作树、暂存区、最新提交之间的差别。  

不妨养成这样一个好习惯：在执行 `git commit`命令之前先执行`git diff HEAD`命令，查看本次提交与上次提交之间有什么差别，等确认完毕后再进行提交。这里的 HEAD 是指向当前分支中最新一次提交的指针。 



### 分支的操作

在进行多个并行作业时，我们会用到分支。在这类并行开发的过程中，往往同时存在多个最新代码状态。比如从 master 分支创建 feature-A 分支和 fix-B 分支后，每个分支中都拥有自己的最新代码。master 分支是 Git 默认创建的分支，因此基本上所有开发都是以这个分支为中心进行的。 

不同分支中，可以同时进行完全不同的作业。等该分支的作业完成之后再与 master 分支合并。比如 feature-A 分支的作业结束后与 master合并通过灵活运用分支，可以让多人同时高效地进行并行开发。

**显示分支一览表**

```shell
$ git branch
* master
```

可以看到 master 分支左侧标有“*”（星号），表示这是我们当前所在的分支。也就是说，我们正在 master 分支下进行开发。结果中没有显示其他分支名，表示本地仓库中只存在 master 一个分支。 

**创建、切换分支** 

如果想以当前的 master 分支为基础创建新的分支，我们需要用到`git checkout -b`命令。 

执行下面的命令，创建名为 feature-A 的分支。

```shell
$ git checkout -b feature-A
Switched to a new branch 'feature-A'
```

实际上，连续执行下面两条命令也能收到同样效果。 

```shell
$ git branch feature-A
$ git checkout feature-A
```

创建 feature-A 分支，并将当前分支切换为 feature-A 分支。这时再来查看分支列表，会显示我们处于 feature-A 分支下。 

```shell
$ git branch
* feature-A
master
```

切换回上一个分支 

```shell
$ git checkout -
```

像上面这样用“-”（连字符）代替分支名，就可以切换至上一个分支。当然，将“-”替换成 feature-A 同样可以切换到 feature-A 分支。 

**特性分支**

Git 与 Subversion（SVN）等集中型版本管理系统不同，创建分支时不需要连接中央仓库，所以能够相对轻松地创建分支。因此，当今大部分工作流程中都用到了特性（Topic）分支。 

特性分支顾名思义，是集中实现单一特性（主题），除此之外不进行任何作业的分支。在日常开发中，往往会创建数个特性分支，同时在此之外再保留一个随时可以发布软件的稳定分支。稳定分支的角色通常由 master 分支担当 。

基于特定主题的作业在特性分支中进行，主题完成后再与 master 分支合并。只要保持这样一个开发流程，就能保证 master 分支可以随时供人查看。这样一来，其他开发者也可以放心大胆地从 master 分支创建新的特性分支。 

**主干分支**

主干分支是刚才我们讲解的特性分支的原点，同时也是合并的终点。通常人们会用 master 分支作为主干分支。主干分支中并没有开发到一半的代码，可以随时供他人查看。

有时我们需要让这个主干分支总是配置在正式环境中，有时又需要用标签 Tag 等创建版本信息，同时管理多个版本发布。拥有多个版本发布时，主干分支也有多个。 

**合并分支**

我们假设 feature-A 已经实现完毕，想要将它合并到主干分支 master 中。首先切换到 master 分支。 

```shell
$ git checkout master
```

然后合并 feature-A 分支。为了在历史记录中明确记录下本次分支合并，我们需要创建合并提交。因此，在合并时加上 `--no-ff`参数。 

```shell
$ git merge --no-ff feature-A
```

随后编辑器会启动，用于录入合并提交的信息 。

默认信息中已经包含了是从 feature-A 分支合并过来的相关内容，所以可不必做任何更改。将编辑器中显示的内容保存，关闭编辑器。 