# Neoantigenproject overview
**This is a Neoantigen database project, Used for data preprocessing**

## 1. Introduction of folder and file:
+ /script is a folder with all .R scripts for data preprocessing and data statistic.
+ /data is a folder with all data related with your package(attention:This is for your package not for your projects, It could be generated automatically).
+ /data-raw is a folder with all datasets you would like to put them in your package
+/data_project is a folder with all data related with you project.
+ /R is for All function I would like to create a R package.
+ /man is a document folder for R function, It could be generated by document() function automatically.
+ DESCRIPTION is a file where contain your package meta information.


## 2. create or update your R package:
1. Change R function under /R folder.
2. Add Documentation comment in front of R code.  

![](./Attachment/Rcomment.png)  

3. Update DESCRIPTION file if their have a need.
4. library(devtools)
5. document() #update all document file.
6. build() # Build all R function again.
7. You can upload your new R package now.  

## 3. Reference(How to Build a Project and make a package together)

1. create a new project use Rstudio version control system.
[restudio-git-github](https://happygitwithr.com/rstudio-git-github.html)
2. Copy DESCRIPTION and NAMESPACE model files under this folder.
[make a Rpackage](https://www.davekleinschmidt.com/r-packages/)
3. mkdir data_project/script folder(This is for your project,all folder or file related with making package could be generate automatically).  

![](./Attachment/projectfile.png)
4. Add all file you need and update with git button.  

![](./Attachment/Rstudiogit.png)

5. You can also make sample data in your package.([Add data in your Rpackage](https://www.davekleinschmidt.com/r-packages/))

**Attention:** Please don't forget to git push after you add your commit information.  

**Contact:** If their are any question, please file free to contact with author:jijunyu140@gmail.com
