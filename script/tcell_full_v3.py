#######################解析tcell_full件，获得入库前文件（ligand、mhcI和II）
#制作含有所需行的小文件
tcell_full = open("Dataset_v2/tcell_full_v3.csv")
head1 = pd.DataFrame(tcell_full.readline().split(","))
head2 =  pd.DataFrame(tcell_full.readline().split(","))
head = head1+"_"+head2
head = head[0].tolist()
tcell_full = pd.read_csv("Dataset_v2/tcell_full_v3.csv",skiprows= 2)
tcell_full.columns=head
print(len(tcell_full))
select_col=["MHC_Allele Name","Epitope_Description","Epitope_Starting Position","Epitope_Ending Position"]
select_col.extend(["Assay Antigen_Antigen Epitope Relation","Assay_Assay Group"])
select_col.extend(["Epitope_Organism Name","Host_Name","Epitope_Antigen Name","Epitope_Antigen IRI","Reference_PubMed ID"])
select_col.extend(["Assay_Units","Assay_Qualitative Measure","Assay_Measurement Inequality","Assay_Quantitative measurement"])
tcell_full=tcell_full[select_col]
tcell_full.to_csv("Dataset_v2/tcell_full_v3_select_col.csv",index=False)


tcell_full = pd.read_csv("Dataset_v2/tcell_full_v3_select_col.csv")
print len(tcell_full)
print tcell_full.columns
#362985
#Index([u'Peptide', u'Start_Position_in_Antigen', u'Epitope_Organism',
#      u'Antigen_Name', u'Epitope_Antigen_IRI', u'Reference', u'Assay_Units',
#       u'Qualitative_Measure', u'Measurement_Inequality',
#       u'Quantitative_measurement', u'MHC_class', u'protein_id', u'Species',
#       u'MHC_Allele', u'MHC_Type', u'Length_of_Peptide'],
#      dtype='object')

#确定如何保留物种
#去除了9万条
tcell_full=tcell_full.dropna(subset=["MHC_Allele Name"])
mhc=tcell_full["MHC_Allele Name"].drop_duplicates().tolist()
print len(mhc)

mhc_1 = mhc[:]
for i in mhc:
    if i.startswith("HLA") or i.startswith("H2") or i.startswith("RT1"):
        mhc_1.remove(i)
print mhc_1

print len(tcell_full) #272825
#通过MHC的拼写选取三个物种
tcell_full=pd.concat([tcell_full[tcell_full["MHC_Allele Name"].str.startswith("HLA")],
                           tcell_full[tcell_full["MHC_Allele Name"].str.startswith("H2")] ,
                           tcell_full[tcell_full["MHC_Allele Name"].str.startswith("RT1")]],)
len(tcell_full)#260049

#页面展示prtein id
tcell_full["protein_id"]=tcell_full["Epitope_Antigen IRI"].str.split("/",expand=True).iloc[:,-1]
#物种根据MHC的拼写规则定义
tcell_full["Species"]="Homo sapiens"
tcell_full.loc[tcell_full["MHC_Allele Name"].str.startswith("H2"),"Species"] ="Mus musculus"
tcell_full.loc[tcell_full["MHC_Allele Name"].str.startswith("RT1"),"Species"] ="Rattus norvegicus"

#定义MHC类型及修正亚类的书写格式
mhc_type = pd.DataFrame(tcell_full["MHC_Allele Name"].drop_duplicates())
mhc_type.index=range(len(mhc_type))
for i in range(len(mhc_type)):
    i_i=mhc_type.loc[i,"MHC_Allele Name"]
    #共有127万数据，有HLA I类45万，HLA II类1万,RT11 class II 2个, 未给出明确亚型
    if i_i=="HLA class I":
        mhc_type.loc[i,"MHC_Allele"] = "HLA class I"
        mhc_type.loc[i,"MHC_Type"] = "HLA class I(Undetermined type)"
    elif i_i=="HLA class II":
        mhc_type.loc[i,"MHC_Allele"] = "HLA class II"
        mhc_type.loc[i,"MHC_Type"] = "HLA class II(Undetermined type)"
    elif i_i=="RT11 class II":
        mhc_type.loc[i,"MHC_Allele"] = "RT11 class II"
        mhc_type.loc[i,"MHC_Type"] = "RT11 class II(Undetermined type)"
    elif i_i.startswith("HLA-"):
        if "*" not in i_i:
           #打印出所有情况，目前存在后两位可能存在数字或无数字，故为字母和数字间加*，其余原来有*的不做改动
            try:
                int(i_i[-2])
            except:
                try:
                    int(i_i[-1])
                except:
                    pass
                else:
                    i_i=i_i[:-1]+"*"+i_i[-1:]
            else:
                i_i=i_i[:-2]+"*"+i_i[-2:]
        mhc_type.loc[i,"MHC_Allele"] = i_i
        mhc_type.loc[i,"MHC_Type"] = i_i[:5]
    elif i_i.startswith("RT1"):
        mhc_type.loc[i,"MHC_Allele"] = i_i
        mhc_type.loc[i,"MHC_Type"] = i_i[:5]
    elif i_i.startswith("H2"):
        #H2类型不太好分，需要于老师后期确定后再进行修改
        mhc_type.loc[i,"MHC_Allele"] = i_i
        mhc_type.loc[i,"MHC_Type"] = "H2"
    else:
        print(i_i)

tcell_full=pd.merge(tcell_full,mhc_type,on="MHC_Allele Name",how="left")


tcell_full["Length_of_Peptide"] = tcell_full["Epitope_Ending Position"]-tcell_full["Epitope_Starting Position"]

#通过蛋白名称获得基因symbol
#该文件通过url去重后共有119347中链接，三类分别为ncbi、uniprot和iedb的ontology，前两种通过网页均能查到基因名称，后一种
#该问题暂放，看是否有其他更佳方式获得基因symbol
url = pd.DataFrame(tcell_full["Epitope_Antigen IRI"].dropna().drop_duplicates())
url.index=range(len(url))  

#写入新增的列["Assay Antigen_Antigen Epitope Relation","Assay_Assay Group"]






#重新规范定义列名,去除不必要的列
tcell_full.rename(columns={'Epitope_Description':'Peptide'},inplace=True)
tcell_full.rename(columns={'Reference_PubMed ID':'Reference'},inplace=True)
tcell_full.rename(columns={'Epitope_Starting Position':'Start_Position_in_Antigen'},inplace=True)
tcell_full.rename(columns={'Epitope_Organism Name':'Epitope_Organism'},inplace=True)
tcell_full.rename(columns={'Epitope_Antigen Name':'Antigen_Name'},inplace=True)
tcell_full.rename(columns={'Assay_Qualitative Measure':'Qualitative_Measure'},inplace=True)
tcell_full.rename(columns={'Assay_Quantitative measurement':'Quantitative_measurement'},inplace=True)
tcell_full.rename(columns={'Assay_Measurement Inequality':'Measurement_Inequality'},inplace=True)
tcell_full.rename(columns={'Epitope_Antigen IRI':'Epitope_Antigen_IRI'},inplace=True)
tcell_full.rename(columns={'MHC_MHC allele class':'MHC_class'},inplace=True)
tcell_full.rename(columns={'Assay Antigen_Antigen Epitope Relation':'Relation'},inplace=True)
tcell_full.rename(columns={'Assay_Assay Group':'Assay_Group'},inplace=True)

tcell_full.drop(["MHC_Allele Name","Epitope_Ending Position","Host_Name"],axis=1,inplace=True)

tcell_full.to_csv("Dataset_v2/tcell_full_CNRD.csv",index=False)



#入库，每类数据一个数据库，加快检索速度
tcell_full = pd.read_csv("Dataset_v2/tcell_full_CNRD.csv")
tcell_full=tcell_full[tcell_full["Relation"]=="Epitope"] #筛选属于Epitope，其余放入新抗原
print len(tcell_full)  #243270

tcell_full_text=[( "Species", "text"),( "Peptide", "text"),( "MHC_Type", "text")]
tcell_full_1 = [( "Species", 1),( "Peptide", 1),( "MHC_Type", 1),( "MHC_Allele", 1),( "Antigen_Name", 1),("type",1),("Assay_Group",1)]
#Tcell_activation库
tcell_full["type"] = "Tcell_activation库"


#处理另外部分入新抗原库
tcell_full = pd.read_csv("Dataset_v2/tcell_full_CNRD.csv")
#另外部分如何确定

