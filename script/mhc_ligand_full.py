#######################解析mhc_ligand_full件，获得入库前文件（ligand、mhcI和II）
#制作含有所需行的小文件
mhc_ligand_full = open("Dataset_v2/mhc_ligand_full.csv")
head1 = pd.DataFrame(mhc_ligand_full.readline().split(","))
head2 =  pd.DataFrame(mhc_ligand_full.readline().split(","))
head = head1+"_"+head2
head = head[0].tolist()
mhc_ligand_full = pd.read_csv("Dataset_v2/mhc_ligand_full.csv",skiprows= 2)
mhc_ligand_full.columns=head
print(len(mhc_ligand_full))
select_col=["MHC_Allele Name","Epitope_Description","Epitope_Starting Position","Epitope_Ending Position"]
select_col.extend(["Epitope_Organism Name","Host_Name","Epitope_Antigen Name","Epitope_Antigen IRI","Reference_PubMed ID"])
select_col.extend(["Assay_Units","Assay_Qualitative Measure","Assay_Measurement Inequality","Assay_Quantitative measurement"])
select_col.extend(["MHC_MHC allele class"])
mhc_ligand_full=mhc_ligand_full[select_col]
mhc_ligand_full.to_csv("mhc_ligand_full_select_col.csv",index=False)


mhc_ligand_full = pd.read_csv("Dataset_v2/mhc_ligand_full_select_col.csv")
#确定如何保留物种
mhc=mhc_ligand_full["MHC_Allele Name"].drop_duplicates().tolist()
print len(mhc)

mhc_1 = mhc[:]
for i in mhc:
    if i.startswith("HLA") or i.startswith("H2") or i.startswith("RT1"):
        mhc_1.remove(i)
print mhc_1

print len(mhc_ligand_full)
#通过MHC的拼写选取三个物种
mhc_ligand_full=pd.concat([mhc_ligand_full[mhc_ligand_full["MHC_Allele Name"].str.startswith("HLA")],
                           mhc_ligand_full[mhc_ligand_full["MHC_Allele Name"].str.startswith("H2")] ,
                           mhc_ligand_full[mhc_ligand_full["MHC_Allele Name"].str.startswith("RT1")]],)

#页面展示prtein id
mhc_ligand_full["protein_id"]=mhc_ligand_full["Epitope_Antigen IRI"].str.split("/",expand=True).iloc[:,-1]
#物种根据MHC的拼写规则定义
mhc_ligand_full["Species"]="Homo sapiens"
mhc_ligand_full.loc[mhc_ligand_full["MHC_Allele Name"].str.startswith("H2"),"Species"] ="Mus musculus"
mhc_ligand_full.loc[mhc_ligand_full["MHC_Allele Name"].str.startswith("RT1"),"Species"] ="Rattus norvegicus"

#定义MHC类型及修正亚类的书写格式
mhc_type = pd.DataFrame(mhc_ligand_full["MHC_Allele Name"].drop_duplicates())
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

mhc_ligand_full=pd.merge(mhc_ligand_full,mhc_type,on="MHC_Allele Name",how="left")


mhc_ligand_full["Length_of_Peptide"] = mhc_ligand_full["Epitope_Ending Position"]-mhc_ligand_full["Epitope_Starting Position"]

#通过蛋白名称获得基因symbol
#该文件通过url去重后共有119347中链接，三类分别为ncbi、uniprot和iedb的ontology，前两种通过网页均能查到基因名称，后一种
#该问题暂放，看是否有其他更佳方式获得基因symbol
url = pd.DataFrame(mhc_ligand_full["Epitope_Antigen IRI"].dropna().drop_duplicates())
url.index=range(len(url))  

#重新规范定义列名,去除不必要的列
mhc_ligand_full.rename(columns={'Epitope_Description':'Peptide'},inplace=True)
mhc_ligand_full.rename(columns={'Reference_PubMed ID':'Reference'},inplace=True)
mhc_ligand_full.rename(columns={'Epitope_Starting Position':'Start_Position_in_Antigen'},inplace=True)
mhc_ligand_full.rename(columns={'Epitope_Organism Name':'Epitope_Organism'},inplace=True)
mhc_ligand_full.rename(columns={'Epitope_Antigen Name':'Antigen_Name'},inplace=True)
mhc_ligand_full.rename(columns={'Assay_Qualitative Measure':'Qualitative_Measure'},inplace=True)
mhc_ligand_full.rename(columns={'Assay_Quantitative measurement':'Quantitative_measurement'},inplace=True)
mhc_ligand_full.rename(columns={'Assay_Measurement Inequality':'Measurement_Inequality'},inplace=True)
mhc_ligand_full.rename(columns={'Epitope_Antigen IRI':'Epitope_Antigen_IRI'},inplace=True)
mhc_ligand_full.rename(columns={'MHC_MHC allele class':'MHC_class'},inplace=True)
mhc_ligand_full.drop(["MHC_Allele Name","Epitope_Ending Position","Host_Name"],axis=1,inplace=True)

mhc_ligand_full.to_csv("Dataset_v2/mhc_ligand_full_CNRD.csv",index=False)



#入库，每类数据一个数据库，加快检索速度
mhc_ligand_full = pd.read_csv("Dataset_v2/mhc_ligand_full_CNRD.csv")
mhc_ligand_full_text=[( "Species", "text"),( "Peptide", "text"),( "MHC_Type", "text")]
mhc_ligand_full_1 = [( "Species", 1),( "Peptide", 1),( "MHC_Type", 1),( "MHC_Allele", 1),( "Antigen_Name", 1),("type",1)]
# print len(mhc_ligand_full)
#Ligands库
mhc_ligand_full["type"] = "Ligands"


#MHCI Binding库
mhcI=mhc_ligand_full[mhc_ligand_full["MHC_class"]=="I"]
# print len(mhcI[pd.isnull(mhcI["Qualitative_Measure"])]) #定性都有结果
mhcI["type"] = "MHCI_Binding"


# MHCII Binding库
mhcII=mhc_ligand_full[mhc_ligand_full["MHC_class"]=="II"]
# print len(mhcII[pd.isnull(mhcII["Qualitative_Measure"])]) #定性都有结果
mhcII["type"] = "MHCII_Binding"
