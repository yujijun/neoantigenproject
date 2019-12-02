from bs4 import BeautifulSoup
import urllib2
import sys
reload(sys)
sys.setdefaultencoding('utf8')
response = urllib2.urlopen("https://caped.icp.ucl.ac.be/Peptide/list")
content = response.read().decode("utf-8")
soup = BeautifulSoup(content,"html.parser")
# soup2 = BeautifulSoup(content,"lxml")

title={"mutationTable":"Mutation"}
title["tumorTable"]="Tumor-specific"
title["differentationTable"]="Differentiation"
title["overexpressedTable"]="Overexpressed"
out = pd.DataFrame()
n=0
for title_i,title_label in title.items():
#     print title_i
    content= soup.find_all("table",{"id":title_i})[0].find_all("tbody")[0]
    for i in content:
        if len(i) != 1:     
            t_n=0
            for t in i:
                if t_n ==  1:
                    Antigen_Name = t.get_text().strip()
                    out.loc[n,"Antigen_Name"] = Antigen_Name
                    try:
                        out.loc[n,"genecard_symbol"] = t.find_all("a",{"target":"_blank"})[0].attrs["href"].split("=")[1]
                    except:
                        pass
                if t_n ==  3:
                    out.loc[n,"Normal_Tumor"] = t.get_text().strip()
                if t_n ==  5:
                    hla = t.get_text().strip()
                    if hla!="":
                        try:
                            int(hla[-2])
                        except:
                            try:
                                int(hla[-1])
                            except:
                                pass
                            else:
                                hla=hla[:-1]+"*"+hla[-1:]
                        else:
                            hla=hla[:-2]+"*"+hla[-2:]
                        out.loc[n,"MHC_Type"] = "HLA-"+hla[0]
                        out.loc[n,"MHC_Allele"] = "HLA-"+hla
                if t_n ==  7:
                    out.loc[n,"HLA_Frequency"] = t.get_text()
                if t_n ==  9:
                    out.loc[n,"Peptide"] =  t.get_text()
                    out.loc[n,"Length_of_Peptide"] =len(t.get_text())
                    try:
                        red = t.find_all("span",{"style":"color:red"})[0].get_text()
                        if red !="":
                            out.loc[n,"Peptide_Type"] = "Mutation type"
                            Peptide = t.get_text("|").split("|")
                            if len(Peptide) ==1:
                                out.loc[n,"Mutated_Location"] = str(1)+"-"+str(len(Peptide[0]))
                            else:
                                if len(red)>1:
                                    out.loc[n,"Mutated_Location"] = str(len(Peptide[0])+1)+"-"+str(len(Peptide[0])+len(Peptide[1]))
                                else:
                                    out.loc[n,"Mutated_Location"] = str(len(Peptide[0])+1)
                    except:
                        out.loc[n,"Peptide_Type"] = "Wild type"                                                     
                if t_n ==  11:
                    loc = t.get_text()
                    if loc !="":
                        if "and" in loc:
                            #<td class="center_td">172-176 and 217-220</td>
                            #<td class="center_td">195-202 and 191 or 192e</td>
                            loc = loc.replace(" and ","-").replace(" or ","-").split("-")
                            out.loc[n,"Start_Position_in_Antigen"] =  loc[0]+"|"+loc[2]
                        #<td class="center_td">-i</td>
                        elif "-i" not in loc:
                            if "-" in loc:
                                if loc.startswith("("):
                                    loc = loc[1:]
                                loc= loc.split("(")[0].split("-")
                                out.loc[n,"Start_Position_in_Antigen"] =  int(loc[0])
                            else:
                                out.loc[n,"Start_Position_in_Antigen"] = loc
                if t_n ==  13:
                    out.loc[n,"Lymphocyte_stimulation"] = t.get_text()
                if t_n ==  17:
                    out.loc[n,"Reference"] = t.get_text().split(":")[-1].split(")")[0]                
                t_n += 1
            out.loc[n,"Species"]  = "Homo sapiens"
            out.loc[n,"type"]  = title_label
#             break
            n += 1
#     print out

#删除没有分型的数据
out.dropna(subset=["MHC_Type"],inplace=True)



out.to_csv("Dataset_v2/Neoantigen_CNRD.csv",index=False)

