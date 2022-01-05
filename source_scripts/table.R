library(gtsummary)
dput(names(CGGA_phe))
CGGA_phe$Gender=ifelse(CGGA_phe$Gender=="Male","Male","Female")
CGGA_phe$status=ifelse(CGGA_phe$status=="1","Dead",
                    ifelse(CGGA_phe$status=="0","Alive",NA))

CGGA_phe$Radio_status=ifelse(CGGA_phe$Radio_status=="1","Yes1",
                          ifelse(CGGA_phe$Radio_status=="0","No1",NA))

CGGA_phe$Chemo_status=ifelse(CGGA_phe$Chemo_status=="1","Yes1",
                          ifelse(CGGA_phe$Chemo_status=="0","No1",NA))

c("Histology",
  "Grade", 
  "Recurrence", 
  "Subtype", 
  "IDH_status",  
  "Radio_status", 
  "Chemo_status", 
  "Age", 
  "Gender", 
  "survival", 
  "status", 
  "Group")

CGGA_phe2 <- CGGA_phe %>% select(
                           "survival", 
                           "status",
                           "Age", 
                           "Gender", 
                           "Grade", 
                           "Histology",
                           "Subtype",
                           "IDH_status",
                           "co_del_1p19q",
                           "Recurrence",
                           "group")
# 绘制
table1 <- tbl_summary(CGGA_phe2)

table1

table2 <-tbl_summary(
  CGGA_phe2,
  by = group, # 分组
  statistic = list(all_continuous() ~ "{mean} ({p25}, {p75})"),  
  missing = "no") %>%
  add_n() %>% # 添加非NA观测值个数
  add_p(test = list(all_continuous() ~ "t.test")) %>% # 添加P值
  add_overall() %>%
  modify_header(label = "**Variable**") %>% # 标签列header
  bold_labels()  #label 粗体

table2



