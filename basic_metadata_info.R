# Daniel Castaneda Mogollon, PhD
# April 23rd, 2026


path = "/Users/danielcm/Desktop/MIGUT/"
setwd(path)
df = read.csv("MIGUT_metadata_in_tab_format.txt", header = TRUE, sep = "\t", row.names = 1)

parkinsons_df = df[df$Group == "PD", ]
healthy_df = df[df$Group == "HC", ]


dopamine_agonists = parkinsons_df[parkinsons_df$name..mgx.tab.!=0,]
dopamine_agonists = dopamine_agonists[!is.na(dopamine_agonists$name..mgx.tab.),]
dopamine_agonists = dopamine_agonists[dopamine_agonists$name..mgx.tab.!="missing",]

dopamine_agonists$name..mgx.tab.=="Domperidone"

1/38
105/124
19/124
15/125
14/113
100/112

t.test()