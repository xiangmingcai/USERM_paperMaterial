
library(USERM)
library(reshape2)
library(ggplot2)
library(MASS)
library(ComplexHeatmap)
library(circlize)
library(grid)

saveDir = "E:/ResidualModel/Coef_vs_mSS"

# Step 1 querySig
Sig_info = querySig()

#Panel (19c) OMIP-099
prefix = "Panel_OMIP-099_adapted (19c)"
fluors_selected = c(Sig_info$id[c(139,140,143,144,
                                  149,150,151,152,
                                  153,154,155,156,
                                  158,159,160,161,
                                  162,201,196,
                                  138)])

#Panel (21c) OMIP-116
prefix = "Panel_OMIP-116_adapted (21c)"
fluors_selected = c(Sig_info$id[c(140,142,144,149,
                                  150,151,152,153,
                                  154,155,156,157,
                                  158,159,160,161,
                                  162,163,196,197,
                                  198,
                                  138)])



#Panel (23c) OMIP-117
prefix = "Panel_OMIP-117_adapted (23c)"
fluors_selected = c(Sig_info$id[c(140,141,149,150,
                                  151,152,153,154,
                                  155,156,157,158,
                                  159,160,161,162,
                                  174,189,201,196,
                                  198,199,200,
                                  138)])

#Panel 1 (40c) BV AURORA
prefix = "Panel_custom_1 (40c)"
fluors_selected = c(Sig_info$id[c(149:163,
                                  140:142,144:148,
                                  190,192:201,
                                  167,168,173,174,181,183,
                                  138)])

#Panel 2 (38c) NF AURORA
prefix = "Panel_custom_2 (38c)"
fluors_selected = c(Sig_info$id[c(145,146,147,148,
                                  193,194,195,
                                  139,140,141,144,
                                  199,201,149:162,
                                  168:170,172:174,177:178,181,183:184,
                                  138)])

#Panel 3 (45c) cFluor AURORA
prefix = "Panel_custom_3 (45c)"
fluors_selected = c(Sig_info$id[c(164,166:188,
                                  144,146:155,157:163,
                                  192,195,201,
                                  138)])

#Panel 4 (40c) BV XENITH
prefix = "Panel_custom_4 (40c)"
fluors_selected = c(Sig_info$id[c(71:85,
                                  66:69,70,
                                  127:137,
                                  115,117,120,124,126,
                                  95,106,108,105,
                                  64)])


#Panel 5 (42c) NF XENITH
prefix = "Panel_custom_5 (42c)"
fluors_selected = c(Sig_info$id[c(112,113,115,117,119:124,126,
                                  67,70,127,129:131,134:137,
                                  71:85,
                                  90,93,103,105,106,107,
                                  64)])


#Panel 6 (42c) cFluor XENITH
prefix = "Panel_custom_6 (42c)"
fluors_selected = c(Sig_info$id[c(86,89:96,99:103,105:110,
                                  66,70,130,134:136,
                                  114,115,118,121,124,126,
                                  71:75,77,78,79,80,85,
                                  64)])
# for (i in 1:length(fluors_selected)) {
#   fluors_selected[i] = strsplit(fluors_selected[i],"_")[[1]][3]
# }
# paste(fluors_selected,sep = ",",collapse = ",")

dev.off()
{
print(fluors_selected)
Sig_mtx  = getSigMtx(ids = fluors_selected)
dim(Sig_mtx)


UsermObj = CreateUserm(A = Sig_mtx)
#add ResObj into UsermObj
for (save_suf in colnames(Sig_mtx)) {
  ResObj = getRes(id = save_suf)
  UsermObj = AddRes2Userm(Res = ResObj, Userm = UsermObj)
}

Similarity_mtx = EstimateSimilarityMtx(A = UsermObj$A)
pdf(file = file.path(saveDir,paste0(prefix,"_Similarity_mtx.pdf")),width = 14,height = 14)
p = Vis_Mtx(mat = Similarity_mtx,mincolor = "white",midcolor = "white", maxcolor = "#95ABDB",
        max = 1,mid = 0.8,min = 0,legend_name = "Cosine",
        title = "Cosine similarity matrix")
print(p)
dev.off()


ssm = EstimateSSM(SSM_fluor = UsermObj$fluors[1:(length(UsermObj$fluors)-1)],
                  Userm = UsermObj,mSS = F,
                  A = UsermObj$A,quiet = T)
# pdf(file = file.path(saveDir,paste0(prefix,"_SSM.pdf")),width = 6,height = 6)
# p = Vis_Mtx(mat = ssm,mincolor = "white",midcolor = "#E3E7F2", maxcolor = "#95ABDB",
#             max = 10,mid = 5,min = 0,legend_name = "SS",text_reg = "%.1f",fontsize = 6,
#             title = "Spillover Spreading Matrix")
# print(p)
# dev.off()

pdf(file = file.path(saveDir,paste0(prefix,"_SSM.pdf")),width = 10,height = 10)
p = Vis_Mtx(mat = ssm,mincolor = "white",midcolor = "#E3E7F2", maxcolor = "#95ABDB",
            max = 10,mid = 5,min = 0,legend_name = "SS",text_reg = "%.1f",
            title = "Spillover Spreading Matrix")
print(p)
dev.off()

Coef_mtx = EstimateCoefMtx(Userm = UsermObj)
Coef_mtx = Coef_mtx[rownames(ssm), colnames(ssm)]
# pdf(file = file.path(saveDir,paste0(prefix,"_Coef_mtx.pdf")),width = 6,height = 6)
# p = Vis_Mtx(mat = Coef_mtx,mincolor = "white",midcolor = "#E3E7F2", maxcolor = "#95ABDB",
#             max = 10,mid = 5,min = 0,legend_name = "Coef",text_reg = "%.1f",fontsize = 6,
#             title = "Coefficient of residual model matrix (Coef Matrix) (row spread into column)")
# print(p)
# dev.off()
pdf(file = file.path(saveDir,paste0(prefix,"_Coef_mtx.pdf")),width = 10,height = 10)
p = Vis_Mtx(mat = Coef_mtx,mincolor = "white",midcolor = "#E3E7F2", maxcolor = "#95ABDB",
            max = 10,mid = 5,min = 0,legend_name = "Coef",text_reg = "%.1f",
            title = "Coefficient of residual model matrix (Coef Matrix) (row spread into column)")
print(p)
dev.off()

ssm_long <- melt(ssm)
Coef_long <- melt(Coef_mtx)

all(ssm_long$Var1 == Coef_long$Var1)
all(ssm_long$Var2 == Coef_long$Var2)
metric_long = Coef_long
colnames(metric_long) = c("row","col","coef")
metric_long$mSS = ssm_long$value
metric_long = metric_long[metric_long$row != metric_long$col,]

cor_test <- cor.test(metric_long$coef, metric_long$mSS)
r_val <- round(cor_test$estimate, 3)
# p_val <- signif(cor_test$p.value)
p_val <- sprintf("%.2e", cor_test$p.value)

pdf(file = file.path(saveDir,paste0(prefix,"_Dotplot_Correlation.pdf")),width = 6,height = 6)
p = ggplot(metric_long, aes(x = coef, y = mSS)) +
  geom_point(size = 3, alpha = 0.8, color = "#2C7BB6") +
  geom_smooth(method = "lm", se = FALSE, color = "#D7191C", linewidth = 1) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Coefficient",
    y = "Spreading Error (SS)",
    title = "Relationship Between Coefficient and SS",
    subtitle =paste0(prefix, paste0("\n","R = ", round(r_val,digits = 3),  "; p = ", p_val)),
    caption = "Each point represents one marker pair"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, color = "gray40"),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
print(p)
dev.off()
}
colnames(ssm)
ss_output = check_ss(f_pos = "SCC3_Cell_cFluorBYG750_CD4",
                     f_neg = "SCC3_Cell_BV510_CD4",
                     SSM_fluor = UsermObj$fluors[1:(length(UsermObj$fluors)-1)],
                     A = UsermObj$A,
                     Userm = UsermObj)

vis_ss_scatter(ss_output)

checkSig_linePlot(id = "SCC3_Cell_cFluorBYG750_CD4")
checkSig_linePlot(id = "SCC3_Cell_BV510_CD4")


fluor = "SCC3_Cell_BV750_CD4"
ssm_obj = readRDS(system.file("ssm", paste0("SSMObj_",fluor,".rds"), package = "USERM"))
export_data = ssm_obj$df_pos
export_A = UsermObj$A[,c("SCC3_Cell_BV750_CD4","SCC3_Cell_BV711_CD4")]

library(flowCore)
ff <- flowFrame(as.matrix(export_data)) #
pData(parameters(ff))$desc <- colnames(export_data)
write.FCS(ff, file.path(saveDir,"BV750.fcs"))

Sig_df = as.data.frame(t(export_A))
Sig_df$Primary = NA
Sig_df$Secondary = NA
Sig_df <- Sig_df[c("Primary", "Secondary",colnames(Sig_df)[1:64])]
for (i in 1:nrow(Sig_df)) {
  Sig_df[i,"Primary"] = strsplit(rownames(Sig_df)[i],split = "_")[[1]][3]
  Sig_df[i,"Secondary"] = strsplit(rownames(Sig_df)[i],split = "_")[[1]][4]
}
write.table(Sig_df,file = file.path(saveDir,"SigMtx_BV750.csv"),row.names = F,quote = F,sep = ",")

corrected_SigMtx = read.csv2(file = file.path(saveDir,"corrected_SigMtx_BV750.csv"),sep = ",")
corrected_SigMtx = t(corrected_SigMtx[,3:ncol(corrected_SigMtx)])

UsermObj$A[,c("SCC3_Cell_BV750_CD4")] = as.numeric(corrected_SigMtx[,1])
UsermObj$A[,c("SCC3_Cell_BV711_CD4")] = as.numeric(corrected_SigMtx[,2])

colnames(A_pinv) = rownames(A)
df_pos = ssm_obj$df_pos
colnames(B_pos)
ggplot(ss_output$B_pos, aes(x = SCC3_Cell_BV750_CD4, y = SCC3_Cell_BV711_CD4)) +
  geom_point(size = 3, alpha = 0.8, color = "#2C7BB6")


B_pos
