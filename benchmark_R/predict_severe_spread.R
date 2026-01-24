
library(USERM)
library(reshape2)
library(ggplot2)
library(MASS)


Sig_info = querySig()

saveDir = "E:/ResidualModel/predict_severe_spread"

#original 31 fluors panel
prefix = "Original_Panel(31color)"
fluors_selected = c(Sig_info$id[c(32:38,41:42,44,46:48,50:54,58,59,60,61,62,8,9,12,14,18,24,25,26,63)])
# checkSig_linePlot(id = "SCC_Bead_TCRVa24Ja18_BV785")
# checkSig_linePlot(id = "SCC_Cell_TCRVa24Ja18_BV785")

#original 31 fluors panel + one overlap fluor
prefix = "AddOne_Panel(31color+cFluorV610)"
fluors_selected = c(Sig_info$id[c(32:38,41:42,44,46:48,50:54,58,59,60,61,62,8,9,12,14,18,24,25,26,
                                  108,##good: 110,108,107,106,124,99,
                                  63)])

#original 31 fluors panel + one overlap fluor
prefix = "AddOne_Panel(31color+cFluorV505)"
fluors_selected = c(Sig_info$id[c(32:38,41:42,44,46:48,50:54,58,59,60,61,62,8,9,12,14,18,24,25,26,
                                  106,##good: 110,108,107,106,124,99,
                                  63)])

#original 31 fluors panel + one overlap fluor
prefix = "AddOne_Panel(31color+NFY660)"
fluors_selected = c(Sig_info$id[c(32:38,41:42,44,46:48,50:54,58,59,60,61,62,8,9,12,14,18,24,25,26,
                                  124,##good: 110,108,107,106,124,99,
                                  63)])

#original 31 fluors panel + one overlap fluor
prefix = "AddOne_Panel(31color+ AF700)"
fluors_selected = c(Sig_info$id[c(32:38,41:42,44,46:48,50:54,58,59,60,61,62,8,9,12,14,18,24,25,26,
                                  65,#65:85
                                  63)])

{
print(fluors_selected)
Sig_mtx  = getSigMtx(ids = fluors_selected)
UsermObj = CreateUserm(A = Sig_mtx)
#add ResObj into UsermObj
for (save_suf in colnames(Sig_mtx)) {
  ResObj = getRes(id = save_suf)
  UsermObj = AddRes2Userm(Res = ResObj, Userm = UsermObj)
}

Similarity_mtx = EstimateSimilarityMtx(A = UsermObj$A)
pdf(file = file.path(saveDir,paste0(prefix,"_Similarity_mtx.pdf")),width = 14,height = 14)
p = Vis_Mtx(mat = Similarity_mtx,mincolor = "white",midcolor = "white", maxcolor = "#95ABDB",
            max = 1,mid = 0.8,min = 0,legend_name = "Cosine",text_reg = "%.2f",
            title = "Cosine similarity matrix")
print(p)
dev.off()

Hotspot_mtx = EstimateHotspotMtx(A = UsermObj$A)
pdf(file = file.path(saveDir,paste0(prefix,"_Hotspot_mtx.pdf")),width = 10,height = 10)
p = Vis_Mtx(mat = Hotspot_mtx,mincolor = "white",midcolor = "white", maxcolor = "#95ABDB",
        max = 2,mid = 1,min = 0,legend_name = "Hotspot",text_reg = "%.2f",
        title = "Hotspot matrix")
print(p)
dev.off()

ssm = EstimateSSM(SSM_fluor = UsermObj$fluors[1:(length(UsermObj$fluors)-1)],
                  Userm = UsermObj,
                  A = UsermObj$A,quiet = T)
pdf(file = file.path(saveDir,paste0(prefix,"_SSM.pdf")),width = 8,height = 8)
p = Vis_Mtx(mat = ssm,mincolor = "white",midcolor = "white", maxcolor = "#95ABDB",
            max = 2,mid = 0.5,min = 0,legend_name = "SS",text_reg = "%.1f",
            title = "Spillover Spreading Matrix")
print(p)
dev.off()

Coef_mtx = EstimateCoefMtx(Userm = UsermObj,A = UsermObj$A)
Coef_mtx = Coef_mtx[rownames(ssm), colnames(ssm)]
pdf(file = file.path(saveDir,paste0(prefix,"_Coef_mtx.pdf")),width = 8,height = 8)
p = Vis_Mtx(mat = Coef_mtx,mincolor = "white",midcolor = "white", maxcolor = "#95ABDB",
            max = 2,mid = 0.5,min = 0,legend_name = "Coef",text_reg = "%.1f",
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

pdf(file = file.path(saveDir,paste0(prefix,"_Dotplot_Correlation.pdf")),width = 6,height = 6)
p = ggplot(metric_long, aes(x = coef, y = mSS)) +
  geom_point(size = 3, alpha = 0.8, color = "#2C7BB6") +
  geom_smooth(method = "lm", se = FALSE, color = "#D7191C", linewidth = 1) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Coefficient",
    y = "Modified Spreading Error (mSS)",
    title = "Relationship Between Coefficient and mSS",
    subtitle = prefix,
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
Sig_df = as.data.frame(t(UsermObj$A * 100))
# rownames(Sig_df)[nrow(Sig_df)-1] = "SCC_Cell_add_new"
Sig_df$Primary = NA
Sig_df$Secondary = NA
Sig_df <- Sig_df[c("Primary", "Secondary",rownames(Sig_mtx))]
for (i in 1:nrow(Sig_df)) {
  Sig_df[i,"Primary"] = strsplit(rownames(Sig_df)[i],split = "_")[[1]][3]
  Sig_df[i,"Secondary"] = strsplit(rownames(Sig_df)[i],split = "_")[[1]][4]
}
Sig_df$Primary[Sig_df$Primary == "TCRrd"] = "TCRgd"

write.table(Sig_df,file = file.path(saveDir,paste0(prefix,"_SigMtx.csv")),row.names = F,quote = F,sep = ",")


ss_output = check_ss(f_pos = "SCC_Cell_CD45_AF532",
                     f_neg = "SCC_Cell_CRTH2_PECy5",
                     SSM_fluor = UsermObj$fluors[1:(length(UsermObj$fluors)-1)],
                     A = UsermObj$A,
                     Userm = UsermObj)

vis_ss_scatter(ss_output)
