
library(ggplot2)
library(ggbreak)

dataDir = "E:/ResidualModel/validation/data"
saveDir = "E:/ResidualModel/validation/output"

# internal validation Xenith BUV805_CD16 ####

SCC_Cell_BUV805_CD16_12 = read.csv2(file = file.path(dataDir,"Internal validation","Xenith","SCC_Cell_BUV805_CD16","SCC_Cell_BUV805_CD16 from Xenith (11 sig unmixing).csv"),sep = ",")
SCC_Cell_BUV805_CD16_12$Parameter = "SCC Cell"
SCC_Cell_BUV805_CD16_12$PanelSize = "11"
SCC_Cell_BUV805_CD16_22 = read.csv2(file = file.path(dataDir,"Internal validation","Xenith","SCC_Cell_BUV805_CD16","SCC_Cell_BUV805_CD16 from Xenith (21 sig unmixing).csv"),sep = ",")
SCC_Cell_BUV805_CD16_22$Parameter = "SCC Cell"
SCC_Cell_BUV805_CD16_22$PanelSize = "21"
SCC_Cell_BUV805_CD16_32 = read.csv2(file = file.path(dataDir,"Internal validation","Xenith","SCC_Cell_BUV805_CD16","SCC_Cell_BUV805_CD16 from Xenith (31 sig unmixing).csv"),sep = ",")
SCC_Cell_BUV805_CD16_32$Parameter = "SCC Cell"
SCC_Cell_BUV805_CD16_32$PanelSize = "31"

SCC_Bead_BUV805_CD16_12 = read.csv2(file = file.path(dataDir,"Internal validation","Xenith","SCC Bead BUV805_CD16","SCC_Bead_BUV805_CD16 from Xenith (11 sig unmixing).csv"),sep = ",")
SCC_Bead_BUV805_CD16_12$Parameter = "SCC Bead"
SCC_Bead_BUV805_CD16_12$PanelSize = "11"
SCC_Bead_BUV805_CD16_22 = read.csv2(file = file.path(dataDir,"Internal validation","Xenith","SCC Bead BUV805_CD16","SCC_Bead_BUV805_CD16 from Xenith (21 sig unmixing).csv"),sep = ",")
SCC_Bead_BUV805_CD16_22$Parameter = "SCC Bead"
SCC_Bead_BUV805_CD16_22$PanelSize = "21"
SCC_Bead_BUV805_CD16_32 = read.csv2(file = file.path(dataDir,"Internal validation","Xenith","SCC Bead BUV805_CD16","SCC_Bead_BUV805_CD16 from Xenith (31 sig unmixing).csv"),sep = ",")
SCC_Bead_BUV805_CD16_32$Parameter = "SCC Bead"
SCC_Bead_BUV805_CD16_32$PanelSize = "31"

data = rbind(SCC_Cell_BUV805_CD16_12,SCC_Cell_BUV805_CD16_22)
data = rbind(data,SCC_Cell_BUV805_CD16_32)
data = rbind(data,SCC_Bead_BUV805_CD16_12)
data = rbind(data,SCC_Bead_BUV805_CD16_22)
data = rbind(data,SCC_Bead_BUV805_CD16_32)
data$Coverage.Rate = as.numeric(data$Coverage.Rate)
data$Parameter = factor(x = data$Parameter,levels = c("SCC Cell","SCC Bead"))
data$PanelSize = factor(x = data$PanelSize,levels = c("11","21","31"))


p = ggplot(data, aes(x = PanelSize, y = Coverage.Rate, fill = Parameter)) +
  geom_boxplot(
    width = 0.65,
    position = position_dodge(width = 0.75),
    outlier.shape = NA,
    color = "black",
    linewidth = 0.5
  ) +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.75
    ),
    size = 1.8,
    alpha = 0.5,
    stroke = 0
  ) +
  geom_hline(
    yintercept = 0.95,
    linetype = "dashed",
    color = "#666666",
    linewidth = 0.6
  ) +
  scale_fill_manual(
    name = "Parameter source",
    values = c("#4C72B0", "#DD8452", "#55A868", "#C44E52")
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(seq(0, 1, 0.05))
  ) +
  labs(
    x = "Panel Size",
    y = "Coverage Rate",
    title = "Internal validation",
    subtitle = "Spread prediction of unmixed SCC cell BUV805-CD16",
    caption = "Each point represents unmixed spread at one axis"
  ) +
  guides(fill = guide_legend(nrow = 1)) +
  scale_y_break(
    c(0.06, 0.89),
    expand = c(0, 0)
  ) +
  theme_classic(base_size = 14) +
  theme(
    text = element_text(family = "Helvetica"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    axis.text = element_text(color = "black"),
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 14, face = "bold"),
    legend.position = "top",
    legend.key.size = unit(0.8, "cm"),
    legend.background = element_blank(),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12),
    plot.caption = element_text(size = 9, color = "gray40"),
    panel.spacing = unit(1, "lines"),
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line.y.right = element_blank()
  )
pdf(file = file.path(saveDir,paste0("Internalvalidation_Xenith_BUV805_CD16.pdf")),width = 6,height = 6,onefile =F)
print(p)
dev.off()


# internal validation Xenith Cells ####

files = dir(file.path(dataDir,"Internal validation","Xenith","SCC_Cell"))

# file = files[1]
for (file in files) {
  if(file == files[1]){
    data_tmp = read.csv2(file = file.path(dataDir,"Internal validation","Xenith","SCC_Cell",file,"coverage_summary.csv"),sep = ",")
    data_tmp = data_tmp[data_tmp$Dim != file,]
    data_tmp = data_tmp[data_tmp$Dim != "AF",]
    for (i in 1:nrow(data_tmp)) {
      data_tmp[i,"Dim"] = paste0(strsplit(data_tmp[i,"Dim"],"_")[[1]][3],"_",strsplit(data_tmp[i,"Dim"],"_")[[1]][4])
    }
    data_tmp$Dim
    data_tmp$SCC = paste0(strsplit(file,"_")[[1]][3],"_",strsplit(file,"_")[[1]][4])

    data = data_tmp
  }else{
    data_tmp = read.csv2(file = file.path(dataDir,"Internal validation","Xenith","SCC_Cell",file,"coverage_summary.csv"),sep = ",")
    data_tmp = data_tmp[data_tmp$Dim != file,]
    data_tmp = data_tmp[data_tmp$Dim != "AF",]
    for (i in 1:nrow(data_tmp)) {
      data_tmp[i,"Dim"] = paste0(strsplit(data_tmp[i,"Dim"],"_")[[1]][3],"_",strsplit(data_tmp[i,"Dim"],"_")[[1]][4])
    }
    data_tmp$Dim
    data_tmp$SCC = paste0(strsplit(file,"_")[[1]][3],"_",strsplit(file,"_")[[1]][4])

    data = rbind(data,data_tmp)
  }
}

data$Coverage.Rate = as.numeric(data$Coverage.Rate)

p = ggplot(data, aes(x = SCC, y = Coverage.Rate)) +
  geom_boxplot(
    width = 0.65,
    position = position_dodge(width = 0.75),
    outlier.shape = NA,
    color = "black",
    fill = "#4C72B0",
    linewidth = 0.5
  ) +
  geom_jitter(
    width = 0.15,
    size = 1.8,
    alpha = 0.5
  ) +
  geom_hline(
    yintercept = 0.95,
    linetype = "dashed",
    color = "#666666",
    linewidth = 0.6
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(seq(0, 1, 0.05))
  ) +
  labs(
    x = "SSC",
    y = "Coverage Rate",
    title = "Internal validation",
    subtitle = "Spread prediction of unmixed SCC cell samples",
    caption = "Each point represents unmixed spread at one axis"
  ) +
  scale_y_break(
    c(0.06, 0.69),
    expand = c(0, 0)
  ) +
  theme_classic(base_size = 14) +
  theme(
    text = element_text(family = "Helvetica"),
    axis.line.y = element_line(color = "black", linewidth = 0.6),  # ← 必须加，确保闪电符号可见
    axis.line.x = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    axis.text = element_text(color = "black"),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 14, face = "bold"),
    legend.position = "top",
    legend.key.size = unit(0.8, "cm"),
    legend.background = element_blank(),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12),
    plot.caption = element_text(size = 9, color = "gray40"),
    panel.spacing = unit(1, "lines"),
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line.y.right = element_blank()
  )
pdf(file = file.path(saveDir,paste0("Internalvalidation_Xenith_SCC_Cell.pdf")),width = 12,height = 6,onefile =F)
print(p)
dev.off()




# internal validation Aurora5L Cells ####

files = dir(file.path(dataDir,"Internal validation","Aurora5L","SCC_Cell"))

# file = files[1]
for (file in files) {
  if(file == files[1]){
    data_tmp = read.csv2(file = file.path(dataDir,"Internal validation","Aurora5L","SCC_Cell",file,"coverage_summary.csv"),sep = ",")
    data_tmp = data_tmp[data_tmp$Dim != file,]
    data_tmp = data_tmp[data_tmp$Dim != "AF",]
    for (i in 1:nrow(data_tmp)) {
      data_tmp[i,"Dim"] = paste0(strsplit(data_tmp[i,"Dim"],"_")[[1]][3],"_",strsplit(data_tmp[i,"Dim"],"_")[[1]][4])
    }
    data_tmp$Dim
    data_tmp$SCC = paste0(strsplit(file,"_")[[1]][3],"_",strsplit(file,"_")[[1]][4])

    data = data_tmp
  }else{
    data_tmp = read.csv2(file = file.path(dataDir,"Internal validation","Aurora5L","SCC_Cell",file,"coverage_summary.csv"),sep = ",")
    data_tmp = data_tmp[data_tmp$Dim != file,]
    data_tmp = data_tmp[data_tmp$Dim != "AF",]
    for (i in 1:nrow(data_tmp)) {
      data_tmp[i,"Dim"] = paste0(strsplit(data_tmp[i,"Dim"],"_")[[1]][3],"_",strsplit(data_tmp[i,"Dim"],"_")[[1]][4])
    }
    data_tmp$Dim
    data_tmp$SCC = paste0(strsplit(file,"_")[[1]][3],"_",strsplit(file,"_")[[1]][4])

    data = rbind(data,data_tmp)
  }
}

data$Coverage.Rate = as.numeric(data$Coverage.Rate)


p = ggplot(data, aes(x = SCC, y = Coverage.Rate)) +
  geom_boxplot(
    width = 0.65,
    position = position_dodge(width = 0.75),
    outlier.shape = NA,
    color = "black",
    fill = "#4C72B0",
    linewidth = 0.5
  ) +
  geom_jitter(
    width = 0.15,
    size = 1.8,
    alpha = 0.5
  ) +
  geom_hline(
    yintercept = 0.95,
    linetype = "dashed",
    color = "#666666",
    linewidth = 0.6
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(seq(0, 1, 0.05))
  ) +
  labs(
    x = "SSC",
    y = "Coverage Rate",
    title = "Internal validation",
    subtitle = "Spread prediction of unmixed SCC cell samples",
    caption = "Each point represents unmixed spread at one axis"
  ) +
  scale_y_break(
    c(0.01, 0.84),
    expand = c(0, 0)
  ) +

  theme_classic(base_size = 14) +
  theme(
    text = element_text(family = "Helvetica"),
    axis.line.y = element_line(color = "black", linewidth = 0.6),  # ← 必须加，确保闪电符号可见
    axis.line.x = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    axis.text = element_text(color = "black"),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 14, face = "bold"),
    legend.position = "top",
    legend.key.size = unit(0.8, "cm"),
    legend.background = element_blank(),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12),
    plot.caption = element_text(size = 9, color = "gray40"),
    panel.spacing = unit(1, "lines"),
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line.y.right = element_blank()
  )
pdf(file = file.path(saveDir,paste0("Internalvalidation_Aurora5L_SCC_Cell.pdf")),width = 18,height = 6,onefile =F)
print(p)
dev.off()


# external validation Xenith Cells ####

files = dir(file.path(dataDir,"External validation","Xenith","SCC_Cell"))

# file = files[1]
for (file in files) {
  if(file == files[1]){
    data_tmp_beads = read.csv2(file = file.path(dataDir,"External validation","Xenith","SCC_Cell",file,"coverage_beads_corrected_summary.csv"),sep = ",")
    data_tmp_beads = data_tmp_beads[data_tmp_beads$Dim != file,]
    data_tmp_beads = data_tmp_beads[data_tmp_beads$Dim != "AF",]
    for (i in 1:nrow(data_tmp_beads)) {
      data_tmp_beads[i,"Dim"] = paste0(strsplit(data_tmp_beads[i,"Dim"],"_")[[1]][3],"_",strsplit(data_tmp_beads[i,"Dim"],"_")[[1]][4])
    }

    data_tmp_beads$SCC = paste0(strsplit(file,"_")[[1]][3],"_","CD4")
    data_tmp_beads$Parameter = "SCC Bead"

    data_tmp_cells = read.csv2(file = file.path(dataDir,"External validation","Xenith","SCC_Cell",file,"coverage_cells_corrected_summary.csv"),sep = ",")
    data_tmp_cells = data_tmp_cells[data_tmp_cells$Dim != file,]
    data_tmp_cells = data_tmp_cells[data_tmp_cells$Dim != "AF",]
    for (i in 1:nrow(data_tmp_cells)) {
      data_tmp_cells[i,"Dim"] = paste0(strsplit(data_tmp_cells[i,"Dim"],"_")[[1]][3],"_",strsplit(data_tmp_cells[i,"Dim"],"_")[[1]][4])
    }

    data_tmp_cells$SCC = paste0(strsplit(file,"_")[[1]][3],"_","CD4")
    data_tmp_cells$Parameter = "SCC Cell"

    data_tmp = rbind(data_tmp_beads,data_tmp_cells)

    data = data_tmp
  }else{
    data_tmp_beads = read.csv2(file = file.path(dataDir,"External validation","Xenith","SCC_Cell",file,"coverage_beads_corrected_summary.csv"),sep = ",")
    data_tmp_beads = data_tmp_beads[data_tmp_beads$Dim != file,]
    data_tmp_beads = data_tmp_beads[data_tmp_beads$Dim != "AF",]
    for (i in 1:nrow(data_tmp_beads)) {
      data_tmp_beads[i,"Dim"] = paste0(strsplit(data_tmp_beads[i,"Dim"],"_")[[1]][3],"_",strsplit(data_tmp_beads[i,"Dim"],"_")[[1]][4])
    }

    data_tmp_beads$SCC = paste0(strsplit(file,"_")[[1]][3],"_","CD4")
    data_tmp_beads$Parameter = "SCC Bead"

    data_tmp_cells = read.csv2(file = file.path(dataDir,"External validation","Xenith","SCC_Cell",file,"coverage_cells_corrected_summary.csv"),sep = ",")
    data_tmp_cells = data_tmp_cells[data_tmp_cells$Dim != file,]
    data_tmp_cells = data_tmp_cells[data_tmp_cells$Dim != "AF",]
    for (i in 1:nrow(data_tmp_cells)) {
      data_tmp_cells[i,"Dim"] = paste0(strsplit(data_tmp_cells[i,"Dim"],"_")[[1]][3],"_",strsplit(data_tmp_cells[i,"Dim"],"_")[[1]][4])
    }

    data_tmp_cells$SCC = paste0(strsplit(file,"_")[[1]][3],"_","CD4")
    data_tmp_cells$Parameter = "SCC Cell"

    data_tmp = rbind(data_tmp_beads,data_tmp_cells)

    data = rbind(data,data_tmp)
  }
}

data$Coverage.Rate = as.numeric(data$Coverage.Rate)
data$Parameter = factor(x = data$Parameter,levels = c("SCC Cell","SCC Bead"))


p = ggplot(data, aes(x = SCC, y = Coverage.Rate, fill = Parameter)) +
  geom_boxplot(
    width = 0.65,
    position = position_dodge(width = 0.75),
    outlier.shape = NA,
    color = "black",
    linewidth = 0.5
  ) +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.75
    ),
    size = 1.8,
    alpha = 0.5,
    stroke = 0
  ) +
  geom_hline(
    yintercept = 0.95,
    linetype = "dashed",
    color = "#666666",
    linewidth = 0.6
  ) +
  scale_fill_manual(
    name = "Parameter source",
    values = c("#4C72B0", "#DD8452", "#55A868", "#C44E52")
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(seq(0, 1, 0.05))
  ) +
  labs(
    x = "SCC",
    y = "Coverage Rate",
    title = "External validation",
    subtitle = "Spread prediction of unmixed SCC cell",
    caption = "Each point represents unmixed spread at one axis"
  ) +
  guides(fill = guide_legend(nrow = 1)) +
  scale_y_break(
    c(0.06, 0.59),
    expand = c(0, 0)
  ) +
  theme_classic(base_size = 14) +
  theme(
    text = element_text(family = "Helvetica"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    axis.text = element_text(color = "black"),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 14, face = "bold"),
    legend.position = "top",
    legend.key.size = unit(0.8, "cm"),
    legend.background = element_blank(),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12),
    plot.caption = element_text(size = 9, color = "gray40"),
    panel.spacing = unit(1, "lines"),
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line.y.right = element_blank()
  )
pdf(file = file.path(saveDir,paste0("Externalvalidation_Xenith_Cells.pdf")),width = 18,height = 6,onefile =F)
print(p)
dev.off()
