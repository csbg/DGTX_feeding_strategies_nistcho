# Figure 2 
ggarrange(
    IVCD_bar_STD, qp_bar_STD, growth_plot, qp_timecourse, qp_gluc, p_lactate_qp,
    labels = c("A", "B", "C", "D", "E", "F"),
    common.legend = TRUE,
    legend = "bottom",
    align = "hv" # Align both horizontally and vertically
)

ggsave("results/IVCD_csr.pdf",
    units = c("cm"),
    bg = "white",
    dpi = 600,
    height = 18,
    width = 26
)

desired_legend <- get_legend(totaltiter_STD)

ggarrange(vcd, via, dia, high_glc_plt, med_glc_plt, low_glc_plt, lactate_plot, titer_plot, totaltiter_STD,
    labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
    common.legend = TRUE,
    legend = "bottom",
    align = "hv", # Align both horizontally and vertically
    legend.grob = desired_legend
)

ggsave("results/VCD_VIA_DIA_IVCD_growth_arranged.pdf",
    units = c("cm"),
    height = 27,
    width = 26,
    bg = "white",
    dpi = 600
)


