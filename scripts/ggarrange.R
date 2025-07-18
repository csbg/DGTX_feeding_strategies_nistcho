


ggarrange(
    IVCD_plot, qp_entire, titer_qp, growth, qp_gluc, lactate_qp,
    labels = c("A", "B", "C", "D", "E", "F"),
    common.legend = TRUE,
    legend = "bottom",
    align = "hv" # Align both horizontally and vertically
)

ggsave("results/IVCD_csr.pdf",
    units = c("cm"),
    bg = "white",
    dpi = 600
)

ggarrange(vcd, via, dia, high_glc_plt, med_glc_plt, low_glc_plt, lactate_plot, titer_plot, totaltiter,
    labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
    common.legend = TRUE,
    legend = "bottom",
    align = "hv" # Align both horizontally and vertically
)

ggsave("results/VCD_VIA_DIA_IVCD_growth_arranged.pdf",
    units = c("cm"),
    height = 27,
    width = 26,
    bg = "white",
    dpi = 600
)


