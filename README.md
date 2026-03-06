# CLASE-06-03-2026-2
base <- read.csv("genes2.csv")
> library(dplyr)

Adjuntando el paquete: ‘dplyr’

The following objects are masked from ‘package:stats’:

    filter, lag

The following objects are masked from ‘package:base’:

    intersect, setdiff, setequal, union
> 
> data <- read.csv("genes2.csv")
> library(dplyr)
> 
> data <- read.csv("genes2.csv")
> promedios <- data %>%
+     group_by(Gene, Tratamiento) %>%
+     summarise(mean_TPM = mean(Expresion_TPM), .groups = "drop")
> control <- promedios %>%
+     filter(Tratamiento == "Control") %>%
+     select(Gene, control_TPM = mean_TPM)
> comparacion <- promedios %>%
+     filter(Tratamiento != "Control") %>%
+     left_join(control, by = "Gene")
> comparacion <- comparacion %>%
+     mutate(log2FC = log2(mean_TPM / control_TPM))
> head(comparacion)
# A tibble: 6 × 5
  Gene     Tratamiento mean_TPM control_TPM log2FC
  <chr>    <chr>          <dbl>       <dbl>  <dbl>
1 Gene_001 CRISPR          7.36        9.36 -0.347
2 Gene_001 DrugA           8.45        9.36 -0.147
3 Gene_001 DrugB          NA           9.36 NA    
4 Gene_001 RNAi            5.80        9.36 -0.689
5 Gene_002 CRISPR          7.79        7.09  0.136
6 Gene_002 DrugA          12.2         7.09  0.782
> modelo <- aov(rendimiento ~ fertilizante, data = genes2)
Error: objeto 'genes2' no encontrado

> modelo <- aov(rendimiento ~ fertilizante, data = base)
Error en eval(predvars, data, env): objeto 'rendimiento' no encontrado

> plot(res_aov, which = 3)
Error: objeto 'res_aov' no encontrado

> plot(base_aov, which = 3)
Error: objeto 'base_aov' no encontrado

> library(dplyr)
> library(ggplot2)
> base <- read.csv("genes2.csv")
> resultados <- base %>%
+     group_by(Gene) %>%
+     summarise(
+         control_mean = mean(Expresion_TPM[Tratamiento == "Control"]),
+         treat_mean = mean(Expresion_TPM[Tratamiento != "Control"]),
+         log2FC = log2(treat_mean / control_mean),
+         pvalue = t.test(
+             Expresion_TPM[Tratamiento != "Control"],
+             Expresion_TPM[Tratamiento == "Control"]
+         )$p.value
+     )
> resultados$negLog10P <- -log10(resultados$pvalue)
> resultados$Significance <- "No significativo"
> 
> resultados$Significance[resultados$log2FC > 1 & resultados$pvalue < 0.05] <- "Upregulated"
> resultados$Significance[resultados$log2FC < -1 & resultados$pvalue < 0.05] <- "Downregulated"
> ggplot(resultados, aes(x = log2FC, y = negLog10P, color = Significance)) +
+     geom_point(size = 3, alpha = 0.8) +
+     scale_color_manual(values = c(
+         "Upregulated" = "red",
+         "Downregulated" = "blue",
+         "No significativo" = "gray"
+     )) +
+     geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
+     geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
+     labs(
+         title = "Volcano Plot de Expresión Génica",
+         x = "Log2 Fold Change",
+         y = "-Log10(p-value)"
+     ) +
+     theme_minimal()
Aviso:
Removed 159 rows containing missing values or values outside the scale range (`geom_point()`). 

> resultados <- base %>%
+     group_by(Gene) %>%
+     summarise(
+         control_mean = mean(Expresion_TPM[Tratamiento == "Control"]),
+         treat_mean = mean(Expresion_TPM[Tratamiento != "Control"]),
+         log2FC = log2((treat_mean + 1) / (control_mean + 1)),
+         pvalue = t.test(
+             Expresion_TPM[Tratamiento != "Control"],
+             Expresion_TPM[Tratamiento == "Control"]
+         )$p.value
+     )
> resultados <- resultados %>%
+     filter(!is.na(log2FC) & !is.na(pvalue))
> View(base)
> sum(is.na(resultados$log2FC))
[1] 0
> sum(is.na(resultados$pvalue))
