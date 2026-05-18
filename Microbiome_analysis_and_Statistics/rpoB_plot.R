library(readxl)
rpoB_nuevos <- read_excel("Documents/Doctorado/Estudios/INCLIVA_2024/rpoB/rpoB_nuevos.xlsx", 
                          sheet = "all")
View(rpoB_nuevos)


library(tidyverse)
library(ggh4x)

# Prepare data
rpoB_nuevos <- rpoB_nuevos %>%
  mutate(
    Patient = str_extract(Muestra, "^[^/]+"),
    Rifaximin = factor(Rifaximin, levels = c("Before", "After")),
    Mutant = factor(Mutant, levels = c("Negative", "Positive"))
  )

# Strip colors
strip_df <- rpoB_nuevos %>%
  distinct(Patient, Grupo)

# Plot
ggplot(rpoB_nuevos,
       aes(x = Rifaximin,
           fill = Mutant)) +
  
  geom_bar(
    position = "dodge",
    alpha = 0.5
  ) +
  
  facet_wrap2(
    ~ Patient,
    strip = strip_themed(
      background_x = elem_list_rect(
        fill = ifelse(
          strip_df$Grupo == "R",
          alpha("#009E73", 0.3),  # responders
          alpha("#CC79A7", 0.3)   # non-responders
        ))
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "Negative" = "#88CCEE",
      "Positive" = "#DDCC77"
    ),
    name = "Mutant status"
  )  +
  
  theme_bw() +
  
  labs(
    title = "Mutant status before and after rifaximin by patient",
    x = "Rifaximin",
    y = ""
  )+
  coord_cartesian(ylim = c(0, 1))+ 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) 

