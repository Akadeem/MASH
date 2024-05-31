# General function to draw the pyramid plots ----
pyramid <- function(age.range,
                    male.morbidity,
                    female.morbidity,
                    male.mortality,
                    female.mortality,
                    y_max,
                    y_label,
                    y_step,
                    legend.pos,
                    legend.jus) {
  df_male <- data.frame(
    age_range = age.range,
    mortality = male.mortality,
    morbidity = male.morbidity
  )
  df_female <- data.frame(
    age_range = age.range,
    mortality = female.mortality,
    morbidity = female.morbidity
  )

  ## Reshape the data frames ----
  df_male_long <- df_male %>%
    pivot_longer(
      cols = -age_range,
      names_to = "category",
      values_to = "value"
    )

  df_female_long <- df_female %>%
    pivot_longer(
      cols = -age_range,
      names_to = "category",
      values_to = "value"
    )


  ## Create the plot ----
  ggplot() +
    geom_bar(data = df_female_long, aes(y = value, x = age_range, fill = interaction("Female", category)), position = "stack", stat = "identity", color = "#7F7F7F", width = 0.62) +
    geom_bar(data = df_male_long, aes(y = value, x = age_range, fill = interaction("Male", category)), position = "stack", stat = "identity", color = "#000000", width = 0.62) +
    coord_flip() +
    scale_fill_manual(
      values = c("Male.morbidity" = "#FFFFFF", "Male.mortality" = "#000000", "Female.morbidity" = "#FFFFFF", "Female.mortality" = "#7F7F7F"),
      labels = c("Male.morbidity" = "Male Morbidity", "Male.mortality" = "Male Mortality", "Female.morbidity" = "Female Morbidity", "Female.mortality" = "Female Mortality")
    ) +
    theme_minimal() +
    labs(x = "Age Group", y = y_label, fill = "") +
    guides(fill = guide_legend(override.aes = list(color = c("#7F7F7F", "#7F7F7F", "#000000", "#000000"), size = 2.5))) +
    theme(
      legend.position = legend.pos, legend.justification = legend.jus, legend.background = element_rect(fill = "white", color = "black"), legend.title = element_blank(),
      panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()
    ) +
    scale_y_continuous(breaks = seq(-y_max, y_max, by = y_step), limits = c(-y_max, y_max), labels = function(x) format(abs(x), big.mark = ",", scientific = F))
}
