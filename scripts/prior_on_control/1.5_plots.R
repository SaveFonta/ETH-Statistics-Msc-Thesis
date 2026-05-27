

library(RBesT)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)
library(stringr)

source("scripts/prior_on_control/00.functions.MC.R")
all_results <- readRDS("data/1.6.rds")


euii <- lapply(all_results, function(dec.func) {
  lapply(dec.func, function(job) {
    compute_euii(job)
  })
})


# trasform "ap.MAP_dp.Robust_0.20" -> list(ap, dp) 
parse_job_name <- function(job_name) {
  ap <- sub("^ap\\.(.+)_dp\\..+$", "\\1", job_name)
  dp <- sub("^ap\\..+_dp\\.(.+)$", "\\1", job_name)
  list(ap = ap, dp = dp)
}


# -- Flatten across all decision lists / jobs / deltas --
flatten_power <- function(all_results) {
  bind_rows(lapply(names(all_results), function(dec_name) {
    dec_res <- all_results[[dec_name]]
    bind_rows(lapply(names(dec_res), function(job_name) {
      parsed <- parse_job_name(job_name)
      job    <- dec_res[[job_name]]        # named list: delta.0, delta.10, ...
      bind_rows(lapply(names(job), function(delta_name) {
        d  <- as.numeric(sub("delta\\.", "", delta_name))
        ov <- job[[delta_name]]$Overall
        data.frame(
          decision_list  = dec_name,
          analysis_prior = parsed$ap,
          design_prior   = parsed$dp,
          delta          = d,
          as.data.frame(t(ov)),
          check.names    = FALSE
        )
      }))
    }))
  }))
}




all_results$decision_list$ap.MAP_dp.MAP$delta.90$Overall  |> t()  |>  as.data.frame()

# all_results
#      -----------decision_list
#                      --------------# ap.MAP_dp.Vague
                     #                      --------------- delta.0
                     #                                         ----------------- > $Overall
                     #                      --------------- delta.1
                     #                                         ----------------- > $Overall
#                        ----              ...  
#      -----------decision_list2
#                         ---------------




# -- Flatten EUII across all decision lists / jobs / prior_H1 values --
flatten_euii <- function(euii) {
  bind_rows(lapply(names(euii), function(dec_name) {
    dec_res <- euii[[dec_name]]
    bind_rows(lapply(names(dec_res), function(job_name) {
      parsed   <- parse_job_name(job_name)
      job_euii <- dec_res[[job_name]]      # named list: "0.01", "0.1", "0.5"
      bind_rows(lapply(names(job_euii), function(prior_h1) {
        job_euii[[prior_h1]] |>
          mutate(
            decision_list  = dec_name,
            analysis_prior = parsed$ap,
            design_prior   = parsed$dp,
            prior_H1       = as.numeric(prior_h1)
          )
      }))
    }))
  }))
}

df_power <- flatten_power(all_results)
df_euii  <- flatten_euii(euii)


# Not use skeptical as analysiis prior since this is what Best did, but 
#also cause we would never use the skeptical as historical
df_power <- df_power  |> filter (analysis_prior != "Skeptical")
df_euii <- df_euii  |> filter (analysis_prior != "Skeptical")




# rename different decisions
dec_labels <- c(
  "decision_list"             = "Classic",
  "decision_list_nofut"       = "No Futility",
  "decision_list_norel"       = "No Relevance",
  "decision_list_nofut.norel" = "No Futility No Relevance"
)

DEC_LEVELS <- c("Classic", "No Futility", "No Relevance", "No Futility No Relevance")

df_power <- df_power |>
  mutate(decision_list = dplyr::recode(decision_list, !!!dec_labels),
         decision_list = factor(decision_list, levels = DEC_LEVELS))  |> 
         select(decision_list, analysis_prior, design_prior, delta, Power)

df_euii <- df_euii |>
  mutate(decision_list = dplyr::recode(decision_list, !!!dec_labels),
         decision_list = factor(decision_list, levels = DEC_LEVELS))

AP_LEVELS  <- c("MAP", "Normal", "Robust_0.20", "Robust_0.5", "Vague")
DP_LEVELS  <- c("MAP", "Normal", "Robust_0.20", "Robust_0.5", "Vague", "Skeptical")


df_power <- df_power |>
  mutate(analysis_prior = factor(analysis_prior, levels = AP_LEVELS),
         design_prior   = factor(design_prior,   levels = DP_LEVELS))

df_euii <- df_euii |>
  mutate(analysis_prior = factor(analysis_prior, levels = AP_LEVELS),
         design_prior   = factor(design_prior,   levels = DP_LEVELS))

cat("Power rows:", nrow(df_power), "\n") 
cat("EUII rows :", nrow(df_euii),  "\n")


df_classic <- df_power  |> filter (
  decision_list == "Classic", delta != 90)


# plot to put in manuscript
df_classic |>
  ggplot(aes(x = design_prior, y = analysis_prior, fill = Power)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", Power)), size = 2.8) +
  scale_fill_gradient2(
    low      = "#2E6DA4",
    mid      = "white",
    high     = "#E04F39",
    midpoint = 0.5,
    labels   = scales::percent_format(accuracy = 0.1),
    name     = "Power"
  ) +
  facet_wrap(~ delta, ncol = 3,
             labeller = labeller(delta = function(x) paste0("Delta = ", x))) +
  labs(  x = "Design Prior",
    y = "Analysis Prior"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 35, hjust = 1),
    panel.grid       = element_blank(),
    strip.background = element_rect(fill = "grey92")
  )


df_classic |>
  ggplot(aes(
    x = delta,
    y = Power,
    colour = analysis_prior,
    group = analysis_prior
  )) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(
    ~ design_prior,
    ncol = 3
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 0.1)
  ) +
  labs(
    x = "Delta",
    y = "Power",
    colour = "Analysis Prior"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92"),
    legend.position = "bottom"
  )







# =============================================================================
# POWER
#    One plot per analysis prior (6 total)
#    facet for decision, 4 panels per plot
#    colour = design_prior
# =============================================================================

plot_power <- function(df_power, fixed_ap) {
  df_power |>
    filter(analysis_prior == fixed_ap) |>
    ggplot(aes(x = delta, y = Power,
               colour = design_prior, group = design_prior)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.5) +
    geom_hline(yintercept = 0.80, linetype = "dashed",
               colour = "grey40", linewidth = 0.5) +
    facet_wrap(~ decision_list, ncol = 2) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_x_continuous(breaks = seq(0, 90, 20)) +
    labs(
      title    = paste("Power \u2014 Analysis prior:", fixed_ap),
      subtitle = "Facets = decision list  |  Colour = design prior",
      x        = "Delta",
      y        = "Power",
      colour   = "Design prior"
                ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position  = "bottom"
    )
}
power_plots <- lapply(AP_LEVELS, function(ap) {
  plot_power(df_power, ap)
}) |> setNames(AP_LEVELS)


# =============================================================================
# EUII
# SAMe as power, but linetype is PR h1
# =============================================================================

SHOW_PRIOR_H1_LINETYPE <- FALSE
PRIOR_H1_REF <- 0.1


plot_euii <- function(df_euii, fixed_ap, show_linetype = SHOW_PRIOR_H1_LINETYPE) {
  df <- df_euii |> filter(analysis_prior == fixed_ap)

  if (!show_linetype) {
    # keep only the reference prior_H1, one line per design prior
    df <- df |> filter(prior_H1 == PRIOR_H1_REF)
  }

  df <- df |> mutate(prior_H1 = factor(prior_H1))

  if (show_linetype) {
    p <- ggplot(df, aes(x = Delta, y = EUII,
                        colour   = design_prior,
                        linetype = prior_H1,
                        group    = interaction(design_prior, prior_H1)))
  } else {
    p <- ggplot(df, aes(x = Delta, y = EUII,
                        colour = design_prior,
                        group  = design_prior))
  }

  p <- p +
    geom_line(linewidth = 0.85) +
    geom_hline(yintercept = 1, linetype = "dashed",
               colour = "grey40", linewidth = 0.5) +
    facet_wrap(~ decision_list, ncol = 2) +
    scale_x_continuous(breaks = seq(0, 90, 20)) +
    labs(
      title  = paste("EUII for Analysis prior:", fixed_ap),
      x      = "Delta",
      y      = "EUII",
      colour = "Design prior"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position  = "bottom",
      legend.box       = "vertical"
    )

  if (show_linetype) {
    p <- p +
      scale_linetype_manual(
        values = c("0.01" = "dotted", "0.1" = "solid", "0.5" = "longdash"),
        labels = c("0.01" = "P(H\u2081) = 0.01",
                   "0.1"  = "P(H\u2081) = 0.10",
                   "0.5"  = "P(H\u2081) = 0.50")
      ) +
      labs(
        subtitle = "Facets = decision list  |  Colour = design prior  |  Linetype = P(H\u2081)",
        linetype = "Prior on H\u2081"
      )
  } else {
    p <- p +
      labs(subtitle = paste0("Facets = decision list  |  Colour = design prior  |  P(H\u2081) = ",
                             PRIOR_H1_REF))
  }

  p
}

cat("Building EUII plots (one per analysis prior)...\n")
euii_plots <- lapply(AP_LEVELS, function(ap) {
  plot_euii(df_euii, ap)
}) |> setNames(AP_LEVELS)


# =============================================================================
# TYPE I ERROR PLOT
# =============================================================================

t1e_df <- df_power |>
  filter(delta == 0) |>
  select(decision_list, analysis_prior, design_prior, Power) |>
  rename(T1E = Power)

plot_t1e <- t1e_df |>
  ggplot(aes(x = analysis_prior, y = T1E,
             colour = design_prior, group = design_prior)) +
  geom_hline(yintercept = 0.05, linetype = "dashed",
             colour = "#ee4a4a", linewidth = 0.4) +
  geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
  facet_wrap(~ decision_list, ncol = 2) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(
    subtitle = "x = analysis prior  |  Colour = design prior  |  Facets = decisions",
    x        = "Analysis prior",
    y        = "Type I Error",
    colour   = "Design prior"  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "bottom",
    axis.text.x      = element_text(angle = 30, hjust = 1)
  )

cat("Building T1E plot...\n")







print(plot_t1e)
readline(prompt = "Press [Enter] for power plots...")

for (ap in AP_LEVELS) {
  cat("Analysis prior:", ap, "\n")
  print(power_plots[[ap]])
  readline(prompt = "Press [Enter] for next plot...")
}

for (ap in AP_LEVELS) {
  cat("Analysis prior:", ap, "\n")
  print(euii_plots[[ap]])
  readline(prompt = "Press [Enter] for next plot...")
}








# diagonal function
matched <- function(df) {
  df |> filter(analysis_prior == design_prior)
}

# extract each decision list 
dec_base   <- df_power |> filter(decision_list == "No Futility No Relevance")
dec_nofut  <- df_power |> filter(decision_list == "No Futility")
dec_norel  <- df_power |> filter(decision_list == "No Relevance")
dec_full   <- df_power |> filter(decision_list == "Classic")

decomp_df <- matched(dec_base) |>
  rename(Power_base = Power) |>
  left_join(
    matched(dec_nofut) |> rename(Power_dc = Power),
    by = c("analysis_prior", "design_prior", "delta")
  ) |>
  left_join(
    matched(dec_norel) |> rename(Power_fut = Power),
    by = c("analysis_prior", "design_prior", "delta")
  ) |>
  left_join(
    matched(dec_full) |> rename(Power_full = Power),
    by = c("analysis_prior", "design_prior", "delta")
  ) |>
  mutate(
    Cost_DC   = Power_base - Power_dc,    # cost of adding DC criterion
    Cost_Fut  = Power_base - Power_fut,   # cost of adding futility
    Cost_Both = Power_base - Power_full   # cost of adding both
  ) |>
  select(analysis_prior, delta, Cost_DC, Cost_Fut, Cost_Both) |>
  pivot_longer(cols = starts_with("Cost"),
               names_to  = "Component",
               values_to = "Power_Loss") |>
  mutate(Component = recode(Component,
    Cost_DC   = "DC criterion",
    Cost_Fut  = "Futility stopping",
    Cost_Both = "Both combined"
  ))

p_decomp <- decomp_df |>
  ggplot(aes(x = delta, y = Power_Loss,
             colour = Component, group = Component)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.6) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey40", linewidth = 0.4) +
  scale_colour_manual(values = c(
    "DC criterion"      = "#E04F39",
    "Futility stopping" = "#2E6DA4",
    "Both combined"     = "#8B4FA8"
  ), name = "Power cost of:") +
  scale_x_continuous(breaks = unique(decomp_df$delta)) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  facet_wrap(~ analysis_prior, ncol = 3) +
  labs(
    subtitle = "Baseline = No Futility No Relevance | Matched prior diagonal",
    x        = "Delta",
    y        = "Power loss vs baseline"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )

print(p_decomp)
