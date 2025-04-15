# libs ----
pacman::p_load(
    tidyverse,
    ggplot2,
    ggtext,
    patchwork
)


# source lickometer library
devtools::source_url("https://github.com/lab-cpl/lickometer-library/blob/main/src/lickometer_functions_compilate.R?raw=TRUE")

# https://stackoverflow.com/questions/47044068/get-the-path-of-current-script
# get path of source file
getCurrentFileLocation <-  function()
{
    this_file <- commandArgs() %>% 
        tibble::enframe(name = NULL) %>%
        tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
        dplyr::filter(key == "--file") %>%
        dplyr::pull(value)
    if (length(this_file)==0)
    {
        this_file <- rstudioapi::getSourceEditorContext()$path
    }
    return(dirname(this_file))
}
# sets path on source file location
script_path <- getCurrentFileLocation()
setwd(script_path)

# publication theme
source("../../PHD_data/publication_theme/publication_theme.R")
fill <- palette.colors(palette = "Okabe-Ito")
color <- palette.colors(palette = "Okabe-Ito")

# get data-set ----
boop_dataset <- read_csv("../datasets/behavioral-tests/boops.csv") %>% 
    group_by(ID) %>% 
    mutate(
        date = as.numeric(as.factor(date)),
        week = as.numeric(as.factor(week))
    ) %>% 
    filter(
        date <= 59
    )

boop_time_dataset <- read_csv("../datasets/behavioral-tests/boops_time.csv") %>% 
    group_by(ID) %>% 
    mutate(
        date = as.numeric(as.factor(date)),
        week = as.numeric(as.factor(week))
    ) %>% 
    filter(
        date <= 59
    )


# model 1: boops per group ----

# data is discrete and positive, so I model conditional expectation
# with either poisson or non-negative binomial
mdl1 <- lme4::glmer.nb(
    data = boop_dataset,
    failed_retrievals ~ group * experimental_phase + date + (1|ID)
)
summary(mdl1)

emm1 <-
    emmeans::emmeans(
    mdl1, 
    pairwise ~ group | experimental_phase + date,
    regrid = "response"
)

emm1_means <- as_tibble(emm1$emmeans)
emm1_contrasts <- as_tibble(emm1$contrasts)
emm1_effsize <- as_tibble(emmeans::eff_size(emm1, edf = df.residual(mdl1), sigma = sigma(mdl1)))

# this is the raw data
emm1_raw <-
    boop_dataset %>% 
    group_by(group, experimental_phase, ID) %>% 
    summarise(
        estimate = (failed_retrievals)
    ) %>% 
    rename(response = estimate)
emm1_raw

annotation_df <- data.frame(
    group = c("Uncertainty"),
    start = c("Base"),
    end = c("Exp"),
    y = c(log10(30)),
    label = c("")
)

add_logticks  <- function (base = 10, sides = "bl", scaled = TRUE, 
                           short = unit(0.1, "cm"), mid = unit(0.2, "cm"),  long = unit(0.3, "cm"), 
                           colour = "black",  size = 0.5, linetype = 1, alpha = 1, color = NULL, 
                           data =data.frame(x = NA),... )   {
    if (!is.null(color)) 
        colour <- color
    layer(geom = "logticks", params = list(base = base, 
                                           sides = sides, scaled = scaled, short = short, 
                                           mid = mid, long = long, colour = colour, size = size, 
                                           linetype = linetype, alpha = alpha, ...), 
          stat = "identity", data = data , mapping = NULL, inherit.aes = FALSE, position = "identity",
          show.legend = FALSE)
}

p1 <-
    emm1_raw %>% 
    mutate(
        group = if_else(group == "No-Uncertainty", "Control", "Uncertainty"),
        experimental_phase = if_else(
            experimental_phase == "Baseline", "Base", "Exp"
    )) %>% 
    group_by(
        ID, group, experimental_phase
    ) %>% 
    summarise(
        response = mean(response)
    ) %>% 
    ggplot(aes(
        experimental_phase, response
    )) +
    geom_boxplot(
        width = 0.5,
        outlier.shape = NA,
        fill = NA,
        lwd = 0.75,
        aes(color = group)
    ) +
    geom_line(aes(group = ID, color = group)) +
    geom_point(
        size = 5,
        shape = 21,
        aes(fill = group),
        color = "black",
        alpha = 0.5
    ) +
    ggsignif::geom_signif(
        data = annotation_df,
        inherit.aes = FALSE,
        manual = TRUE,
        tip_length = c(0.75, 0.15),
        comparisons = list(c("Base", "Exp")),
        aes(xmin = start, xmax = end, annotations = label, y_position = y),
        map_signif_level = TRUE,
        size = 1,
        textsize = 6,
        vjust = 0.5,
        color = "black"
    ) +
#    ylab(latex2exp::TeX(r'(Failed retrievals $(\log_{10}) )')) +
    ylab(latex2exp::TeX(r'(Failed ret.)')) +
    xlab("") +
    ggpubr::theme_pubr() +
    theme(
        text = element_text(size = 30),
        axis.text=element_text(size=24),
        legend.position = "none",
        strip.background = element_blank(),
        strip.placement = "outside"
    ) +
    scale_y_log10() +
    add_logticks(side = 'l', data = data.frame(x= NA, group = 'Control')) +
    scale_fill_manual(values = c("black", "orange")) +
    scale_color_manual(values = c("black", "orange")) +
    facet_wrap(~group)
p1

# model 2: boops times ----
boop_time_ret <- boop_time_dataset %>% 
        filter(
            ret > 0, ret < 1800,
            delay %in% c(15, 60, 120, 180, 240, 300)
            ) %>% 
        group_by(ID, group, experimental_phase) %>% 
        mutate(
            date = scales::rescale(date, to = c(0, 1)),
            week = scales::rescale(week, to = c(0, 1)),
            ret = abs(ret)
            )

mdl2 <- lme4::glmer(
    data = boop_time_ret %>% 
        filter(experimental_phase == "Experimental"),
    ret ~ group * week + (1 | ID),
    family = Gamma(link = "log")
)
summary(mdl2)

emm2 <- emmeans::emtrends(
    mdl2,
    pairwise ~ group,
    var = "week",
    regrid = "response"
) %>% summary(., infer = TRUE)
emm2

emm2_plot <- emmeans::emmeans(
    mdl2,
    ~ group * week,
    regrid = "response",
    at = list(week = seq(0, 1, length.out = 8))
) %>% 
    as_tibble()
emm2_plot

ret_time_daily_mean <- boop_time_ret %>%
    filter(experimental_phase == "Experimental") %>% 
    group_by(ID, group, week) %>% 
    summarise(ret = mean(ret)) 
ret_time_daily_mean



p2 <- emm2_plot %>% 
    ggplot(aes(
        week*8, response, color = group, fill = group
    )) +
    geom_point(
        data = ret_time_daily_mean,
        aes(week*8, ret),
        alpha = 0.5,
        shape = 21,
        size = 5,
        color = "black"
    ) +
    geom_ribbon(
        aes(ymin=asymp.LCL, ymax=asymp.UCL, fill = group),
        alpha = 0.75,
        color = "black"
    ) +
    geom_line(color = "black") +
    annotate(
        geom = "text",
        x = 4,
        y = 250,
        label = "*",
        size = 10
    ) +
    ggpubr::theme_pubr() +
    guides(
        color = guide_legend(override.aes = list(color = NA, fill = NA)),
        fill = "none"
    ) +
    ylab(latex2exp::TeX(r'(Ret. latency)')) +
    xlab("Weeks") +
    scale_y_log10() +
    annotation_logticks(side = "l") +
    theme(
        text = element_text(size = 30),
        axis.text=element_text(size=24),
        legend.position = "none",
        strip.background = element_blank(),
        strip.placement = "outside"
    ) +
    scale_fill_manual(values = c("gray", "orange")) +
    scale_color_manual(values = c("gray", "orange"))
p2

free(p1) + p2 +
    plot_layout(widths = c(1, 1)) +
    plot_annotation(tag_levels = "A")

