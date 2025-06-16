# lib load ----
pacman::p_load(
    tidyverse,
    ggplot2,
    patchwork,
    robustlmm,
    magick,
    figpatch,
    DiagrammeR,
    ggforce,
    ggthemes,
    mgcv
)
setwd(this.path::here())

# Plot stuff ----
theme_uncertainty <-
    theme_par() +
    update_geom_defaults("point", list(size = 5, alpha = 0.5, shape = 21)) +
    update_geom_defaults("pointrange", list(linewidth = 1.5, size = 1.5)) +
    theme(
        text = element_text(size = 24),
        axis.text = element_text(size = 14),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        legend.position = "none"
    )

boxplot_sig_bracket <- function(group1, group2) {
    ggsignif::geom_signif(
        comparisons = list(c(group1, group2)),
        map_signif_level = TRUE,
        textsize = 0
    )
}


# Weight (w) ----
wd <- read_csv("../datasets/weights/weights.csv") %>%
    group_by(ID) %>%
    mutate(
        rel_date = as.numeric(date - min(date))
    )

# same but weekly aggregate
weekly_weight <- wd %>%
    ungroup() %>%
    group_by(ID) %>%
    mutate(
        rel_date = cumsum(replace_na(as.numeric(date - lag(date)), 0)),
        days_from_end = as.numeric(abs(date - max(date))),
        phase = if_else(days_from_end <= 56, "uncertainty", "baseline"),
        week = as.numeric(as.factor(lubridate::isoweek(date))),
        week = if_else(phase == "uncertainty", week + 1, week)
    ) %>%
    ungroup() %>%
    group_by(ID, phase, week, group) %>%
    summarise(
        weight = mean(weight)
    ) %>%
    ungroup() %>%
    group_by(ID) %>%
    mutate(
        rel_weight = weight - weight[week == 1]
    )

## mdls ----

### weight @ 10 weeks ----
w_mdl1 <- lme4::lmer(
    data = weekly_weight %>% filter(phase == "uncertainty"),
    rel_weight ~ group * week + (1 | ID)
)

w_mdl1_emm <- emmeans::emmeans(
    w_mdl1,
    pairwise ~ group | week,
    type = "response",
    at = list(week = c(10))
)
w_mdl1_emm
w_mdl1_emm$contrasts %>% broom.mixed::tidy(., conf.int = TRUE)

### weight trends ----
w_mdl1_trend <- emmeans::emtrends(
    w_mdl1,
    pairwise ~ group * week,
    var = "week"
)
w_mdl1_trend
w_mdl1_trend$emtrends %>% broom.mixed::tidy(., conf.int = TRUE)



## p::weekly weight ----
wp1 <- weekly_weight %>%
    filter(week <= 10) %>%
    ggplot(aes(
        week, rel_weight,
        group = ID
    )) +
    geom_vline(xintercept = 2.5, linetype = "dashed") +
    geom_line(alpha = 0.25, aes(color = group)) +
    stat_summary(
        fun.data = "mean_se",
        geom = "pointrange",
        aes(group = interaction(group, week), color = group)
    ) +
    stat_summary(
        fun.data = "mean_se",
        geom = "line",
        aes(group = group, color = group)
    ) +
    scale_x_continuous(breaks = seq(1, 10, 1)) +
    theme_uncertainty +
    scale_y_continuous(
        breaks = seq(-10, 10, 2),
        limits = c(-10, 10),
        expand = c(0, 0)
    ) +
    ylab(expression(Delta * " Weekly weight (gr.)")) +
    xlab("Weeks") +
    scale_fill_manual(values = c("black", "orange")) +
    scale_color_manual(values = c("black", "orange"))
wp1

## p::weight @ week10 ----
wp2 <- weekly_weight %>%
    filter(week == 10) %>%
    ggplot(aes(
        group, weight,
        color = group
    )) +
    geom_boxplot(outlier.shape = NA, width = 0.5, aes(color = group)) +
    geom_point(aes(fill = group), color = "black") +
    ggsignif::geom_signif(
        comparisons = list(c(1, 2)),
        textsize = 0,
        y_position = 32
    ) +
    theme_uncertainty +
    scale_x_discrete(labels = c("LU", "HU")) +
    scale_y_continuous(
        breaks = seq(20, 35, 5),
        limits = c(20, 35),
        expand = c(0, 0)
    ) +
    ylab("Final weight (gr.)") +
    xlab("Weeks") +
    scale_fill_manual(values = c("black", "orange")) +
    scale_color_manual(values = c("black", "orange"))
wp2

## p:weight trends
wp3 <- w_mdl1_trend$emtrends %>%
    broom.mixed::tidy(., conf.int = TRUE) %>%
    ggplot(aes(
        group, week.trend,
        ymin = conf.low, ymax = conf.high
    )) +
    geom_pointrange(aes(color = group), size = 1) +
    geom_hline(yintercept = 0, color = "gray") +
    boxplot_sig_bracket(1, 2) +
    scale_fill_manual(values = c("black", "orange")) +
    scale_color_manual(values = c("black", "orange")) +
    theme_uncertainty +
    scale_x_discrete(labels = c("LU", "HU")) +
    scale_y_continuous(
        breaks = seq(0, 0.6, 0.1),
        limits = c(0, 0.6),
        expand = c(0, 0)
    ) +
    ylab(latex2exp::TeX(r"($\hat{\mu}_{weight/weeks}$)")) +
    xlab("")
wp3

# Intake (I) ----
id <- read_csv("../datasets/intake/intake.csv")

# data for metabolic efficiency, weight interpolated
id_wd <- left_join(id %>% ungroup(), wd %>% ungroup()) %>%
    filter(experimental_phase == "Experimental") %>%
    ungroup() %>%
    group_by(ID) %>%
    arrange(date, .by_group = TRUE) %>%
    mutate(
        rel_date = as.numeric(date - min(date)),
        weeks = cumsum(as.numeric(rel_date %% 7 == 0)),
        interp_weight = zoo::na.approx(weight, rel_date, na.rm = FALSE),
        cum_intake = cumsum(intake) - head(intake, n = 1)
    ) %>%
    drop_na(interp_weight) %>%
    mutate(
        delta_weight = interp_weight - head(interp_weight, n = 1)
    )
id_wd

## mdls ----

### pellet intake ----
pellet_mdl <- lmerTest::lmer(
    data = id_wd,
    intake ~ rel_date * group + (rel_date | ID),
    control = lme4::lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
    )
)
summary(pellet_mdl)

pellet_emm <- emmeans::emmeans(
    pellet_mdl,
    pairwise ~ group | rel_date,
    at = list(rel_date = c(56))
)
pellet_emm
broom.mixed::tidy(pellet_emm$emmeans, conf.int = TRUE)
broom.mixed::tidy(pellet_emm$contrasts, conf.int = TRUE)

### energy efficiency ----
LU_mdl <- lmerTest::lmer(
    data = id_wd %>% filter(group == "No-Uncertainty"),
    delta_weight ~ cum_intake * rel_date + (rel_date | ID),
    control = lme4::lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
    )
)
summary(LU_mdl)
performance::icc(LU_mdl)

# get data set with only HU
HU_data <- id_wd %>%
    filter(
        group == "Uncertainty"
    )

# fit LU into HU data
HU_vs_LU_mdl <- HU_data %>%
    ungroup() %>%
    mutate(
        preds = predict(LU_mdl, newdata = HU_data, re.form = NA),
        resid_hu = delta_weight - preds
    )
HU_vs_LU_mdl

# model analysis
resid_mdl <- lmerTest::lmer(
    data = HU_vs_LU_mdl %>% filter(experimental_phase == "Experimental"),
    resid_hu ~ rel_date + (rel_date | ID),
    control = lme4::lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
    )
)
summary(resid_mdl)
broom.mixed::tidy(resid_mdl, conf.int = TRUE)
emmeans::emmeans(resid_mdl, ~rel_date,
    at = list(rel_date = c(56)),
    type = "response"
) %>%
    broom.mixed::tidy(conf.int = TRUE)

### failed retrievals ----
frd <- read_csv("../datasets/behavioral-tests/boops.csv") %>%
    group_by(ID) %>%
    mutate(
        rel_date = as.numeric(date - min(date))
    )

# this model takes a lot of time so I saved it in rds format
# frd_mdl <- lme4::glmer.nb(
#    data = frd,
#    failed_retrievals ~ group * rel_date * experimental_phase + (rel_date|ID),
#    control = lme4::glmerControl(
#        optimizer = "bobyqa",
#        optCtrl = list(maxfun = 2e5)
#    )
# )
frd_mdl <- read_rds("../datasets/frd_mdl.rds")
summary(frd_mdl)

frd_emm <- emmeans::emmeans(
    frd_mdl,
    pairwise ~ group | rel_date * experimental_phase,
    at = list(rel_date = c(56), experimental_phase = c("Experimental")),
    type = "response"
)
frd_emm
broom.mixed::tidy(frd_emm$emmeans, conf.int = TRUE)
broom.mixed::tidy(frd_emm$contrasts, conf.int = TRUE)

frd_emm_trend <- emmeans::emtrends(
    frd_mdl,
    pairwise ~ group | rel_date * experimental_phase,
    at = list(experimental_phase = c("Experimental"), rel_date = c(28)),
    var = "rel_date",
    type = "response"
)
frd_emm_trend
broom.mixed::tidy(frd_emm_trend$emtrends, conf.int = TRUE)
broom.mixed::tidy(frd_emm_trend$contrasts, conf.int = TRUE)

### retrieval latency ----
rld <- read_csv("../datasets/behavioral-tests/retrieval_latency.csv") %>%
    filter(ret > 0, ret <= 1378) %>%
    group_by(ID) %>%
    mutate(
        rel_date = as.numeric(date - min(date))
    )

rld_mdl <- lme4::glmer(
    data = rld,
    ret ~ group * experimental_phase * rel_date + (1 | ID),
    family = Gamma(link = "log"),
    control = lme4::glmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
    )
)
summary(rld_mdl)

rld_emm <- emmeans::emmeans(
    rld_mdl,
    revpairwise ~ group | experimental_phase * rel_date,
    at = list(experimental_phase = c("Experimental"), rel_date = c(56)),
    type = "response"
)
rld_emm
broom.mixed::tidy(rld_emm$emmeans, conf.int = TRUE)
broom.mixed::tidy(rld_emm$contrasts, conf.int = TRUE)

rld_emm_trend <- emmeans::emtrends(
    rld_mdl,
    pairwise ~ group | experimental_phase * rel_date,
    at = list(experimental_phase = c("Experimental"), rel_date = c(30)),
    var = "rel_date",
    type = "response"
)
rld_emm_trend
broom.mixed::tidy(rld_emm_trend$emtrends, conf.int = TRUE)
broom.mixed::tidy(rld_emm_trend$contrasts, conf.int = TRUE)

### demand curve ----
devtools::source_url("https://github.com/lab-cpl/lickometer-library/blob/main/src/lickometer_functions_compilate.R?raw=TRUE")

lickometer_data <- load_experiment(
    metadataFileName = "../datasets/lickometer/metadata/lickometer_metadata.csv",
    data_directory_path = "../datasets/lickometer/raw"
)

write_csv(x = lickometer_data, file = "../datasets/behavioral-tests/lickometer_dc.csv")



alphaq0d <- read_csv("../datasets/behavioral-tests/alpha_q0.csv")
dcd <- read_csv("../datasets/behavioral-tests/demand_curve_data.csv")

alpha_mdl <- glm(
    data = alphaq0d %>% filter(term == "alpha"),
    estimate ~ group,
    family = Gamma(link = "log")
)
summary(alpha_mdl)

alpha_emm <- emmeans::emmeans(
    alpha_mdl,
    pairwise ~ group,
    type = "response"
)
alpha_emm
broom.mixed::tidy(alpha_emm$emmeans, conf.int = TRUE)
broom.mixed::tidy(alpha_emm$contrasts, conf.int = TRUE)

q0_mdl <- glm(
    data = alphaq0d %>% filter(term == "q0"),
    estimate ~ group,
    family = Gamma(link = "log")
)
summary(q0_mdl)

q0_emm <- emmeans::emmeans(
    q0_mdl,
    pairwise ~ group,
    type = "response"
)
q0_emm
broom.mixed::tidy(q0_emm$emmeans, conf.int = TRUE)
broom.mixed::tidy(q0_emm$contrasts, conf.int = TRUE)

# experimental this is the correlation between alpha and q0
alpha_q0_corr <- glm(
    data = alphaq0d %>%
        ungroup() %>%
        select(ID, group, term, estimate) %>%
        pivot_wider(values_from = "estimate", names_from = "term"),
    alpha ~ q0 * group,
    family = Gamma(link = "log")
)
summary(alpha_q0_corr)


alpha_q0_emtrend <- emmeans::emtrends(
    alpha_q0_corr,
    pairwise ~ group | q0,
    var = "q0",
    type = "response"
)
alpha_q0_emtrend
broom::tidy(alpha_q0_emtrend$emtrends, conf.int = TRUE)
broom::tidy(alpha_q0_emtrend$contrasts, conf.int = TRUE)

### RNA-seq ----
# rnaseq_rank is the variable importance derived from random forest
# the question it answer is: how importance is gene X in determining the class (uncertainty)
# of any given sample, we used mean decrease in gini as "importance".
# The model was asked this 10 times (iteration column is the index), then
# column rank has the actual rank obtained, all other variables are descriptors of gene X
rnaseq_rank <- read_csv("../datasets/rna-seq/data_rank.csv")
# rnaseq_pdp are the results of the partial dependence plots, basically this is a simulation
# asking random forest how the probability of a sample being uncertainty group
# moves when the gene expression changes, the x-axis is within the actuall values present
# in all samples and the y-axis (yhat) is a probability
rnaseq_pdp <- read_csv("../datasets/rna-seq/data_rank_pdp.csv")

# generate lists of genes
# each list is an ordered list with the rank of genes
rank_list <- rnaseq_rank %>%
    select(symbol, rank, iteration) %>%
    group_by(iteration) %>%
    arrange(rank, .by_group = TRUE) %>%
    group_split() %>%
    map(., function(it) {
        return(c(it$symbol))
    }) %>%
    set_names(paste0("iteration", unique(rnaseq_rank$iteration)))
rank_list

# this procedure ask the question if given the iterations, is the typical rank
# similar to random, the lower the p-value the better, the higher the rank
rra_mdl <- RobustRankAggreg::aggregateRanks(
    glist = rank_list,
    N = 26
)
rra_mdl

rra_res <- tibble(rra_mdl) %>%
    filter(Score != 1)

# get slopes from pdp
pdp_slopes <- rnaseq_pdp %>%
    filter(symbol %in% rra_res$Name) %>%
    group_by(symbol) %>%
    group_split() %>%
    map_dfr(., function(X) {
        lm(data = X, yhat ~ x + I(x^3)) %>%
            broom::tidy(., conf.int = TRUE) %>%
            filter(term == "x") %>%
            mutate(gene = X$symbol[1])
    })
pdp_slopes

### Clarity ----

#### intensity ----

# this is the complete clarity data with all classifications
clarity_zpos <- read_csv("../datasets/clarity/z_pos_readings.csv")

# pos_z_lim is the main feature to make the anterior-posterior axis comparable
all_neuron_reading <- clarity_zpos %>%
    mutate(
        pos_z_lim = scales::rescale(pos_z, to = c(0, 1))
    ) %>%
    group_by(animal_id, filename) %>%
    mutate(
        rel_distance = pos_z - min(pos_z),
        scaled_rel_distance = scale(rel_distance),
        scaled_int = scale(mean_intensity),
        mean_int_rel = mean_intensity / max(mean_intensity),
        integrated_in = mean_intensity * Area / max(mean_intensity * Area)
    )

# filter only neurons
neuron_reading <- clarity_zpos %>%
    filter(predictions == "neuron") %>%
    mutate(
        pos_z_lim = scales::rescale(pos_z, to = c(0, 1))
    ) %>%
    group_by(animal_id, filename) %>%
    mutate(
        rel_distance = pos_z - min(pos_z),
        scaled_rel_distance = scale(rel_distance),
        scaled_int = scale(mean_intensity),
        mean_int_rel = mean_intensity / max(mean_intensity),
        integrated_in = mean_intensity * Area / max(mean_intensity * Area)
    )

# mean intensity is highly non-linear across the anterior-posterior axis
# here's is the statistically justification of fitting a polynomial model
# going from the null model to higher complexity of models
# null model
clarity_mdl0 <- lme4::lmer(
    data = neuron_reading,
    mean_int_rel ~ 1 + (1 | animal_id),
    REML = TRUE,
    control = lme4::lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
    )
)
# null model with all random effects
clarity_mdl1 <- lme4::lmer(
    data = neuron_reading,
    mean_int_rel ~ 1 + (pos_z_lim | animal_id),
    REML = TRUE,
    control = lme4::lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
    )
)
# base model
clarity_mdl2 <- lme4::lmer(
    data = neuron_reading,
    mean_int_rel ~ pos_z_lim * group + (pos_z_lim | animal_id),
    REML = TRUE,
    control = lme4::lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
    )
)
# adding squared term
clarity_mdl3 <- lme4::lmer(
    data = neuron_reading,
    mean_int_rel ~ pos_z_lim * I(pos_z_lim^2) * group + (pos_z_lim | animal_id),
    REML = TRUE,
    control = lme4::lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
    )
)
# complete model
clarity_mdl4 <- lme4::lmer(
    data = neuron_reading,
    mean_int_rel ~ pos_z_lim * I(pos_z_lim^2) * group + (pos_z_lim | filename) + (pos_z_lim | animal_id),
    REML = TRUE,
    control = lme4::lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
    )
)
# complete model cubed term
clarity_mdl5 <- lme4::lmer(
    data = neuron_reading,
    mean_int_rel ~ pos_z_lim * I(pos_z_lim^3) * group + (pos_z_lim | filename) + (pos_z_lim | animal_id),
    REML = TRUE,
    control = lme4::lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
    )
)

# compare models with AIC
clarity_mdl_AIC <- AIC(
    clarity_mdl0,
    clarity_mdl1,
    clarity_mdl2,
    clarity_mdl3,
    clarity_mdl4,
    clarity_mdl5
)
clarity_mdl_AIC

# with f-test
anova(clarity_mdl0, clarity_mdl1)
anova(clarity_mdl1, clarity_mdl2)
anova(clarity_mdl2, clarity_mdl3)
anova(clarity_mdl3, clarity_mdl4)
anova(clarity_mdl4, clarity_mdl5)

# fit best model
clarity_mdl_opt <- lme4::lmer(
    data = neuron_reading,
    mean_int_rel ~ pos_z_lim * I(pos_z_lim^2) * group + (pos_z_lim | filename) + (pos_z_lim | animal_id),
    control = lme4::lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
    )
)

clarity_emm <- emmeans::emmeans(
    clarity_mdl_opt,
    revpairwise ~ group | pos_z_lim * I(pos_z_lim^2),
    at = list(pos_z_lim = seq(0, 1, 0.05))
)
clarity_emm

clarity_anterior_posterior <- emmeans::emmeans(
    clarity_mdl_opt,
    revpairwise ~ pos_z_lim * I(pos_z_lim^2) * group,
    at = list(pos_z_lim = c(0.2715, 0.539)) # 1st and 3rd quartile
)$contrasts %>% broom.mixed::tidy(conf.int = TRUE)
clarity_anterior_posterior


clarity_anterior_posterior_emm <- emmeans::emmeans(
    clarity_mdl_opt,
    revpairwise ~ group * pos_z_lim * I(pos_z_lim^2),
    at = list(pos_z_lim = c(0.2715, 0.539)) # 1st and 3rd quartile
)$contrasts %>% broom.mixed::tidy(conf.int = TRUE)
clarity_anterior_posterior_emm

clarity_emtrend <- emmeans::emtrends(
    clarity_mdl_opt,
    revpairwise ~ group | pos_z_lim * I(pos_z_lim^2),
    var = "pos_z_lim",
    at = list(pos_z_lim = seq(0, 1, 0.01))
)
clarity_emtrend

clarity_emtrend_q <- emmeans::emtrends(
    clarity_mdl_opt,
    revpairwise ~ group | pos_z_lim * I(pos_z_lim^2),
    var = "pos_z_lim",
    at = list(pos_z_lim = summary(neuron_reading$pos_z_lim)[c(2, 3, 5)])
)$contrasts %>% broom.mixed::tidy(conf.int = TRUE)
clarity_emtrend_q

clarity_emtrend_ov <- emmeans::emtrends(
    clarity_mdl_opt,
    revpairwise ~ group | pos_z_lim * I(pos_z_lim^2),
    var = "pos_z_lim"
)
clarity_emtrend_ov

#### counts ----

# this tries to do the same procedure made for intensity but for counts.
# note that the probability of a given neuron to be detected as oxa is going
# to be some function over its intensity, that is, neurons brighter (with the oxa color)
# than background are considered as a count, so controlling for this makes sense

total_counts <- neuron_reading %>%
    group_by(animal_id, group) %>%
    summarise(
        oxa_count = n()
    )

# not necessary but I did not like super long filenames
mean_position <- neuron_reading %>%
    mutate(
        hem = str_extract(filename, pattern = "(izq|der)")
    )

# the model is asking for differences in the mean pos_z_lim, which is the same as saying
# where's is the mass of neurons situated, controlling for mean intensity, which we know
# from the previous model is heterogenous among groups along the antero-posterior axis
pos_mdl <- lme4::lmer(
    data = mean_position,
    pos_z_lim ~ group + mean_int_rel + (1 | animal_id),
    control = lme4::lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
    )
)
summary(pos_mdl)

pos_emm <- emmeans::emmeans(
    pos_mdl,
    revpairwise ~ group + mean_int_rel,
    type = "response"
)
pos_emm

pos_emm_p <- broom.mixed::tidy(pos_emm$emmeans, conf.int = TRUE)
broom.mixed::tidy(pos_emm$contrasts, conf.int = TRUE)


## p::clarity total counts ----

clarityp4 <- total_counts %>%
    ggplot(aes(
        group, oxa_count
    )) +
    geom_boxplot(outlier.shape = NA, aes(color = group), width = 0.5) +
    geom_point(aes(fill = group)) +
    scale_x_discrete(
        labels = c("LU", "HU")
    ) +
    theme_uncertainty +
    scale_fill_manual(values = c("black", "orange")) +
    scale_color_manual(values = c("black", "orange")) +
    scale_y_continuous(
        breaks = seq(1500, 3000, 500),
        limits = c(1500, 3000),
        expand = c(0, 0)
    ) +
    xlab("") +
    ylab(latex2exp::TeX(r"($OXA^{+} counts$)"))
clarityp4

## p::clarity count distribution ----

clarityp5 <- neuron_reading %>%
    ggplot(aes(
        pos_z_lim
    )) +
    geom_density(aes(color = group, fill = group), linewidth = 1, alpha = 0.1) +
    geom_vline(xintercept = 0.356, linetype = "dashed", color = "black") +
    geom_vline(xintercept = 0.481, linetype = "dashed", color = "orange") +
    geom_point(
        data = pos_emm_p,
        inherit.aes = FALSE,
        aes(x = estimate, y = 2.25, fill = group)
    ) +
    geom_errorbarh(
        data = pos_emm_p,
        inherit.aes = FALSE,
        aes(
            y = 2.25, xmin = conf.low, xmax = conf.high,
            color = group
        ),
        height = 0.25,
        linewidth = 1.5
    ) +
    geom_segment(aes(x = 0.356, y = 2.425, xend = 0.481, yend = 2.425),
        inherit.aes = FALSE
    ) +
    theme_uncertainty +
    scale_fill_manual(values = c("black", "orange")) +
    scale_color_manual(values = c("black", "orange")) +
    xlab(latex2exp::TeX(r"($Anterior \rightarrow \ Posterior$)")) +
    ylab(latex2exp::TeX(r"($OXA^{+} counts \ density$)")) +
    scale_y_continuous(
        breaks = seq(0, 2.5, 0.5),
        limits = c(0, 2.5),
        expand = c(0, 0)
    ) +
    scale_x_continuous(
        breaks = seq(0, 1, 0.25),
        limits = c(0, 1),
        expand = c(0, 0)
    )
clarityp5


## p::clarity slopes ----

clarityp1 <- clarity_emtrend$emtrends %>%
    broom.mixed::tidy(conf.int = TRUE) %>%
    ggplot(aes(
        pos_z_lim, pos_z_lim.trend
    )) +
    geom_ribbon(aes(
        ymin = conf.low,
        ymax = conf.high,
        fill = group
    ), alpha = 0.25) +
    geom_line(aes(color = group)) +
    scale_x_continuous(breaks = seq(0, 1, 0.1)) +
    theme_uncertainty +
    scale_y_continuous(
        breaks = seq(-3, 2.5, 0.5),
        limits = c(-3, 2.5),
        expand = c(0, 0)
    ) +
    ylab(latex2exp::TeX(r"($\hat{\beta_}_{oxa \ intensity}$)")) +
    xlab(latex2exp::TeX(r"($Anterior \rightarrow \ Posterior$)")) +
    scale_fill_manual(values = c("black", "orange")) +
    scale_color_manual(values = c("black", "orange"))
clarityp1

clarityp2 <- clarity_emm$emmeans %>%
    broom.mixed::tidy(conf.int = TRUE) %>%
    ggplot(aes(
        pos_z_lim, estimate
    )) +
    geom_point(
        data = sample_n(neuron_reading, size = 100),
        aes(pos_z_lim, mean_int_rel, fill = group),
        alpha = 0.5,
        size = 1.5
    ) +
    geom_ribbon(aes(
        ymin = conf.low,
        ymax = conf.high,
        fill = group
    ), alpha = 0.25) +
    geom_line(aes(color = group)) +
    geom_vline(
        xintercept = summary(neuron_reading$pos_z_lim)[c(2, 5)],
        linetype = "dashed"
    ) +
    scale_x_continuous(breaks = seq(0, 1, 0.1)) +
    theme_uncertainty +
    scale_y_continuous(
        breaks = seq(0, 1.4, 0.2),
        limits = c(0, 1.4),
        expand = c(0, 0)
    ) +
    ylab(latex2exp::TeX(r"($\hat{\mu_}_{oxa \ intensity}$)")) +
    xlab(latex2exp::TeX(r"($Anterior \rightarrow \ Posterior$)")) +
    scale_fill_manual(values = c("black", "orange")) +
    scale_color_manual(values = c("black", "orange"))
clarityp2

clarityp3 <- clarity_anterior_posterior_emm %>%
    slice(c(2, 5)) %>%
    ggplot(aes(
        contrast, estimate
    )) +
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_uncertainty +
    scale_y_continuous(
        breaks = seq(-0.2, 0.4, 0.1),
        limits = c(-0.2, 0.4),
        expand = c(0, 0)
    ) +
    scale_x_discrete(labels = c(
        "LU",
        "HU"
    )) +
    ylab(latex2exp::TeX(r"($\hat{\mu_}_{posterior - anterior}$)")) +
    xlab("") +
    scale_fill_manual(values = c("black", "orange")) +
    scale_color_manual(values = c("black", "orange"))
clarityp3



## p::metabolic efficiency ----
ip1 <- HU_vs_LU_mdl %>%
    filter(experimental_phase == "Experimental") %>%
    ggplot(aes(
        rel_date, resid_hu
    )) +
    stat_summary(
        fun.data = "mean_se",
        geom = "ribbon",
        aes(group = 1),
        fill = "orange"
    ) +
    stat_summary(
        fun.data = "mean_se",
        geom = "line",
        aes(group = 1)
    ) +
    geom_line(aes(group = ID), alpha = 0.25) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    scale_x_continuous(breaks = seq(0, 56, 7)) +
    theme_uncertainty +
    scale_y_continuous(
        breaks = seq(-10, 4, 2),
        limits = c(-10, 4),
        expand = c(0, 0)
    ) +
    ylab(latex2exp::TeX(r"($HU_{\Delta wt./\Delta int.} - LU_{\Delta wt./\Delta int.}$)")) +
    xlab("Days")
ip1




## p::pellet intake over weeks ----
ip2 <- id_wd %>%
    ggplot(aes(
        weeks, intake,
        group = ID
    )) +
    geom_vline(xintercept = 2.5, linetype = "dashed") +
    stat_summary(
        fun.data = "mean_se",
        geom = "point",
        aes(group = interaction(ID, weeks), color = group)
    ) +
    stat_summary(
        fun.data = "mean_se",
        geom = "pointrange",
        aes(group = interaction(group, weeks), color = group)
    ) +
    stat_summary(
        fun.data = "mean_se",
        geom = "line",
        aes(group = group, color = group)
    ) +
    scale_x_continuous(breaks = seq(1, 10, 1)) +
    theme_uncertainty +
    scale_y_continuous(
        breaks = seq(0, 100, 10),
        limits = c(0, 100),
        expand = c(0, 0)
    ) +
    ylab("Weekly # pellets") +
    xlab("Weeks") +
    scale_fill_manual(values = c("black", "orange")) +
    scale_color_manual(values = c("black", "orange"))
ip2

### p::pellet intake @ week 10
ip3 <- id_wd %>%
    ungroup() %>%
    group_by(ID, group) %>%
    summarise(
        intake = mean(intake)
    ) %>%
    ggplot(aes(
        group, intake,
        color = group
    )) +
    geom_boxplot(outlier.shape = NA, width = 0.5, aes(color = group)) +
    geom_point(aes(fill = group), color = "black") +
    ggsignif::geom_signif(
        comparisons = list(c(1, 2)),
        textsize = 0,
        y_position = 87
    ) +
    theme_uncertainty +
    scale_x_discrete(labels = c("LU", "HU")) +
    scale_y_continuous(
        breaks = seq(0, 100, 10),
        limits = c(0, 100),
        expand = c(0, 0)
    ) +
    ylab("# Pellets") +
    xlab("") +
    scale_fill_manual(values = c("black", "orange")) +
    scale_color_manual(values = c("black", "orange"))
ip3

## p::failed retrievals ----
frp1 <- broom.mixed::tidy(frd_emm$emmeans, conf.int = TRUE) %>%
    ggplot(aes(
        group, response
    )) +
    geom_pointrange(size = 1, aes(color = group, ymin = conf.low, ymax = conf.high)) +
    geom_point(
        aes(fill = group),
        color = "black",
        alpha = 0.25,
        data = frd %>%
            group_by(ID, group) %>%
            summarise(response = mean(failed_retrievals))
    ) +
    ggsignif::geom_signif(
        comparisons = list(c(1, 2)),
        textsize = 0,
        y_position = 19
    ) +
    theme_uncertainty +
    scale_x_discrete(labels = c("LU", "HU")) +
    scale_y_continuous(
        breaks = seq(0, 20, 5),
        limits = c(0, 20),
        expand = c(0, 0)
    ) +
    ylab("# Failed retrievals") +
    xlab("") +
    scale_color_manual(values = c("black", "orange")) +
    scale_fill_manual(values = c("black", "orange"))
frp1

## p::retrieval latency ----
broom.mixed::tidy(rld_emm$emmeans, conf.int = TRUE)
## p::failed retrievals ----
rlp1 <- broom.mixed::tidy(rld_emm$emmeans, conf.int = TRUE) %>%
    ggplot(aes(
        group, response
    )) +
    geom_pointrange(size = 1, aes(color = group, ymin = conf.low, ymax = conf.high)) +
    geom_point(
        data = rld %>%
            group_by(ID, group) %>%
            summarise(ret = mean(ret)),
        aes(group, ret, fill = group),
        color = "black", alpha = 0.25
    ) +
    ggsignif::geom_signif(
        comparisons = list(c(1, 2)),
        textsize = 0,
        y_position = 460
    ) +
    theme_uncertainty +
    scale_x_discrete(labels = c("LU", "HU")) +
    scale_y_continuous(
        breaks = seq(0, 500, 100),
        limits = c(0, 500),
        expand = c(0, 0)
    ) +
    ylab("Retrieval latency (sec.)") +
    xlab("") +
    scale_color_manual(values = c("black", "orange")) +
    scale_fill_manual(values = c("black", "orange"))
rlp1

## p::example dc ----
dcp1 <- dcd %>%
    group_by(ID) %>%
    mutate(
        s = as.numeric(row_number() %% 100 == 0)
    ) %>%
    filter(s == 1) %>%
    ggplot(aes(
        x, y,
        group = ID
    )) +
    geom_line(aes(color = group)) +
    stat_summary(
        fun.data = "mean_se",
        geom = "ribbon",
        alpha = 0.25,
        aes(group = group, fill = group)
    ) +
    stat_summary(
        fun.data = "mean_se",
        geom = "line",
        linewidth = 2,
        aes(group = group, color = group)
    ) +
    theme_uncertainty +
    scale_x_continuous(
        breaks = c(5, 10, 20, 40, 80, 120),
        transform = "log"
    ) +
    scale_y_continuous(
        transform = "log",
        limits = c(exp(1)^-0.5, exp(1)^3),
        expand = c(0, 0),
        breaks = scales::trans_breaks("log", function(x) exp(1)^x),
        labels = scales::trans_format("log", scales::math_format(e^.x))
    ) +
    annotation_logticks() +
    ylab("# Rewards") +
    xlab("Cost") +
    scale_color_manual(values = c("black", "orange")) +
    scale_fill_manual(values = c("black", "orange"))
dcp1

## p::demand curve alpha ----
alphaq0p1 <- alphaq0d %>%
    filter(term == "alpha") %>%
    ggplot(aes(
        group, estimate
    )) +
    geom_boxplot(outlier.shape = NA, width = 0.5, aes(color = group)) +
    geom_point(size = 5, shape = 21, aes(fill = group), color = "black", alpha = 0.5) +
    ggsignif::geom_signif(
        comparisons = list(c(1, 2)),
        textsize = 0
    ) +
    scale_y_continuous(
        transform = "log",
        limits = c(exp(1)^-9, exp(1)^-6),
        expand = c(0, 0),
        breaks = scales::trans_breaks("log", function(x) exp(1)^x),
        labels = scales::trans_format("log", scales::math_format(e^.x))
    ) +
    scale_x_discrete(labels = c("LU", "HU")) +
    annotation_logticks() +
    ggpubr::theme_pubr() +
    theme_uncertainty +
    ylab(expression(alpha)) +
    xlab(" ") +
    scale_color_manual(values = c("black", "orange")) +
    scale_fill_manual(values = c("black", "orange"))
alphaq0p1

## p::demand curve q0 ----
alphaq0p2 <- alphaq0d %>%
    filter(term == "q0") %>%
    ggplot(aes(
        group, estimate
    )) +
    geom_boxplot(outlier.shape = NA, width = 0.5, aes(color = group)) +
    geom_point(size = 5, shape = 21, aes(fill = group), color = "black", alpha = 0.5) +
    ggsignif::geom_signif(
        comparisons = list(c(1, 2)),
        textsize = 0,
        y_position = 25
    ) +
    scale_y_continuous(
        breaks = seq(0, 30, 5),
        limits = c(0, 30),
        expand = c(0, 0)
    ) +
    scale_x_discrete(labels = c("LU", "HU")) +
    ggpubr::theme_pubr() +
    theme_uncertainty +
    ylab(latex2exp::TeX(r"($Q_{0}$)")) +
    xlab(" ") +
    scale_color_manual(values = c("black", "orange")) +
    scale_fill_manual(values = c("black", "orange"))
alphaq0p2

## p::demand curve dq0/dalpha ----

alphaq0p4d <- alphaq0d %>%
    select(-c(std.error, p.value)) %>%
    pivot_wider(
        names_from = term,
        values_from = estimate
    )
alphaq0p4d

alpha_mid <- summary(alphaq0p4d$alpha)["Median"]
q0_mid <- summary(alphaq0p4d$q0)["Median"]

alphaq0p4 <- alphaq0p4d %>%
    ggplot(aes(
        alpha, q0
    )) +
    geom_point(size = 5, shape = 21, aes(fill = group), color = "black", alpha = 0.5) +
    geom_vline(
        xintercept = alpha_mid,
        linetype = "dashed"
    ) +
    geom_hline(
        yintercept = q0_mid,
        linetype = "dashed"
    ) +
    scale_y_continuous(
        breaks = seq(0, 30, 5),
        limits = c(0, 30),
        expand = c(0, 0)
    ) +
    scale_x_continuous(
        transform = "log",
        limits = c(exp(1)^-9, exp(1)^-6),
        expand = c(0, 0),
        breaks = scales::trans_breaks("log", function(x) exp(1)^x),
        labels = scales::trans_format("log", scales::math_format(e^.x))
    ) +
    ggpubr::theme_pubr() +
    theme_uncertainty +
    ylab(latex2exp::TeX(r"($Q_{0}$)")) +
    xlab(latex2exp::TeX(r"($\alpha$)")) +
    scale_color_manual(values = c("black", "orange")) +
    scale_fill_manual(values = c("black", "orange"))
alphaq0p4

alphaq0p3 <- broom::tidy(alpha_q0_emtrend$emtrends, conf.int = TRUE) %>%
    ggplot(aes(
        group, q0.trend
    )) +
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high, color = group), size = 1) +
    boxplot_sig_bracket(1, 2) +
    scale_y_continuous(
        breaks = seq(-0.2, 0.5, 0.1),
        limits = c(-0.2, 0.5),
        expand = c(0, 0)
    ) +
    scale_x_discrete(labels = c("LU", "HU")) +
    ggpubr::theme_pubr() +
    theme_uncertainty +
    ylab(latex2exp::TeX(r"($\hat{\mu}_{\frac{\Delta Q_{0}}{\Delta \alpha}}$)")) +
    xlab(" ") +
    scale_color_manual(values = c("black", "orange")) +
    scale_fill_manual(values = c("black", "orange"))
alphaq0p3

## p::rna-seq pdp ----
tmp <- rnaseq_pdp %>%
    filter(symbol %in% rra_res$Name) %>%
    group_by(symbol) %>%
    mutate(
        x = scales::rescale(x, to = c(0, 1)),
        yhat = scales::rescale(yhat, to = c(0, 1)),
        symbol = as.factor(symbol)
    ) %>%
    select(symbol, x, yhat) %>%
    drop_na()
tmp

mdl_tmp <- gam(
    x ~ s(yhat, by = symbol),
    data = tmp
)

pred_grid <- expand_grid(
    yhat = seq(0, 1, 0.25),
    x = seq(0, 1, 0.25),
    symbol = unique(tmp$symbol)
) %>%
    mutate(
        .pred = predict(mdl_tmp, newdata = .)
    )

rnap1 <- pred_grid %>%
    ggplot(aes(
        yhat, symbol,
        fill = .pred
    )) +
    geom_tile() +
    scale_x_continuous(
        breaks = seq(0, 1, 0.25),
        expand = c(0, 0)
    ) +
    scale_y_discrete(
        expand = c(0, 0)
    ) +
    scale_fill_gradient(
        low = "white",
        high = "orange",
        limits = c(0, 1),
        breaks = c(0.1, 0.9),
        labels = c("Low", "High")
    ) +
    theme_uncertainty +
    theme(
        legend.position = "right"
    ) +
    labs(
        fill = "Expression"
    ) +
    coord_fixed(1 / 4) +
    xlab(latex2exp::TeX(r"($\hat{\mu}_{\frac{\Delta Pr(Group = HU)}{\Delta Gene \ expression}}$)")) +
    ylab("")
rnap1

## p::rna-seq pdp-slopes ----
rnap2 <- pdp_slopes %>%
    mutate(gene = fct_reorder(gene, desc(estimate))) %>%
    ggplot(aes(
        gene, estimate
    )) +
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
        color = "black",
        fill = "black",
        alpha = 0.5,
        size = 1,
        shape = 21
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_y_continuous(
        breaks = seq(-0.05, 0.02, 0.01),
        limits = c(-0.05, 0.02),
        expand = c(0, 0)
    ) +
    xlab("") +
    ylab(latex2exp::TeX(r"($\hat{\mu}_{\frac{\Delta Pr(Group = HU)}{\Delta Gene \ expression}}$)")) +
    theme_uncertainty +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
rnap2

## p::rna-seq gene rank ----
rnap3 <- rra_res %>%
    filter(Score < 0.05) %>%
    mutate(Name = fct_reorder(Name, desc(Score))) %>%
    ggplot(aes(
        Name, Score
    )) +
    geom_segment(aes(Name, xend = Name, y = exp(1)^1, yend = Score), color = "grey") +
    geom_point(fill = "black") +
    scale_y_continuous(
        transform = trans_reverser("log"),
        limits = c(exp(1)^1, exp(1)^-20),
        expand = c(0, 0),
        breaks = scales::trans_breaks("log", function(x) exp(1)^x),
        labels = scales::trans_format("log", scales::math_format(e^.x))
    ) +
    theme_uncertainty +
    xlab("") +
    ylab(latex2exp::TeX(r"($\rho^{-1}_{Gene \ rank}$)")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
rnap3


# figures ----

## fig1 ----
exp_setup_fig1 <- fig(
    path = "../figures/fig1_experimental_setup.png",
    aspect.ratio = "free"
) +
    theme(text = element_text(size = 24))

fig1 <-
    exp_setup_fig1 + theme(aspect.ratio = 1) +
    wp3 + theme(aspect.ratio = 1) +
    ip3 + theme(aspect.ratio = 1) +
    ip1 + theme(aspect.ratio = 1) +
    frp1 + theme(aspect.ratio = 1) +
    rlp1 + theme(aspect.ratio = 1) +
    plot_layout(
        ncol = 3,
        nrow = 2,
        widths = 1,
        heights = 1
    ) +
    plot_annotation(tag_levels = c("A"))
fig1
ggsave(
    plot = fig1,
    filename = "../figures/fig1.png",
    width = 5120,
    height = 2880,
    units = "px",
    dpi = 300
)

## fig2 ----

exp_setup_fig2 <- fig(
    path = "../figures/fig2_experimental_setup.png",
    aspect.ratio = "free"
) +
    theme(text = element_text(size = 24))

fig2 <-
    exp_setup_fig2 + theme(aspect.ratio = 1) +
    dcp1 + theme(aspect.ratio = 1) +
    alphaq0p1 + theme(aspect.ratio = 1) +
    alphaq0p2 + theme(aspect.ratio = 1) +
    alphaq0p4 + theme(aspect.ratio = 1) +
    alphaq0p3 + theme(aspect.ratio = 1) +
    plot_layout(
        ncol = 3,
        nrow = 2,
        widths = 1,
        heights = 1
    ) +
    plot_annotation(tag_levels = c("A"))
ggsave(
    plot = fig2,
    filename = "../figures/fig2.png",
    width = 5120,
    height = 2880,
    units = "px",
    dpi = 300
)

## fig3 ----

exp_setup_fig3 <- fig(
    path = "../figures/rna_seq_diag.png",
    aspect.ratio = "free"
) +
    theme(
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        text = element_text(size = 24)
    )

fig3 <-
    (exp_setup_fig3 + theme(aspect.ratio = 2)) +
    (rnap3 + theme(aspect.ratio = 1)) /
        (rnap2 + theme(aspect.ratio = 1)) +
    (rnap1 + theme(aspect.ratio = 1)) +
    plot_layout(
        ncol = 3,
        nrow = 1,
        widths = c(1, 1, 1, 1),
        heights = 1
    ) +
    plot_annotation(tag_levels = c("A"))
fig3
ggsave(
    plot = fig3,
    filename = "../figures/fig3.png",
    width = 5120,
    height = 2880,
    units = "px",
    dpi = 300
)

## fig4 ----

exp_setup_fig4 <- fig(
    path = "../figures/clarity_example_annotated.PNG",
    aspect.ratio = "free"
) +
    theme(text = element_text(size = 24))

fig4 <-
    exp_setup_fig4 +
    clarityp2 + theme(aspect.ratio = 1) +
    clarityp1 + theme(aspect.ratio = 1) +
    clarityp3 + theme(aspect.ratio = 1) +
    clarityp4 + theme(aspect.ratio = 1) +
    clarityp5 + theme(aspect.ratio = 1) +
    plot_layout(
        ncol = 3,
        nrow = 2,
        widths = 1,
        heights = 1
    ) +
    plot_annotation(tag_levels = c("A"))
ggsave(
    plot = fig4,
    filename = "../figures/fig4.png",
    width = 5120,
    height = 2880,
    units = "px",
    dpi = 300
)
