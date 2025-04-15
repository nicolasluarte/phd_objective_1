# lib load ----
pacman::p_load(
    tidyverse,
    ggplot2,
    patchwork
)
setwd(this.path::here())

# Plot stuff ----
theme_uncertainty <- ggpubr::theme_pubr() +
    update_geom_defaults("point", list(size = 5, alpha = 0.5, shape = 21)) +
    theme(
        text = element_text(size = 24),
        axis.text=element_text(size=14),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        legend.position = "none"
    )

boxplot_sig_bracket <- function(group1, group2){
    ggsignif::geom_signif(
        comparisons = list(c(group1, group2)),
        map_signif_level = TRUE,
        textsize = 0,
        tip_length = 0
    )
}


# Weight (w) ----
wd <- read_csv("../datasets/weights/weights.csv") %>% 
    group_by(ID) %>% 
    mutate(
        rel_date = as.numeric(date - min(date))
    )

# same but weekly aggregate
wd_weekly <- wd %>% 
    ungroup() %>% 
    group_by(ID) %>% 
    mutate(
        rel_date = cumsum(replace_na(as.numeric(date - lag(date)),0)),
        days_from_end = as.numeric(abs(date - max(date))),
        phase = if_else(days_from_end <= 56, "uncertainty", "baseline"),
        week = as.numeric(as.factor(lubridate::isoweek(date))),
        week = if_else(phase == "uncertainty", week+1, week)
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
    data = weekly_weight %>% filter(phase=="uncertainty"),
    rel_weight ~ group * week + (1|ID)
)

w_mdl1_emm <- emmeans::emmeans(
    w_mdl1,
    pairwise ~ group|week,
    type = "response",
    at = list(week=c(10))
)
w_mdl1_emm
w_mdl1_emm$contrasts %>% broom.mixed::tidy(., conf.int=TRUE)

### weight trends ----
w_mdl1_trend <- emmeans::emtrends(
    w_mdl1,
    pairwise ~ group*week,
    var = "week"
)
w_mdl1_trend
w_mdl1_trend$emtrends %>% broom.mixed::tidy(., conf.int=TRUE)



## p:weekly weight ----
wp1 <- weekly_weight %>% 
    filter(week <= 10) %>% 
    ggplot(aes(
        week, rel_weight, group = ID
    )) +
    geom_vline(xintercept = 2.5, alpha = 0.5) +
    geom_line(alpha = 0.25, aes(color = group))  +
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
    scale_y_continuous(breaks = seq(-10,10, 2), 
                       limits = c(-10,10), 
                       expand = c(0,0)) +
    ylab(expression(Delta*" Weekly weight (gr.)")) +
    xlab("Weeks") +
    scale_fill_manual(values = c("black", "orange")) +
    scale_color_manual(values = c("black", "orange")) 
wp1

## p:weight @ week10 ----
wp2 <- weekly_weight %>% 
    filter(week==10) %>% 
    ggplot(aes(
        group, weight, color=group
    )) +
    geom_boxplot(outlier.shape = NA, width=0.5, aes(color=group)) +
    geom_point(aes(fill=group), color="black") +
    boxplot_sig_bracket(1,2) +
    theme_uncertainty + 
    scale_x_discrete(labels = c("Ctrl", "Unc")) +
    scale_y_continuous(breaks = seq(20,35,5), 
                       limits = c(20,35), 
                       expand = c(0,0)) +
    ylab("Final weight (gr.)") +
    xlab("Weeks") +
    scale_fill_manual(values = c("black", "orange")) +
    scale_color_manual(values = c("black", "orange")) 
wp2

## p:weight trends
wp3 <- w_mdl1_trend$emtrends %>% 
    broom.mixed::tidy(., conf.int=TRUE) %>% 
    ggplot(aes(
        group, week.trend, ymin=conf.low, ymax=conf.high
    )) +
    geom_pointrange(aes(color=group), size=1) +
    geom_hline(yintercept = 0, color="gray") +
    boxplot_sig_bracket(1,2) +
    scale_fill_manual(values = c("black", "orange")) +
    scale_color_manual(values = c("black", "orange")) +
    theme_uncertainty + 
    scale_x_discrete(labels = c("Ctrl", "Unc")) +
    scale_y_continuous(breaks = seq(0,0.6,0.1), 
                       limits = c(0,0.6), 
                       expand = c(0,0)) +
    ylab(latex2exp::TeX(r"($\beta_{weight/weeks}$)")) +
    xlab("")
wp3

# Intake (I) ----
id <- read_csv("../datasets/intake/intake.csv")

# data for metabolic efficiency, weight interpolated
id_wd <- left_join(id %>% ungroup, wd %>% ungroup) %>% 
    filter(experimental_phase == "Experimental") %>% 
    ungroup() %>% 
    group_by(ID) %>% 
    arrange(date, .by_group = TRUE) %>% 
    mutate(
        rel_date = as.numeric(date - min(date)),
        weeks = cumsum(as.numeric(rel_date%%7==0)),
        interp_weight = zoo::na.approx(weight, rel_date, na.rm = FALSE),
        cum_intake = cumsum(intake)-head(intake,n=1)
    ) %>% 
    drop_na(interp_weight) %>% 
    mutate(
        delta_weight = interp_weight - head(interp_weight, n=1)
    )
id_wd

## mdls ----

### pellet intake ----
pellet_mdl <- lmerTest::lmer(
    data = id_wd,
    intake ~ rel_date * group + (rel_date|ID),
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
broom.mixed::tidy(pellet_emm$emmeans, conf.int=TRUE)
broom.mixed::tidy(pellet_emm$contrasts, conf.int=TRUE)

### energy efficiency ----
LU_mdl <- lmerTest::lmer(
    data = id_wd %>% filter(group=="No-Uncertainty"),
    delta_weight ~ cum_intake * rel_date + (rel_date|ID),
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
        preds = predict(LU_mdl, newdata = HU_data, re.form=NA),
        resid_hu = delta_weight - preds
    )
HU_vs_LU_mdl

# model analysis
resid_mdl <- lmerTest::lmer(
    data = HU_vs_LU_mdl %>% filter(experimental_phase=="Experimental"),
    resid_hu ~ rel_date + (rel_date|ID),
    control = lme4::lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
    )
)
summary(resid_mdl)
broom.mixed::tidy(resid_mdl, conf.int=TRUE)
emmeans::emmeans(resid_mdl, ~rel_date, at=list(rel_date=c(56)),
                 type = "response") %>% 
    broom.mixed::tidy(conf.int=TRUE)

### failed retrievals ----
frd <- read_csv("../datasets/behavioral-tests/boops.csv") %>% 
    group_by(ID) %>% 
    mutate(
        rel_date = as.numeric(date - min(date))
    )

frd_mdl <- lme4::glmer.nb(
    data = frd,
    failed_retrievals ~ group * rel_date * experimental_phase + (rel_date|ID),
    control = lme4::glmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
    )
)
summary(frd_mdl)

frd_emm <- emmeans::emmeans(
    frd_mdl,
    pairwise ~ group | rel_date * experimental_phase,
    at = list(rel_date=c(56), experimental_phase=c("Experimental")),
    type = "response"
)
frd_emm
broom.mixed::tidy(frd_emm$emmeans, conf.int=TRUE)
broom.mixed::tidy(frd_emm$contrasts, conf.int=TRUE)

frd_emm_trend <- emmeans::emtrends(
    frd_mdl,
    pairwise ~ group | rel_date * experimental_phase,
    at = list(experimental_phase=c("Experimental"), rel_date = c(28)),
    var = "rel_date",
    type = "response"
)
frd_emm_trend
broom.mixed::tidy(frd_emm_trend$emtrends, conf.int=TRUE)
broom.mixed::tidy(frd_emm_trend$contrasts, conf.int=TRUE)

### retrieval latency ----
rld <- read_csv("../datasets/behavioral-tests/retrieval_latency.csv") %>% 
    filter(ret>0, ret<=1378) %>% 
    group_by(ID) %>% 
    mutate(
        rel_date = as.numeric(date - min(date))
    )

rld_mdl <- lme4::glmer(
    data=rld, 
    ret ~ group * experimental_phase * rel_date + (1|ID),
    family = Gamma(link="log"),
    control = lme4::glmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
    )
)
summary(rld_mdl)

rld_emm <- emmeans::emmeans(
    rld_mdl,
    pairwise ~ group | experimental_phase * rel_date,
    at = list(experimental_phase=c("Experimental"), rel_date=c(56)),
    type = "response"
)
rld_emm
broom.mixed::tidy(rld_emm$emmeans, conf.int=TRUE)
broom.mixed::tidy(emmeans::contrast(rld_emm, "revpairwise"), conf.int=TRUE)

rld_emm_trend <- emmeans::emtrends(
    rld_mdl,
    pairwise ~ group | experimental_phase * rel_date,
    at = list(experimental_phase=c("Experimental"), rel_date = c(30)),
    var = "rel_date",
    type = "response"
)
rld_emm_trend
broom.mixed::tidy(rld_emm_trend$emtrends, conf.int=TRUE)
broom.mixed::tidy(rld_emm_trend$contrasts, conf.int=TRUE)

### demand curve ----
alphaq0d <- read_csv("../datasets/behavioral-tests/alpha_q0.csv")
dcd <- read_csv("../datasets/behavioral-tests/demand_curve_data.csv")

alpha_mdl <- glm(
    data=alphaq0d %>% filter(term=="alpha"),
    estimate ~ group,
    family = Gamma(link="log")
)
summary(alpha_mdl)

alpha_emm <- emmeans::emmeans(
    alpha_mdl,
    pairwise ~ group,
    type = "response"
)
alpha_emm
broom.mixed::tidy(alpha_emm$emmeans, conf.int=TRUE)
broom.mixed::tidy(alpha_emm$contrasts, conf.int=TRUE)

q0_mdl <- glm(
    data=alphaq0d %>% filter(term=="q0"),
    estimate ~ group,
    family = Gamma(link="log")
)
summary(q0_mdl)

q0_emm <- emmeans::emmeans(
    q0_mdl,
    pairwise ~ group,
    type = "response"
)
q0_emm
broom.mixed::tidy(q0_emm$emmeans, conf.int=TRUE)
broom.mixed::tidy(q0_emm$contrasts, conf.int=TRUE)

# experimental this is the correlation between alpha and q0
alpha_q0_corr <- glm(
    data = alphaq0d %>%
        ungroup() %>% 
        select(ID, group, term, estimate) %>% 
        pivot_wider(values_from = "estimate", names_from = "term"),
    alpha ~ q0 * group,
    family = Gamma(link="log")
)
summary(alpha_q0_corr)


alpha_q0_emtrend <- emmeans::emtrends(
    alpha_q0_corr,
    pairwise ~ group | q0,
    var = "q0",
    type = "response"
)
alpha_q0_emtrend
broom::tidy(alpha_q0_emtrend$emtrends, conf.int=TRUE)
broom::tidy(alpha_q0_emtrend$contrasts, conf.int=TRUE)

## RNA-seq ----
rnaseq_rank <- read_csv("../datasets/rna-seq/data_rank.csv")
rnaseq_pdp <- read_csv("../datasets/rna-seq/data_rank_pdp.csv")

# generate lists of genes
rank_list <- rnaseq_rank %>% 
    select(symbol, rank, iteration) %>% 
    group_by(iteration) %>% 
    arrange(rank, .by_group = TRUE) %>% 
    group_split() %>% 
    map(., function(it){
        return(c(it$symbol))
    }) %>% 
    set_names(paste0("iteration", unique(rnaseq_rank$iteration)))
rank_list

rra_mdl <- RobustRankAggreg::aggregateRanks(
    glist = rank_list,
    N = 26
)
rra_mdl

rra_res <- tibble(rra_mdl) %>% 
    filter(Score!=1)

# get slopes from pdp
pdp_slopes <- rnaseq_pdp %>% 
    filter(symbol%in%rra_res$Name) %>% 
    group_by(symbol) %>% 
    group_split() %>% 
    map_dfr(., function(X){
        lm(data=X, yhat~x+I(x^3)) %>% 
            broom::tidy(., conf.int=TRUE) %>% 
            filter(term=="x") %>% 
            mutate(gene=X$symbol[1])
    })
pdp_slopes


## p::metabolic efficiency ----
ip1 <- HU_vs_LU_mdl %>% 
    filter(experimental_phase=="Experimental") %>% 
    ggplot(aes(
        rel_date, resid_hu
    )) +
    geom_line(aes(group=ID), alpha=0.5) +
    geom_point(color="black", fill="orange") +
    geom_hline(yintercept = 0, color="gray") +
    scale_x_continuous(breaks = seq(0, 56, 7)) +
    theme_uncertainty + 
    scale_y_continuous(breaks = seq(-10,2, 2), 
                       limits = c(-10,2), 
                       expand = c(0,0)) +
    ylab(latex2exp::TeX(r"($HU_{\Delta wt./\Delta int.} - LU_{\Delta wt./\Delta int.}$)")) +
    xlab("Days")
ip1




## p::pellet intake over weeks ----
ip2 <-id_wd %>% 
    ggplot(aes(
        weeks, intake, group = ID
    )) +
    geom_vline(xintercept = 2.5, alpha = 0.5) +
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
    scale_y_continuous(breaks = seq(0,100, 10), 
                       limits = c(0,100), 
                       expand = c(0,0)) +
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
        group, intake, color=group
    )) + 
    geom_boxplot(outlier.shape = NA, width=0.5, aes(color=group)) +
    geom_point(aes(fill=group), color="black") +
    boxplot_sig_bracket(1,2) +
    theme_uncertainty + 
    scale_x_discrete(labels = c("Ctrl", "Unc")) +
    scale_y_continuous(breaks = seq(0,100,10), 
                       limits = c(0,100), 
                       expand = c(0,0)) +
    ylab("# Pellets") +
    xlab("") +
    scale_fill_manual(values = c("black", "orange")) +
    scale_color_manual(values = c("black", "orange")) 
ip3

## p::failed retrievals ----
frp1 <- broom.mixed::tidy(frd_emm$emmeans, conf.int=TRUE) %>%
    ggplot(aes(
        group, response
    )) +
    geom_pointrange(size=1, aes(color=group, ymin=conf.low, ymax=conf.high)) +
    geom_point(
        aes(fill=group),
        color="black",
        alpha=0.25,
        data=frd %>% 
                     group_by(ID, group) %>% 
                     summarise(response=mean(failed_retrievals))) +
    boxplot_sig_bracket(1,2) +
    theme_uncertainty + 
    scale_x_discrete(labels = c("Ctrl", "Unc")) +
    scale_y_continuous(breaks = seq(0,20,5), 
                       limits = c(0,20), 
                       expand = c(0,0)) +
    ylab("# Failed retrievals") +
    xlab("") +
    scale_color_manual(values = c("black", "orange")) +
    scale_fill_manual(values = c("black", "orange")) 
frp1

## p::retrieval latency ----
broom.mixed::tidy(rld_emm$emmeans, conf.int=TRUE)
## p::failed retrievals ----
rlp1 <- broom.mixed::tidy(rld_emm$emmeans, conf.int=TRUE) %>%
    ggplot(aes(
        group, response
    )) +
    geom_pointrange(size=1, aes(color=group, ymin=conf.low, ymax=conf.high)) +
    geom_point(data=rld %>%
                   group_by(ID,group) %>% 
                   summarise(ret=mean(ret)),
               aes(group, ret, fill=group),
                 color="black",alpha=0.25) +
    boxplot_sig_bracket(1,2)+
    theme_uncertainty + 
    scale_x_discrete(labels = c("Ctrl", "Unc")) +
    scale_y_continuous(breaks = seq(0,500, 100), 
                       limits = c(0,500), 
                       expand = c(0,0)) +
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
    filter(s==1) %>% 
    ggplot(aes(
        x, y, group = ID
    )) +
    geom_line(aes(color=group)) +
    stat_summary(
        fun.data = "mean_se", 
        geom = "ribbon",
        alpha=0.25,
        aes(group=group, fill=group)
    ) +
    stat_summary(
        fun.data = "mean_se", 
        geom = "line",
        linewidth = 2,
        aes(group=group, color=group)
    ) +
    theme_uncertainty + 
    scale_x_continuous(breaks = c(5,10,20,40,80,120),
                       transform = "log") +
    scale_y_continuous(transform = "log",
                       limits = c(exp(1)^0.1, exp(1)^3),
                       expand = c(0,0),
                       breaks = scales::trans_breaks("log", function(x) exp(1)^x),
                       labels = scales::trans_format("log", scales::math_format(e^.x))) +
    annotation_logticks() +
    ylab("# Rewards") +
    xlab("Cost") +
    scale_color_manual(values = c("black", "orange")) +
    scale_fill_manual(values = c("black", "orange"))
dcp1

##p::demand curve alpha ----
alphaq0p1 <- alphaq0d %>% 
    filter(term=="alpha") %>% 
    ggplot(aes(
        group, estimate
    )) +
    geom_boxplot(outlier.shape = NA, width = 0.5, aes(color=group)) +
    geom_point(size = 5, shape = 21, aes(fill=group), color = "black", alpha=0.5) +
    boxplot_sig_bracket(1, 2) +
    scale_y_continuous(transform = "log",
                       limits = c(exp(1)^-9, exp(1)^-6),
                       expand = c(0,0),
                       breaks = scales::trans_breaks("log", function(x) exp(1)^x),
                       labels = scales::trans_format("log", scales::math_format(e^.x))) +
    scale_x_discrete(labels = c("Ctrl", "Unc")) +
    annotation_logticks() +
    ggpubr::theme_pubr() +
    theme_uncertainty +
    ylab(expression(alpha)) +
    xlab(" ") +
    scale_color_manual(values = c("black", "orange")) +
    scale_fill_manual(values = c("black", "orange"))
alphaq0p1

##p::demand curve q0 ----
alphaq0p2 <- alphaq0d %>% 
    filter(term=="q0") %>% 
    ggplot(aes(
        group, estimate
    )) +
    geom_boxplot(outlier.shape = NA, width = 0.5, aes(color=group)) +
    geom_point(size = 5, shape = 21, aes(fill=group), color = "black", alpha=0.5) +
    boxplot_sig_bracket(1, 2) +
    scale_y_continuous(breaks = seq(0,30, 5), 
                       limits = c(0,30), 
                       expand = c(0,0)) +
    scale_x_discrete(labels = c("Ctrl", "Unc")) +
    ggpubr::theme_pubr() +
    theme_uncertainty +
    ylab(latex2exp::TeX(r"($Q_{0}$)")) +
    xlab(" ") +
    scale_color_manual(values = c("black", "orange")) +
    scale_fill_manual(values = c("black", "orange"))
alphaq0p2

##p::demand curve dq0/dalpha
alphaq0p3 <- broom::tidy(alpha_q0_emtrend$emtrends, conf.int=TRUE) %>% 
    ggplot(aes(
        group, q0.trend
    )) +
    geom_pointrange(aes(ymin=conf.low,ymax=conf.high, color=group), size=1) +
    boxplot_sig_bracket(1, 2) +
    scale_y_continuous(breaks = seq(-0.2,0.5,0.1), 
                       limits = c(-0.2,0.5), 
                       expand = c(0,0)) +
    scale_x_discrete(labels = c("Ctrl", "Unc")) +
    ggpubr::theme_pubr() +
    theme_uncertainty +
    ylab(latex2exp::TeX(r"($\beta_{\frac{\Delta Q_{0}}{\Delta \alpha}}$)")) +
    xlab(" ") +
    scale_color_manual(values = c("black", "orange")) +
    scale_fill_manual(values = c("black", "orange"))
alphaq0p3

## p::rna-seq pdp ----
rnap1 <- rnaseq_pdp %>% 
    filter(symbol %in% rra_res$Name) %>% 
    ggplot(aes(
        x, yhat
    )) +
    geom_smooth(
        method = "lm",
        formula = y ~ poly(x, 3, raw = TRUE),
        color = "black"
    ) +
    theme_uncertainty +
    theme(strip.background = element_rect(fill="white")) +
    ylab(latex2exp::TeX(r"($Pr(Group = HU)$)")) +
    xlab("Gene expression") +
    facet_wrap(~symbol, scale="free") 
rnap1

## p::rna-seq pdp-slopes ----
rnap2 <- pdp_slopes %>% 
    mutate(gene = fct_reorder(gene, desc(estimate))) %>% 
    ggplot(aes(
        gene, estimate
    )) +
    geom_pointrange(aes(ymin=conf.low, ymax=conf.high), color="black", shape=21) +
    geom_hline(yintercept = 0, color = "gray") +
    scale_y_continuous(breaks = seq(-0.05,0.02,0.01), 
                       limits = c(-0.05,0.02), 
                       expand = c(0,0)) +
    xlab("") +
    ylab(latex2exp::TeX(r"($\beta_{\frac{\Delta Pr(Group = HU)}{\Delta Gene \ expression}}$)")) +
    theme_uncertainty
rnap2

## p::rna-seq gene rank ----
rnap3 <- rra_res %>% 
    filter(Score<0.05) %>% 
    mutate(Name = fct_reorder(Name, desc(Score))) %>% 
    ggplot(aes(
        Name, Score
    )) +
    geom_col(fill=NA, color="black") +
    scale_y_continuous(transform = "log",
                       limits = c(exp(1)^-17, exp(1)^1),
                       expand = c(0,0),
                       breaks = scales::trans_breaks("log", function(x) exp(1)^x),
                       labels = scales::trans_format("log", scales::math_format(e^.x))) +
    theme_uncertainty +
    xlab("") +
    ylab(latex2exp::TeX(r"($\rho_{Gene \ rank}$)"))
rnap3


# figures ----

(wp1 + wp2 + wp3 + ip2 + ip3 + ip1) +
    plot_annotation(tag_levels = c("A"))

(frp1 + rlp1 + dcp1 + alphaq0p1 + alphaq0p2 + alphaq0p3) +
    plot_annotation(tag_levels = c("A"))

(rnap3 + rnap2 + rnap1) +
    plot_annotation(tag_levels = c("A"))

