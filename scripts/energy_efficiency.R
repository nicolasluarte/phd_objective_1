# libs ----

pacman::p_load(
    tidyverse,
    ggplot2
)

# Set path to script location ----

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

# data
intake_data <- read_csv("../datasets/intake/intake_uncertainty_a_fix.csv") %>% 
    group_by(animal, protocol, date) %>% 
    summarise(
        daily_intake = n()
    )
weight_data <- read_csv('../datasets/weights/objective_1_weights.csv') %>% 
    select(date, animal, weight)
metadata <- read_csv("../datasets/metadata/general_metadata.csv") %>% 
    rename(animal = ID)

intake_weight_data <- left_join(intake_data, weight_data) %>% 
    left_join(., metadata, by = "animal")

# get relative date
rel_intake_weight_data <- intake_weight_data %>% 
    ungroup() %>% 
    group_by(animal) %>% 
    mutate(
        rel_dat = as.numeric(date - min(date)),
        interp_weight = zoo::na.approx(weight, rel_dat, na.rm = FALSE),
        cum_intake = cumsum(daily_intake),
        rel_weight = interp_weight - head(interp_weight, n=1)
    )
rel_intake_weight_data

rel_intake_weight_data %>% 
    ggplot(aes(
        cum_intake, rel_weight, color = group
    )) +
    geom_point()

LU_mdl <- lme4::lmer(
    data = rel_intake_weight_data %>% filter(group=="control"),
    rel_weight ~ cum_intake + (1|animal)
)

performance::icc(LU_mdl)

HU_data <- rel_intake_weight_data %>% 
    filter(
        group == "treatment"
    )

HU_vs_LU_mdl <- HU_data %>% 
    ungroup() %>% 
    mutate(
        preds = predict(LU_mdl, newdata = HU_data, re.form = NA),
        resid_hu = rel_weight - preds
    )
HU_vs_LU_mdl


resid_mdl <- lmerTest::lmer(
    data = HU_vs_LU_mdl,
    resid_hu ~ rel_dat + (1|animal)
)
summary(resid_mdl)
broom.mixed::tidy(resid_mdl, conf.int = TRUE)
emmeans::emmeans(resid_mdl, ~rel_dat, at=list(rel_dat=56))





