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

# data ----

data <- read_csv("../datasets/intake/intake_uncertainty_a_fix.csv")

# interval between pellet intake

pellet_intake_intervals <- data %>% 
    group_by(animal, date) %>% 
    mutate(
        interval = (time - lag(time)) - delay,
        t = row_number()
    ) %>% 
    filter(
        interval < 3600/2, interval > 0
    )

pellet_intake_intervals %>% 
    ggplot(aes(
        t, interval, group = animal
    )) +
    geom_point() +
    geom_smooth(aes(group = animal), method = "gam")

pp <- pellet_intake_intervals %>% 
    group_by(animal, date) %>% 
    group_split() %>% 
    map_dfr(
        ., function(X){
            enframe(acf(X$interval, plot = FALSE)$acf) %>% 
                rename(
                    lag = name,
                    ACF = value
                ) %>% 
                mutate(
                    animal = unique(X$animal),
                    date = unique(X$date)
                )
        }
    )
pp

pp %>% 
    ggplot(aes(
        lag, ACF, group = lag
    )) +
    geom_hline(yintercept = 0) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(shape = 21) +
    ggpubr::theme_classic2()

pellet_intake_intervals %>% 
    ggplot(aes(
        t, interval, group = animal
    )) +
    geom_line()


