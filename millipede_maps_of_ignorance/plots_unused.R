#############################
## Plots for millipedes SDMs - unused plots
## 
## This script is meant to be sourced from within 'main_millipede_maps.R'.
##
## author: Willson Gaul willson.gaul@ucdconnect.ie
## created: 11 June 2020
## last modified: 11 June 2020
##############################
t_size <- 12

### plot results using 100km block CV -----------------------------------------
ggplot(data = evals[evals$metric == "AUC" & 
                      as.character(evals$block_cv_range) == "1e+05", ], 
       aes(x = factor(train_data, 
                      levels = c("raw", "spat_subsamp"), 
                      labels =  c("raw", "spatially\nundersampled")), 
           y = value, 
           color = factor(
             model, 
             levels = c("day_ll_rf", "env_ll_rf", "spat_ll_rf", 
                        "env_spat_ll_rf"), 
             labels = c("\nDay of Year\n(DOY) +\nList Length\n", 
                        "\nEnvironment +\nDOY +\nList Length\n", 
                        "\nLat + Lon +\nDOY +\nList Length\n", 
                        "\nEnvironment + \nLat + Long +\nDOY +\nList Length")))) + 
  geom_boxplot() + 
  facet_wrap(~species + factor(
    test_data, 
    levels = c("raw", "spat_subsamp"), 
    labels = c("test data - raw", "test data -\nspatially undersampled"))) + 
  xlab("Training Data") + 
  ylab("AUC\n(block Cross-Validated)") + 
  ggtitle("100km block CV") + 
  scale_color_viridis_d(name = "Model", #option = "magma", 
                        begin = 0.1, end = 0.8) + 
  theme_bw() + 
  theme(text = element_text(size = t_size))

ggplot(data = evals[evals$metric == "Kappa" & 
                      as.character(evals$block_cv_range) == "1e+05", ], 
       aes(x = factor(train_data, 
                      levels = c("raw", "spat_subsamp"), 
                      labels =  c("raw", "spatially\nundersampled")), 
           y = value, 
           color = factor(
             model, 
             levels = c("day_ll_rf", "env_ll_rf", "spat_ll_rf", 
                        "env_spat_ll_rf"), 
             labels = c("\nDay of Year (DOY) +\nList Length\n", 
                        "\nEnvironment + DOY +\nList Length\n", 
                        "\nLat + Lon + DOY +\nList Length\n", 
                        "\nEnvironment + \nLat + Long +\nDOY +\nList Length")))) + 
  geom_boxplot() + 
  facet_wrap(~species + factor(
    test_data, 
    levels = c("raw", "spat_subsamp"), 
    labels = c("test data - raw", "test data -\nspatially undersampled"))) + 
  xlab("Training Data") + 
  ylab("Kappa\n(block Cross-Validated)") + 
  ggtitle("100km block CV") + 
  scale_color_viridis_d(name = "Model", #option = "magma", 
                        begin = 0.1, end = 0.8) + 
  theme_bw() + 
  theme(text = element_text(size = t_size))

# plot sensitivity
ggplot(data = evals[evals$metric == "sensitivity" & 
                      as.character(evals$block_cv_range) == "1e+05", ], 
       aes(x = factor(train_data, 
                      levels = c("raw", "spat_subsamp"), 
                      labels =  c("raw", "spatially\nundersampled")), 
           y = value, 
           color = factor(
             model, 
             levels = c("day_ll_rf", "env_ll_rf", "spat_ll_rf", 
                        "env_spat_ll_rf"), 
             labels = c("\nDay of Year (DOY) +\nList Length\n", 
                        "\nEnvironment + DOY +\nList Length\n", 
                        "\nLat + Lon + DOY +\nList Length\n", 
                        "\nEnvironment + \nLat + Long +\nDOY +\nList Length")))) + 
  geom_boxplot() + 
  facet_wrap(~species + factor(
    test_data, 
    levels = c("raw", "spat_subsamp"), 
    labels = c("test data - raw", "test data -\nspatially undersampled"))) + 
  xlab("Training Data") + 
  ylab("Sensitivity\n(block Cross-Validated)") + 
  ggtitle("100km block CV") + 
  scale_color_viridis_d(name = "Model", #option = "magma", 
                        begin = 0.1, end = 0.8) + 
  theme_bw() + 
  theme(text = element_text(size = t_size))


# plot specificity
ggplot(data = evals[evals$metric == "specificity" & 
                      as.character(evals$block_cv_range) == "1e+05", ], 
       aes(x = factor(train_data, 
                      levels = c("raw", "spat_subsamp"), 
                      labels =  c("raw", "spatially\nundersampled")), 
           y = value, 
           color = factor(
             model, 
             levels = c("day_ll_rf", "env_ll_rf", "spat_ll_rf", 
                        "env_spat_ll_rf"), 
             labels = c("\nDay of Year (DOY) +\nList Length\n", 
                        "\nEnvironment + DOY +\nList Length\n", 
                        "\nLat + Lon + DOY +\nList Length\n", 
                        "\nEnvironment + \nLat + Long +\nDOY +\nList Length")))) + 
  geom_boxplot() + 
  facet_wrap(~species + factor(
    test_data, 
    levels = c("raw", "spat_subsamp"), 
    labels = c("test data - raw", "test data -\nspatially undersampled"))) + 
  xlab("Training Data") + 
  ylab("Specificity\n(block Cross-Validated)") + 
  ggtitle("100km block CV") + 
  scale_color_viridis_d(name = "Model", #option = "magma", 
                        begin = 0.1, end = 0.8) + 
  theme_bw() + 
  theme(text = element_text(size = t_size))

# plot specificity
ggplot(data = evals[evals$metric == "Brier" & 
                      as.character(evals$block_cv_range) == "1e+05", ], 
       aes(x = factor(train_data, 
                      levels = c("raw", "spat_subsamp"), 
                      labels =  c("raw", "spatially\nundersampled")), 
           y = value, 
           color = factor(
             model, 
             levels = c("day_ll_rf", "env_ll_rf", "spat_ll_rf", 
                        "env_spat_ll_rf"), 
             labels = c("\nDay of Year (DOY) +\nList Length\n", 
                        "\nEnvironment + DOY +\nList Length\n", 
                        "\nLat + Lon + DOY +\nList Length\n", 
                        "\nEnvironment + \nLat + Long +\nDOY +\nList Length")))) + 
  geom_boxplot() + 
  facet_wrap(~species + factor(
    test_data, 
    levels = c("raw", "spat_subsamp"), 
    labels = c("test data - raw", "test data -\nspatially undersampled"))) + 
  xlab("Training Data") + 
  ylab("Brier score\n(block Cross-Validated)") + 
  ggtitle("100km block CV") + 
  scale_color_viridis_d(name = "Model", #option = "magma", 
                        begin = 0.1, end = 0.8) + 
  theme_bw() + 
  theme(text = element_text(size = t_size))
### end 100km block CV ---------------------------------------------------------
