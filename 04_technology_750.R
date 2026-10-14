# Run from the Four_Real_Data_Analyses project root.
# Default: fit both models; completed fits and completed EM iterations resume.
# To reconstruct saved historical results instead:
# Sys.setenv(REAL_DATA_MODE = "report")
source("R/real_data.R")
run_real_analysis("technology_750", mode = Sys.getenv("REAL_DATA_MODE", "fit"))
