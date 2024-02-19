# Read and wrangle MH's Mplus model parameters for comparison

library(MplusAutomation)

mplus <- readModels("data/mh/Output Files")
mplus <- map(mplus, ~pluck(., "parameters") |> pluck("unstandardized"))
names(mplus) <- str_extract(names(mplus), "equation\\.[0-9]")

mplus$equation.4 <- mplus$equation.4[c(3, 5, 6, 4, 2, 1, 7, 9, 10, 8), ]
mplus$equation.5 <- mplus$equation.5[c(11, 7, 8, 12, 3, 4, 13, 5, 6, 10, 9, 2, 1, 15, 16, 17, 14), ]
mplus$equation.6 # This is missing parameters
mplus$equation.7

mplus <- map(mplus, ~rename(., q50 = est, q2.5 = lower_2.5ci, q97.5 = upper_2.5ci))

mplus <- map(
  mplus, 
  ~mutate(
    .,
    across(where(is.numeric), ~number(., .01)),
    mplus = str_glue("{q50} [{q2.5}, {q97.5}]")
  )
)
