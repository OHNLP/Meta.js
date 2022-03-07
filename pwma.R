library('meta')

# load all data
all_data <- read.csv('sample.csv')


######################################################################
# Primary Analysis
######################################################################

######################################################################
# Single round for just 1,000 outcomes
######################################################################
# set the n for number of tests
n = 1000
start_time <- Sys.time()
for (i in 1:n) {
  # get the dataf.rame
  ocn = sprintf('M%04d', i)
  df <- all_data[all_data$outcome == ocn, c('study', 'Et', 'Nt', 'Ec', 'Nc')]
  # run meta-analysis
  results <- metabin(
    Et, 
    Nt, 
    Ec, 
    Nc, 
    data = df, 
    studlab = study, 
    comb.fixed = TRUE, 
    comb.random = FALSE, 
    sm = "OR", 
    method = "Inverse", 
    method.tau = "DL", 
    prediction = FALSE, 
    hakn = FALSE
  )
}
end_time <- Sys.time()
print(paste0(n, ',', end_time - start_time))

# visualize the last one
forest.meta(results)



######################################################################
# 1,000 to 10,000 outcomes
######################################################################
nn = (1:10)*1000

for (n in nn) {
  start_time <- Sys.time()
  for (i in 0:(n-1)) {
    # get the dataf.rame
    ocn = sprintf('M%04d', i)
    df <- all_data[all_data$outcome == ocn, c('study', 'Et', 'Nt', 'Ec', 'Nc')]
    # run meta-analysis
    results <- metabin(
      Et, 
      Nt, 
      Ec, 
      Nc, 
      data = df, 
      studlab = study, 
      comb.fixed = TRUE, 
      comb.random = FALSE, 
      sm = "OR", 
      method = "Inverse", 
      method.tau = "DL", 
      prediction = FALSE, 
      hakn = FALSE
    )
  }
  end_time <- Sys.time()
  print(paste0(n, ',', end_time - start_time))
}

######################################################################
# 10 rounds for random 1,000 outcomes
######################################################################
for (r in 1:50) {
  # get 1,000 outcomes
  n = sample(0:9999, size=1000)
  # begin test
  start_time <- Sys.time()
  for (i in n) {
    # get the dataf.rame
    ocn = sprintf('M%04d', i)
    df <- all_data[all_data$outcome == ocn, c('study', 'Et', 'Nt', 'Ec', 'Nc')]
    # run meta-analysis
    results <- metabin(
      Et, 
      Nt, 
      Ec, 
      Nc, 
      data = df, 
      studlab = study, 
      comb.fixed = TRUE, 
      comb.random = FALSE, 
      sm = "OR", 
      method = "Inverse", 
      method.tau = "DL", 
      prediction = FALSE, 
      hakn = FALSE
    )
  }
  end_time <- Sys.time()
  print(paste0(r, ',', end_time - start_time))
}





######################################################################
# Incidence Analysis
######################################################################

######################################################################
# Single round for just 1,000 outcomes
######################################################################
# set the n for number of tests
n = 1000
start_time <- Sys.time()
for (i in 1:n) {
  # get the dataf.rame
  ocn = sprintf('M%04d', i)
  df <- all_data[all_data$outcome == ocn, c('study', 'Et', 'Nt')]
  # run meta-analysis
  results <- metaprop(
    Et,
    Nt,
    data = df,
    studlab = study,
    sm = 'PLOGIT',
    comb.fixed = TRUE,
    comb.random = FALSE,
    method = 'Inverse'
  )
}
end_time <- Sys.time()
print(paste0(n, ',', end_time - start_time))

# visualize the last one
forest.meta(results)



######################################################################
# 10 rounds for 1,000 to 10,000 outcomes
######################################################################
nn = (1:10)*1000

for (n in nn) {
  start_time <- Sys.time()
  for (i in 0:(n-1)) {
    # get the dataf.rame
    ocn = sprintf('M%04d', i)
    df <- all_data[all_data$outcome == ocn, c('study', 'Et', 'Nt')]
    # run meta-analysis
    results <- metaprop(
      Et,
      Nt,
      data = df,
      studlab = study,
      sm = 'PLOGIT',
      comb.fixed = TRUE,
      comb.random = FALSE,
      method = 'Inverse'
    )
  }
  end_time <- Sys.time()
  print(paste0(n, ',', end_time - start_time))
}

######################################################################
# 10 rounds for random 1,000 outcomes
######################################################################
for (r in 1:50) {
  # get 1,000 outcomes
  n = sample(0:9999, size=100)
  # begin test
  start_time <- Sys.time()
  for (i in n) {
    # get the dataf.rame
    ocn = sprintf('M%04d', i)
    df <- all_data[all_data$outcome == ocn, c('study', 'Et', 'Nt')]
    # run meta-analysis
    results <- metaprop(
      Et,
      Nt,
      data = df,
      studlab = study,
      sm = 'PLOGIT',
      comb.fixed = TRUE,
      comb.random = FALSE,
      method = 'Inverse'
    )
  }
  end_time <- Sys.time()
  print(paste0(r, ',', end_time - start_time))
}