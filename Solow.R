# Solow model calculated in total units

library(data.table)
library(ggplot2)

# Can use any production function that is linearly homogenous and satisfies Inada conditions
# I will use a Cobb-Douglas function in which the technological parameter affects only L
cobb_douglas = function (A, L, K, alpha=1/3){
  return (K^alpha*(A*L)^(1-alpha))
}

# Basic structure of the state_variables and parameters
# Past states will be stored to document the path
sample_state <- data.frame(A = 2, L = 100, K = 100)

sample_parameters <- list(Y = cobb_douglas, s = 0.2, delta = 0.1, n = 0.03, g = 0.02)


# Calculates 1 time step using the formula K_dot = s*Y(A,L,K) - delta*K
iterate = function(state_variables, parameters){
  last_state = state_variables[nrow(state_variables),]
  next_state = data.frame(A=last_state$A*(1+parameters$g), L=last_state$L*(1+parameters$n),
                          K=last_state$K-parameters$delta*last_state$K
                          +parameters$s*parameters$Y(last_state$A,last_state$L,last_state$K))
  return (rbind(state_variables, next_state))
}


# Iterates the model until the growth rate of capital stock converges
calculate_to_equilibrium = function(state_variables, parameters, cutoff = 1e-6) {
  growth = function(){
    n = nrow(state_variables)
    if (n==1) {
      return (100)
    }
    return (state_variables[n, "K"]/state_variables[n-1, "K"])
  }
  last_growth = 1
  while (abs(last_growth - growth()) > cutoff) {
    last_growth = growth()
    state_variables = iterate(state_variables, parameters)
  }
  return (state_variables)
}

sample_path = calculate_to_equilibrium(sample_state, sample_parameters)

# just iterate a fixed number of times
calculate = function (state_variables, parameters, n=40) {
  for (i in (1:n)) {
    state_variables = iterate(state_variables, parameters)
  }
  return (state_variables)
}


# The data on the calculated convergence path can be expanded with derived 
# values we might be interested in
expand_data = function(path, Y=cobb_douglas) {
  path$t=c(1:nrow(path))
  path$Y=Y(path$A, path$L, path$K)
  path$Y_dot=shift(path$Y, -1)-path$Y
  path$Y_hat=path$Y_dot/path$Y
  path$y=path$Y/path$L
  path$y_dot=shift(path$y, -1)-path$y
  path$y_hat=path$y_dot/path$y
  path$y_tilde=path$Y/path$L/path$A
  path$y_tilde_dot=shift(path$y_tilde, -1)-path$y_tilde
  path$y_tilde_hat=path$y_tilde_dot/path$y_tilde
  path$K_dot=shift(path$K, -1)-path$K
  path$K_hat=path$K_dot/path$K
  path$k=path$K/path$L
  path$k_dot=shift(path$k, -1)-path$k
  path$k_hat=path$k_dot/path$k
  path$k_tilde=path$K/path$L/path$A
  path$k_tilde_dot=shift(path$k_tilde, -1)-path$k_tilde
  path$k_tilde_hat=path$k_tilde_dot/path$k_tilde
  return (path)
}

sample_expanded = expand_data(sample_path)



# Plot no-growth-model (s=0.3, delta=0.4)

ng_plot1 = ggplot() + theme_minimal()
ng_plot1 = ng_plot1 + labs(title="No-Growth equilibrium")
ng_plot1 = ng_plot1 + geom_function(fun=function(k){return (k^0.33*100)}, size=1.4)
ng_plot1 = ng_plot1 + geom_function(fun=function(k){return (k^0.33*0.3*100)}, col="green", size=1.4)
ng_plot1 = ng_plot1 + geom_function(fun=function(k){return (k*0.4*100)}, col="red", size=1.4)

ng_plot1

# Calculate no-growth-model (s=0.3, delta=0.4)

ng_parameters = list(Y = cobb_douglas, s = 0.3, delta = 0.4, n = 0, g = 0)

ng_state1 = data.frame(A = 1, L = 100, K = 0.001)
ng_expanded1 = expand_data(calculate(ng_state1, ng_parameters))
ng_expanded1$Scenario = "Low starting capital"
ng_state2 = data.frame(A = 1, L = 100, K = 200)
ng_expanded2 = expand_data(calculate(ng_state2, ng_parameters))
ng_expanded2$Scenario = "High starting capital"
ng_state3 = data.frame(A = 1, L = 100, K = (0.3/0.4)^(3/2)*100)
ng_expanded3 = expand_data(calculate(ng_state3, ng_parameters))
ng_expanded3$Scenario = "Equilibrium starting capital"

ng_allPaths = rbind(ng_expanded1, ng_expanded2, ng_expanded3)


# Plot Y over time

ng_plot3 = ggplot(ng_allPaths, aes(t,Y, group=Scenario, colour=Scenario)) + theme_minimal()
ng_plot3 = ng_plot3 + labs(title="No-Growth output paths")
ng_plot3 = ng_plot3 + geom_line()

ng_plot3

# Plot k_dot over k

ng_plot2 = ggplot() + theme_minimal()
ng_plot2 = ng_plot2 + labs(title="No-Growth net capital accumulation", subtitle="Calculated using discrete-time modeling")
ng_plot2 = ng_plot2 + geom_line(aes(k, k_dot), ng_expanded1)
ng_plot2 = ng_plot2 + geom_line(aes(k, k_dot), ng_expanded2)
ng_plot2 = ng_plot2 + geom_hline(aes(yintercept=0))

ng_plot2

# Calculate and plot growth-model with variable starting capital

g_parameters = list(Y = cobb_douglas, s = 0.3, delta = 0.4, n = 0.05, g = 0.02)

g_state1 = data.frame(A = 1, L = 100, K = 300)
g_expanded1 = expand_data(calculate(g_state1, g_parameters))
g_expanded1$Scenario = "High starting capital"
g_state2= data.frame(A = 1, L = 100, K = 1)
g_expanded2 = expand_data(calculate(g_state2, g_parameters))
g_expanded2$Scenario = "Low starting capital"
g_state3 = data.frame(A = 1, L = 100, K = (0.3/0.4)^(3/2)*100)
g_expanded3 = expand_data(calculate(g_state3, g_parameters))
g_expanded3$Scenario = "Equilibrium starting capital"

g_allPaths = rbind(g_expanded1, g_expanded2, g_expanded3)


g_plot1 = ggplot(g_allPaths, aes(t,Y, group=Scenario, colour=Scenario)) + theme_minimal()
g_plot1 = g_plot1 + labs(title="Growth-model output paths")
g_plot1 = g_plot1 + geom_line()

g_plot1

# What follows is stupid and can be solved a million times more efficiently
# I am iterating the model until equilibrium in order to deduce the effect of the savings rate

golden_rule_parameters = list(Y = cobb_douglas, s = 0, delta = 0.2, n = 0, g = 0)
golden_rule_starting_state = data.frame(A = 1, L = 100, K = 100)

golden_rule_data = data.frame(s=NULL, C=NULL, delta=NULL)

for(delta in c(0.1, 0.2, 0.3)){
  golden_rule_parameters$delta=delta
  for(s in c((1:100)/100)){
    golden_rule_parameters$s=s
    path=calculate_to_equilibrium(golden_rule_starting_state, golden_rule_parameters, 1e-4)
    lastY = golden_rule_parameters$Y(path$A, path$L, path$K)[nrow(path)]
    golden_rule_data = rbind(golden_rule_data, data.frame(s=s,C=lastY*(1-s),delta=golden_rule_parameters$delta))
  }
}


gr_plot = ggplot(golden_rule_data, aes(s,C, group=delta, colour=delta)) + theme_minimal()
gr_plot = gr_plot + labs(title="No-Growth optimal savings rate")
gr_plot = gr_plot + geom_line()

gr_plot
