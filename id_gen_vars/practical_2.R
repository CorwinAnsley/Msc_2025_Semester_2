
y = rnorm(1000,5,1)
hist(y)

var_1 = rnorm(1000,4,1)
plot(var_1, y)

summary(lm(y~var_1))

summary(lm(y~var_1))$coefficients

y = rnorm(1000,5,1)
var_1 = rnorm(1000,4,1)
summary(lm(y~var_1))

plot(var_1, y)

plessthan0.05 <- rep(0,100)
for (i in (1:100))
{
  y2 <- rnorm(1000,5,1)
  Var2 <- rnorm(1000,4,1)
  mod <- lm(y2~Var2)
  if (summary(mod)$coefficients[2,4] < 0.05) {
    plessthan0.05[i] <- 1
  }	
}
sum(plessthan0.05)/100

test_data = data.frame(y)

for (i in (1:100))
{
  Var_x <- rnorm(1000,4,1)
  test_data[,ncol(test_data) + 1] = Var_x
}

generate_test = function(iterations = 15) {
  y_t = rnorm(1000,5,1)
  test_data = data.frame(y_t)
  
  for (i in (1:iterations))
  {
    Var_x <- rnorm(1000,4,1)
    test_data[,ncol(test_data) + 1] = Var_x
  }
  
  for (j in (1:iterations))
  {
    if (summary(lm(y_t~test_data[,j+1],data=test_data))$coefficients[2,4] < 
        0.05) {
      plessthan0.05[j] <- 1
    }
  }
  return(plessthan0.05)
}

plessthan0.05 = generate_test()
sum(plessthan0.05)/100

p_vals <- rep(0,100)

for (i in (1:100)) {
  plessthan0.05 = generate_test()
  
  p_vals[i] = sum(plessthan0.05)/100
}

sum(p_vals)/100
