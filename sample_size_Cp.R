# ===========================================================================
# 
# Функция определения объёма выборки, необходимого для оценки индексов
# воспроизводимости процесса с требуемой точностью [1].
# 
# Аргументы функции:
#    type - строка символов, указывающая тип индекса воспроизводимости,
#           доступные значения: "Cp", "Cpk" и "Cpm", по умолчанию, 
#           используется значение "Cp";
#   ratio - отношение нижнего доверительного предела к выборочное оценке
#           соответствующего индекса, по умолчанию используется значение
#           0.8, то есть нижний доверительный предел будет составлять 80%
#           от найденной оценки индекса;
#   conf.level - доверительный уровень в виде десятичной дроби, по умолчанию
#                используется значение 0.95;
#    Cpk - ожидаемое значение оценки меньшего индекса воспроизводимости, 
#          необходимое для оопредедения объёма выборки для индекса данного
#          типа, по умолчанию используется значение 1;
#   bias - смещение текущего среднего уровня процесса относительно целевого
#          значения в единицах стандартного отклонения - используеться при
#          определении объёма выборки для оценки индекса Cpm, по умолчанию
#          используется значение 0.
#          
# Возвращаемое значение - значение требуемого обёма выборки n.
# 
# Аргумент type должен быть только строкой символов, остальные аргументы могут
# быть как скалярами, так и векторами. В последнем случае функция возвращает
# дата-фрейм со значениями объёмов выборок, в котором строки соответствую 
# каждой комбинации значений аргументов.
#
# Литература
#   [1] Franklin LeRoy A (1999) Sample Size Determination for Lower 
#   Confidence Limits for Estimating Process Capability Indices.
#   Computers & Industrial Engineering, Vol. 36, pp. 603 – 614.
# 
# ==========================================================================

sample_size_Cp <- function(type = c("Cp", "Cpk", "Cpm"), ratio = 0.8,
                           conf.level = 0.95, Cpk = 1, bias = 0) {
  if (length(type) > 1) {
    type = "Cp"
  }
  switch(type, 
         Cp = {if (length(ratio) > 1 & length(conf.level) > 1) {
           df <- expand.grid(ratio = ratio, conf.level = conf.level)
           df$z <- qnorm(1 - df$conf.level)
           df$n <- ceiling(1 + (2/9)*1/(df$z/2 + sqrt(1 + (df$z/2)^2 - 
                                                        (df$ratio)^(2/3)))^2)
           return(df[, -3])
         } else {
           z <- qnorm(1 - conf.level)
           n <- 1 + (2/9)*
             1/(z/2 + sqrt(1 + (z/2)^2 - (ratio)^(2/3)))^2
           if (length(ratio) > 1) {
             names(n) <- ratio
           } else if (length(conf.level) > 1) {
             names(n) <- conf.level
           }
           return(ceiling(n))
         }},
         Cpk = {if (sum(c(length(Cpk),length(ratio),
                          length(conf.level)) > 1) > 1) {
           df <- expand.grid(Cpk = Cpk, ratio = ratio,
                             conf.level = conf.level)
           df$z <- qnorm(1 - df$conf.level)
           df$n <- ceiling(df$z^2*(1/(9*df$Cpk^2) + 1/2)/(1 - df$ratio)^2)
           return(df[, -4])
         } else {
           z <- qnorm(1 - conf.level)
           n <- z^2*(1/(9*Cpk^2) + 1/2)/(1 - ratio)^2
           if (length(Cpk) > 1) {
             names(n) <- Cpk
           } else if (length(ratio) > 1) {
             names(n) <- ratio
           } else if (length(conf.level) > 1) {
             names(n) <- conf.level
           }
           return(ceiling(n))
         }},
         Cpm = {if (sum(c(length(bias), length(ratio),
                          length(conf.level)) > 1) > 1) {
           df <- expand.grid(bias = bias,
                             ratio = ratio, conf.level = conf.level)
           df$z <- qnorm(1 - df$conf.level)
           df$n <- ceiling((1 + 2*df$bias^2)/(1 + df$bias^2)^2*(2/9)*
                          1/(df$z/2 + sqrt(1 + (df$z/2)^2 - df$ratio^(2/3)))^2)
           return(df[, -4])
         } else {
           z <- qnorm(1 - conf.level)
           n <- (1 + 2*bias^2)/(1 + bias^2)^2*(2/9)*
             1/(z/2 + sqrt(1 + (z/2)^2 - ratio^(2/3)))^2
           if (length(bias) > 1) {
             names(n) <- bias
           } else if (length(ratio) > 1) {
             names(n) <- ratio
           } else if (length(conf.level) > 1) {
             names(n) <- conf.level
           }
           return(ceiling(n))
         }})
}
