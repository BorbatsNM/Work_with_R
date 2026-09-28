# =========================================================================
# 
# Функция обработки данных выборочного контроля по количественному признаку 
# с использованием одноступенчатого плана при известном стандартном отклонении
# контролируемой характеристики, имеющей нормальное распределение.
# Аргументами функции являются:
#  data – числовой вектор данных контроля, длина вектора должна
#  соответствовать объёму выборки в соответствии с используемым
#  планом контроля.
#  limits – двухэлементный список границ допуска (спецификации), первый элемент
#  соответствует нижней границе L, а второй элемент – верхней границе U.
#  В случае наличия единственной границы, второму элементу должно быть
#  присвоено значение NULL.
#  sigma – числовое значение известного стандартного отклонения контролируемой
#  характеристики качества.
#  k – числовое значение контрольного норматива в соответствии
#  с используемым планом контроля.
#  na.rm – логический аргумент, если задано значение FALSE (по умолчанию)
#  отсутствующие значения в векторе данных data удаляются с предупреждением,
#  при значении TRUE пропущенные значения удаляются без предупреждения.
# Результат функции представляет собой сообщение, выводимое в консоль,
# и содержащее решение относительно партии и всю необходимую информацию.
# В тоже время, если присвоить результат функции переменной,
# то она будет содержать список значений исходных параметров 
# и всех вычисленных величин.

result_SV_sigma <- function(data, limits, sigma, k, na.rm = FALSE) {
  stopifnot(exprs = {
    is.vector(data, mode = "numeric")
    length(limits) == 2
    length(sigma) == 1
    length(k) == 1
    length(na.rm) == 1
    is.numeric(sigma)
    is.numeric(k)
    is.logical(na.rm)
  })
  
  if (all(sapply(limits, is.null))) {
    stop("Должна быть задана хотя бы одна граница допуска.")
  }
  
  n <- length(data)
  
  if (any(is.na(data))) {
    avg <- mean(data, na.rm = TRUE)
    if (!na.rm) {
      warning(paste0("Удалены отсутствующие значения.", "\n",
                    "Число фактических наблюдений меньше объёма выборки!"))
    }
  } else {
    avg <- mean(data)
  }
  
  p_max <- pnorm(k*sqrt(n/(n - 1)), lower.tail = FALSE)
  
  # обработка результатов при едиственной границе допуска
  if (any(sapply(limits, is.null))) {
    if (is.null(limits[[1]])) {
      border <- limits[[2]]
      info <- " (верхняя) U:"
      A <- border - k*sigma
      if (avg <= A) {
        decision <- "Партию следует принять."
      } else {
        decision <- "Партию следует отклонить."
      }
      Q <- (border - avg)/sigma * sqrt(n/(n - 1))
      p_est <- pnorm(Q, lower.tail = FALSE)
    } else {
      border <- limits[[1]]
      info <- " (нижняя) L:"
      A <- border + k*sigma
      if (avg >= A) {
        decision <- "Партию следует принять."
      } else {
        decision <- "Партию следует отклонить."
      }
      Q <- (avg - border)/sigma * sqrt(n/(n - 1))
      p_est <- pnorm(Q, lower.tail = FALSE)
    }
    
    result <- list(limit = border, k = k, sigma = sigma, Q = Q, p_max = p_max,
                   p = p_est, mean = avg, A = A, n = n)
    
    cat("\n\t", decision, "\n\n")
    cat("Объём выборки:\t\t\t\t", n, "\n", sep = "")
    cat("Контрольный норматив:\t\t\t", k, "\n", sep = "")
    cat("Известное стандартное отклонение:\t", sigma, "\n", sep = "")
    cat("Граница допуска", info, "\t\t", border, "\n", sep = "")
    cat("Среднее по выборке:\t\t\t", format(avg, digits = 3, nsmall = 3),
        "\n", sep = "")
    cat("Граница приёмки:\t\t\t", format(A, digits = 3, nsmall = 3),
        "\n", sep = "")
    cat("Максимально допустимое значение\nуровня несоответствий: \t\t\t",
        format(p_max*100, digits = 3, nsmall = 3),
        "% \n", sep = "")
    cat("Оценка уровня несоответствий\nпо выборке: \t\t\t\t",
        format(p_est*100, digits = 3, nsmall = 3),
        "%\n\n", sep = "")
    
  } else {  # обработка результатов при обеих границах допуска
    A_L <- limits[[1]] + k*sigma
    A_U <- limits[[2]] - k*sigma
    Q_L <- (avg - limits[[1]])/sigma * sqrt(n/(n - 1))
    p_L <- pnorm(Q_L, lower.tail = FALSE)
    Q_U <- (limits[[2]] - avg)/sigma * sqrt(n/(n - 1))
    p_U <- pnorm(Q_U, lower.tail = FALSE)
    p_est <- p_L + p_U
    
    if (p_est <= p_max) {
      decision <- "Партию следует принять."
    } else {
      decision <- "Партию следует отклонить."
    }
    
    result <- list(limits = limits, k = k, sigma = sigma, Q_L = Q_L,
                   Q_U = Q_U, p_L = p_L, p_U = p_U, p = p_est,
                   p_max = p_max, mean = avg, n = n, A_L = A_L, A_U = A_U)
    
    cat("\n\t", decision, "\n\n")
    cat("Объём выборки:\t\t\t\t", n, "\n", sep = "")
    cat("Контрольный норматив:\t\t\t", k, "\n", sep = "")
    cat("Известное стандартное отклонение:\t", sigma, "\n", sep = "")
    cat("Нижняя граница допуска (L):\t\t", limits[[1]], "\n", sep = "")
    cat("Верхняя граница допуска (U):\t\t", limits[[2]], "\n", sep = "")
    cat("Среднее по выборке:\t\t\t", format(avg, digits = 3, nsmall = 3),
        "\n", sep = "")
    cat("Нижняя граница приёмки:\t\t\t", format(A_L, digits = 3, nsmall = 3),
        "\n", sep = "")
    cat("Верхняя граница приёмки:\t\t", format(A_U, digits = 3, nsmall = 3),
        "\n", sep = "")
    cat("Максимально допустимое значение\nуровня несоответствий: \t\t\t",
        format(p_max*100, digits = 3, nsmall = 3),
        "% \n", sep = "")
    cat("Оценка уровня несоответствий\nпо выборке: \t\t\t\t",
        format(p_est*100, digits = 3, nsmall = 3),
        "%\n\n", sep = "")
  }
  return(invisible(result))
}
