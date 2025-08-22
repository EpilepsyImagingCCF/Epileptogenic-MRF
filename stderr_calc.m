function stderr = stderr_calc(data)

stderr = std(data)/sqrt(length(data));        