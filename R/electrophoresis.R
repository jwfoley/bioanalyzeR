#' Combine multiple electrophoresis objects
#'
#' This function combines multiple \code{electrophoresis} objects into one so you can analyze and graph multiple batches together.
#'
#' All data frames are combined by \code{\link{rbind}} and lists are combined by \code{\link{c}}. Factor levels are expanded to the union of all inputs. The \code{data$peak} column is renumbered to match the new peak table.
#'
#' @param ... Two or more objects of class \code{electrophoresis}.
#'
#' @return A new \code{electrophoresis} object containing all the data from the previous ones in the provided order.
#'
#' @export
rbind.electrophoresis <- function(...) {
	arg.list <- list(...)
	
	# increment the peak indexes in the data so they'll match the new table
	for (i in 1:(length(arg.list) - 1)) for (j in (i + 1):length(arg.list)) arg.list[[j]]$data$peak <- arg.list[[j]]$data$peak + nrow(arg.list[[i]]$peaks)
	
	structure(list(
		data = do.call(rbind, lapply(arg.list, function(x) x$data)),
		assay.info = do.call(c, lapply(arg.list, function(x) x$assay.info)),
		samples = do.call(rbind, lapply(arg.list, function(x) x$samples)),
		wells.by.ladder = do.call(c, lapply(arg.list, function(x) x$wells.by.ladder)),
		peaks = do.call(rbind, lapply(arg.list, function(x) x$peaks)),
		regions = do.call(rbind, lapply(arg.list, function(x) x$regions)),
		mobility.functions = do.call(c, lapply(arg.list, function(x) x$mobility.functions)),
		mass.coefficients = do.call(c, lapply(arg.list, function(x) x$mass.coefficients))
	), class = "electrophoresis")
}

