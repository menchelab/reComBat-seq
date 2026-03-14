# Adapted from edgeR 4.0.1 (DGEList.R)
# Original authors: edgeR team (see edgeR package).
# .isAllZero is not exported by edgeR and is included here
# as it is required internally by reComBatseq.

.isAllZero <- function(y) 
# Check whether all counts are zero.
# Also checks and stops with an informative error message if negative, NA or infinite counts are present.
{
	if (!length(y)) return(FALSE)
	check.range <- range(y)
	if (is.na(check.range[1])) stop("NA counts not allowed", call.=FALSE)
	if (check.range[1] < 0) stop("Negative counts not allowed", call.=FALSE)
	if (is.infinite(check.range[2])) stop("Infinite counts not allowed", call.=FALSE)
	check.range[2]==0
}
