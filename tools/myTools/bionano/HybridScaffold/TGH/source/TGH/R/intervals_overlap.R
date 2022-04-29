#' Check if a point overlaps an/some interval(s)
#' @description Return True or False if a point overlaps an interval or
#' any interval in an interval list.
#' If more than one interval, can also return the indices of the interval list where
#' the result is True.
#'
#' @param locus A scalar, i.e. point
#' @param intervals A vector of length 2, i. e. interval; or a matrix/dataframe/interval
#' of an interval::Interval class that each row is an interval
#' @param idxFlag If T, return the index of the interval list instead of T/F. Default is F.
#' @return T/F or a vector of indices
#' @examples
#' point_overlap_interval(13569, c(5e3, 6e3))
#' intervalList <- matrix(sort(runif(100, 0, 3e7)), ncol=2)
#' point_overlap_interval(3e6, intervalList, idxFlag = T)
#' @export
point_overlap_interval <- function(locus, intervals, idxFlag = F) {

  flag=F; i=0; is<-NULL

  if (class(intervals)=="numeric") {
    intervals <- matrix(intervals, ncol = 2)
  }

  if(nrow(intervals) > 0 & length(locus) > 0) {
    while(i<nrow(intervals)) {
      i=i+1
      interval <- intervals[i, ]
      if (locus > min(interval) & locus < max(interval)) {
        flag=T
        is <- c(is, i)
      }
    } # while
  } # if there are intervals

  if(idxFlag){
    return(is)
  }

  return(flag)
}

#' Check if an interval overlaps a/some point(s)
#' @description Return True or False if an interval overlaps a point or
#' any point in a vector of points.
#' If more than one point, can also return the indices of the interval list where
#' the result is True.
#'
#' @param loci A scalar, i.e. point, or a vector of points.
#' @param interval A vector of length 2, i. e. interval; or a matrix/dataframe/interval
#' of an interval::Interval class of one row and two columns
#' @param idxFlag If T, return the index of the point vector instead of T/F. Default is F.
#' @return T/F or a vector of indices
#' @examples
#' interval_overlap_point(13569, c(5e3, 6e3))
#' loci <- sort(runif(100, 0, 3e7))
#' interval_overlap_point(loci, c(5e3, 6e3), indexFlag = T)
#' @export
interval_overlap_point <- function(loci, interval, idxFlag=F) {

  i = 0; flag = F; is <- NULL
  if (length(interval) > 0 & length(loci) > 0) {
    while(i<length(loci)) {
      i=i+1
      locus <- loci[ i ]
      if (locus > min(interval) & locus < max(interval)) {
        flag=T
        is <- c(is, i)
      }
    } # while
  } # sanity check

  if(idxFlag) {return(is)} else { return(flag) }
}

#' Return T/F if a list of points overlap a list of intervals
#' @description Return True or False for each point overlaps any interval in an interval list.
#' @seealso \code{\link{point_overlap_interval}}
#'
#' @param loci A vector of points
#' @param intervalList A matrix/dataframe/interval
#' of an interval::Interval class that each row is an interval
#' @return T/F or a vector of indices
#' @examples
#' loci <- sort(runif(100, 0, 3e7))
#' intervalList <- matrix(sort(runif(100, 0, 3e7)), ncol=2)
#' point_overlap_wrapper(loci, intervalList)
#' @export
point_overlap_wrapper <- function(loci, intervalList) {

  selector <- lapply(1:length(loci), function(i) {
    locus <- loci[ i ]
    flag <- point_overlap_interval(locus, intervalList)
    return(flag)
  })

  selector <- unlist(selector)
  return(selector)
}

#' Check if two intervals overlap
#' @description Return the T/F or overlap portion for 2 intervals
#'
#' @param block1 A vector of length 2
#' @param block2 A vector of length 2
#' @param boundryFlag If True, return the boundry of the overlap portion. Default is F.
#' @param sizeFlag If True, return the overlap size. Default is F.
#' @return T/F, overlap boundry, or overlap size
#' @examples
#' intervals_overlap(c(0, 10), c(5, 8))
#' intervals_overlap(c(0, 10), c(5, 8), sizeFlag=T)
#' intervals_overlap(c(0, 10), c(5, 8), boundryFlag=T)
#' @export
intervals_overlap <- function(block1, block2, boundryFlag=F, sizeFlag=F) {
  block1 <- range(block1)
  block2 <- range(block2)

  flag_overlap <- ((block1[1] < block2[2]) & #check for overlap
                     (block1[2] > block2[1]))

  if (flag_overlap) {
    tmp <- sort(c(block1, block2))
    bdr <- c(tmp[2], tmp[3])
    size <- tmp[3] - tmp[2] + 1
  } else {

    bdr <- NULL
    size <- 0
  }

  if (boundryFlag) {
    return(bdr)
  } else if (sizeFlag) {
    return(size)
  } else {
    return(flag_overlap)
  }
} # intervals_overlap

#' Check if an interval overlap with another or a list of interval(s)
#' @description Return T/F, the idx, overlap interval, or overlap portion of block1
#' overlap a list of block2
#' @seealso \code{\link{intervals_overlap}}
#'
#' @inheritParams intervals_overlap
#' @param block2 A matrix/dataframe/interval::Interval class form that each row
#' is an interval
#' @param idxFlag If True, return the indices of block2 instead of T/F. Default is T.
#'
#' @return a list of T/F, overlap boundry, or overlap size
#' @examples
#' ll <- matrix(sort(runif(50, 0, 3e7)), ncol=2)
#' lll <- matrix(sort(runif(100, -3e5, 3e5)), ncol=2)
#' interval_overlap_wrapper(ll, lll)
#' interval_overlap_wrapper(ll, lll, idxFlag = F)
#' interval_overlap_wrapper(ll, lll, boundryFlag = T)
#' @export
interval_overlap_wrapper <- function(block1, block2, closeInterval=F, idxFlag=T, boundryFlag=F, sizeFlag=F) {
  block1 <- range(block1)

  # if (!is.null(nrow(block2)) && nrow(block2) > 1) { # interval, matrix, dataframe >= 2 rows
  checkO <- lapply(1:nrow(block2), function(i) {
    block2 <- as.vector(block2[i,])

    if (!idxFlag) {
      boundryFlag = F; sizeFlag = F
    }
    checkO <- intervals_overlap(block1, block2, boundryFlag, sizeFlag)
    return(checkO)
  })
  names(checkO) <- 1:nrow(block2)

  if(!idxFlag) {
    checkO <- ifelse(length(which(unlist(checkO)))==0, F, T)
  } else if (idxFlag & !boundryFlag & !sizeFlag) {
    checkO <- as.vector(which(unlist(checkO)))
  } else if(boundryFlag | sizeFlag) {
    checkO <- checkO[lapply(checkO, length)>0]
  }
  # } else { # a vector
  #
  #     checkO <- intervals_overlap(block1, block2, boundryFlag, sizeFlag)
  # }

  return(checkO)
} # interval_overlap_wrapper

#' Check if a set of point(s) or interval(s) overlap with another set of point(s) or interval(s)
#' @description Return T/F, the indices, intersections, or overlap sizes
#' @param block1 A vector/matrix/dataframe/interval::Interval class form that each row
#' is an interval
#' @param block2 A vector/matrix/dataframe/interval::Interval class form that each row
#' is an interval
#' @param closeInterval If True then the interval is closed. Otherwise it is open.
#' @param idxFlag If True then return the indices of block2 instead of T/F. Default is True.
#' @return a list of T/F, overlap boundries, or overlap sizes
#' @examples
#' ll <- matrix(sort(runif(20, 0, 100)), ncol=2)
#' lll <- matrix(sort(runif(100, -100, 100)), ncol=2)
#' intervals_overlap_wrapper(ll, lll)
#' intervals_overlap_wrapper(ll, lll, idxFlag = F)
#' intervals_overlap_wrapper(ll, lll, boundryFlag = T)
#' @export
intervals_overlap_wrapper <- function(block1, block2, closeInterval=F, idxFlag=T,
                                      boundryFlag=F, sizeFlag=F) {

  types <- ifelse(all(as.integer(as.matrix(block1))==unlist(block1)) &
                    all(as.integer(as.matrix(block2))==unlist(block2)), "Z", "R") # R doesn't matter close or open

  for (i in c("block1", "block2")) { # transform to class Interval
    block <- get( i )

    if(class(block)!="Intervals" || class(block)!="Intervals_full") {

      if(is.null(nrow(block))) { # vector
        block <- cbind(block, block)
        closeInterval <- T
      }
      block <- matrix(as.matrix(block), ncol = 2)
      block <- cbind(apply(block, 1, min), apply(block, 1, max))
      block <- intervals::Intervals(block, close = closeInterval, type = types)
      assign(i, block)
    }
  }

  if(intervals::type(block1)!=intervals::type(block2)) {
    intervals::type(block1) <- types
    intervals::type(block2) <- types
  }

  rownames(block1) <- 1:nrow(block1)
  overlaps <- intervals::interval_overlap(block1, block2, check_valid = T)
  #out <- which(sapply(overlaps, length) > 0)
  out <- overlaps

  if(!idxFlag) {
    out <- sapply(overlaps, length) > 0
  }

  if(boundryFlag | sizeFlag) {
    out <- lapply(1:length(overlaps), function( i ) {
      idcs <- overlaps[[ i ]]
      if(length(idcs)==0) {return(NULL)}

      eachIdx <- lapply(idcs, function( idx ) {
        intersection <- intervals::interval_intersection(block1[i,], block2[idx,])
        out <- c("start"=intersection[,1], "end"=intersection[,2])

        if (sizeFlag) {
          siza <- intervals::size(intersection)
          out <- c("size"=siza)
        }
        return(out)
      })
      eachIdx <- cbind(idcs, do.call(rbind, eachIdx))

      return(eachIdx)
    })
    names(out) <- names(overlaps)
  } # if

  return(out)
}
