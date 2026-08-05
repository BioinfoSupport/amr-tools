
library(GenomicAlignments)
library(tidyverse)


gaCigarBedpe <- function(aln) {
	stopifnot("aln must be named"=!is.null(names(aln)))
	
	ref_rg <- cigarRangesAlongReferenceSpace(
		cigar(aln), pos = start(aln),
		ops = c("M", "=", "X"), reduce.ranges = FALSE,f=as.factor(seqnames(aln))
	) |>
		GRanges(seqinfo = seqinfo(aln))
	
	
	qry_len <- cigarWidthAlongQuerySpace(cigar(aln), before.hard.clipping = TRUE)
	qry_rg <- cigarRangesAlongQuerySpace(
		cigar(aln),
		ops = c("M", "=", "X"), reduce.ranges = FALSE,
		before.hard.clipping = TRUE
	) |>
		setNames(names(aln))
	aln_idx <- rep(seq_along(aln),lengths(qry_rg))
	
	qry_rg <- GRanges(
		qry_rg,
		strand = rep(strand(aln),lengths(qry_rg)),
		seqinfo = data.frame(seqnames=names(aln),seqlengths=qry_len,stringsAsFactors = FALSE) |> distinct() |> as("Seqinfo")
	)
	ranges(qry_rg)[strand(qry_rg)=="-"] <- reflect(
		ranges(qry_rg)[strand(qry_rg)=="-"],
		IRanges(1,seqlengths(qry_rg)[as.integer(seqnames(qry_rg[strand(qry_rg)=="-"]))])
	)
	
	DataFrame(aln_idx,ref = ref_rg,query = qry_rg)
}



bedpePolygons <- function(bedpe) {
	bind_rows(
		"r" = as_tibble(bedpe$ref) |> tibble::rowid_to_column("polygon_id") |> bind_cols(aln_idx = bedpe$aln_idx),
		"q" = as_tibble(bedpe$query) |> tibble::rowid_to_column("polygon_id") |> bind_cols(aln_idx = bedpe$aln_idx),
		.id = "source"
	) |>
		select(-width) |>
		dplyr::rename("1"="start","2"="end") |>
		pivot_longer(cols = c("1","2"),names_to = "order",names_transform = list(order=as.integer)) |>
		mutate(order = case_when(
			source=="r" ~ order,
			source=="q" & strand!="-" ~ 5 - order,
			source=="q" & strand=="-" ~ 2 + order,
			TRUE ~ NA_integer_
		)) |>
		select(-strand,-source) |>
		arrange(polygon_id,order)	
}





rotate_queries <- function(qry_fa,aln_bam) {
	aln <- readGAlignments(aln_bam,use.names=TRUE)
	
	# For each query, determine the ref it cover most
	qry_map <- split(aln,names(aln)) |> 
		lapply(\(e) sum(coverage(e)>0)) |>
		enframe(name="query_id") |>
		unnest_longer(value,indices_to="ref_id",values_to="coverage") |>
		group_by(query_id) |>
		slice_max(coverage,n=1,with_ties=FALSE)
	
	# Filter BAM to keep best contig pairs, and longest alignments
	aln <- subset(aln,paste(names(aln),seqnames(aln)) %in% paste(qry_map$query_id,qry_map$ref_id))
	aln <- aln[order(width(aln),decreasing=TRUE)]
	aln <- aln[!duplicated(names(aln))]
	
	# Use selected alignment to rotate queries
	rot <- range(gaCigarBedpe(aln)$query) |> resize(width=0L)
	fa <- readDNAStringSet(qry_fa)
	names(fa) <- str_replace(names(fa)," .*","")
	
	FA <- fa
	FA[seqnames(rot)] <- xscat(subseq(FA[seqnames(rot)],start=start(rot)),subseq(FA[seqnames(rot)],end=end(rot)))
	FA[seqnames(rot)][strand(rot)=="-"] <- reverseComplement(FA[seqnames(rot)][strand(rot)=="-"])
	#names(FA) <- paste0(names(fa),"_",start(rot),"_",strand(rot))
	FA
}


