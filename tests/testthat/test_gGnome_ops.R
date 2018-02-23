
library(gGnome)
library(testthat)
library(gUtils)

context('testing gGnome')


message("Toy segments: ", system.file('extdata', 'testing.segs.rds', package="gGnome"))
test_segs = readRDS(system.file('extdata', 'testing.segs.rds', package="gGnome"))

message("Toy edges: ", system.file('extdata', 'testing.es.rds', package="gGnome"))
test_es = readRDS(system.file('extdata', 'testing.es.rds', package="gGnome"))

jab = system.file('extdata', 'jabba.simple.rds', package="gGnome")
message("JaBbA result: ", jab)

prego = system.file('extdata', 'intervalFile.results', package='gGnome')
message("PREGO results: ", prego)

weaver = system.file('extdata', 'weaver', package='gGnome')
message("Weaver results: ", weaver)


jabba = readRDS(jab)

junctions = jabba$junctions
## 

gr = GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c("A","B","C"))
gr2 = GRanges(1, IRanges(c(1,9), c(6,14)), strand=c('+','-'), seqinfo=Seqinfo("1", 25), field=c(1,2))
dt = data.table(seqnames=1, start=c(2,5,10), end=c(3,8,15))






test_that('junctions works', {

    expect_error(junctions(data.frame()))
    expect_error(junctions(GRanges()))
    expect_equal(length(junctions(GRangesList())), 0)
    expect_equal(length(junctions(grl1)), 250)

})




test_that('ra.duplicated works', {

    expect_false(all(ra.duplicated(junctions(grl1))))
    ##expect_equal(ra.duplicated(GRangesList()), logical(0))

})



## gGraph = R6::R6Class("gGraph"

## initialize = function(tile=NULL, junctions=NULL, cn = FALSE, jabba=NULL,
##     weaver=NULL, prego=NULL, segs=NULL, es=NULL, ploidy=NULL, purity=NULL, regular=TRUE)


## segstats, edges, junctions, G, adj, A, parts, seqinfo, purity, ploidy, td, win, ig


test_that('gGraph works', {

    ggnew = gGraph$new()
    expect_true(is(ggnew, 'gGraph'))
    foobar = gGraph$new(segs = test_segs, es=test_es)
    expect_true(is(foobar, 'gGraph'))
    foojab = gGraph$new(jabba = jab, segs = test_segs, es=test_es)
    expect_true(is(foojab, 'gGraph'))
    fooweaver = gGraph$new(weaver=weaver, segs = test_segs, es=test_es)
    expect_true(is(fooweaver, 'gGraph'))
    fooprego = gGraph$new(prego=prego, segs = test_segs, es=test_es)
    expect_true(is(fooprego, 'gGraph'))
    foocn = gGraph$new(segs = test_segs, es=test_es, cn=TRUE)
    expect_true(is(foocn, 'gGraph'))
    fooregular = gGraph$new(segs = test_segs, es=test_es, regular=FALSE)
    expect_true(is(fooregular, 'gGraph'))
    foopurityploidy = gGraph$new(segs = test_segs, es=test_es, ploidy=3334, purity=233432)
    expect_true(is(foopurityploidy, 'gGraph'))  
    ##
    ##
    added_junctions = foojab$addJuncs(readRDS(jab)$junc)
    expect_true(is(added_junctions, 'gGraph')) 

})


## FUCNTIONS: initialize, set.seqinfo, nullGGraph, simpleGraph, dipGraph, addJuncs, addSegs, karyograph, simplify,
##     decouple, add, jabb2gg, wv2gg, pr2gg, print, plot, window, layout, summary, gg2td, son, html, gg2js, components, 
##     subgraph, filling, 

## ACTIVE BINDINGS: segstats, edges, junctions, G, adj, A, parts, seqinfo, purity, ploidy, td, win, ig

test_that('gGraph works, default', {

    ggnew = gGraph$new()
    expect_true(is(ggnew, 'gGraph'))
    ## ACCESS ACTIVE BINDINGS
    expect_equal(length(ggnew$segstats), 0)
    expect_equal(dim(ggnew$edges)[1], 0)
    expect_equal(length(ggnew$junctions), 0)
    expect_error(ggnew$G, NA)  ## check it works; IGRAPH 84fc0c4 D--- 0 0 -- + edges from 84fc0c4:
    expect_equal(length(ggnew$adj), 0)
    expect_equal(length(ggnew$A), 0)
    expect_equal(ggnew$parts, NULL)
    expect_equal(length(ggnew$seqinfo), 0)
    expect_equal(ggnew$purity, NULL)
    expect_equal(ggnew$ploidy, NULL)
    expect_equal(ggnew$td, NULL)
    expect_equal(length(ggnew$win), 0)
    expect_equal(ggnew$ig, NULL)
    ## FUNCTIONS
    ## set.seqinfo = function(genome=NULL, gname=NULL, drop=FALSE)
    ggnew_setseq = ggnew$set.seqinfo()
    expect_true(is(ggnew_setseq, 'gGraph'))
    expect_equal(length(ggnew_setseq$segstats), 0)
    expect_equal(dim(ggnew_setseq$edges)[1], 0)
    expect_equal(length(ggnew_setseq$junctions), 0)
    expect_error(ggnew_setseq$G, NA)  ## check it works
    expect_equal(length(ggnew_setseq$adj), 0)
    expect_equal(length(ggnew_setseq$A), 0)
    expect_equal(ggnew_setseq$parts, NULL)
    expect_equal(length(ggnew_setseq$seqinfo), 0)
    expect_equal(ggnew_setseq$purity, NULL)
    expect_equal(ggnew_setseq$ploidy, NULL)
    expect_equal(ggnew_setseq$td, NULL)
    expect_equal(length(ggnew_setseq$win), 0)
    expect_equal(ggnew_setseq$ig, NULL)
    ## set.seqinfo, drop = TRUE
    ggnew_setseq_drop = ggnew$set.seqinfo(gname = 'foobar', drop=TRUE)
    expect_true(is(ggnew_setseq_drop, 'gGraph'))
    expect_equal(length(ggnew_setseq_drop$segstats), 0)
    expect_equal(dim(ggnew_setseq_drop$edges)[1], 0)
    expect_equal(length(ggnew_setseq_drop$junctions), 0)
    expect_error(ggnew_setseq_drop$G, NA)   ## check it works
    expect_equal(length(ggnew_setseq_drop$adj), 0)
    expect_equal(length(ggnew_setseq_drop$A), 0)
    expect_equal(ggnew_setseq_drop$parts, NULL)
    expect_equal(length(ggnew_setseq_drop$seqinfo), 0)
    expect_equal(ggnew_setseq_drop$purity, NULL)
    expect_equal(ggnew_setseq_drop$ploidy, NULL)
    expect_equal(ggnew_setseq_drop$td, NULL)
    expect_equal(length(ggnew_setseq_drop$win), 0)
    expect_equal(ggnew_setseq_drop$ig, NULL)
    ## set.seqinfo, genome != NULL, gname != NULL
    ggnew_setseq_hg = ggnew$set.seqinfo(genome = hg_seqlengths(), gname = 'foobar', drop = TRUE)
    expect_true(is(ggnew_setseq_hg, 'gGraph'))
    expect_equal(length(ggnew_setseq_hg$segstats), 0)
    expect_equal(dim(ggnew_setseq_hg$edges)[1], 0)
    expect_equal(length(ggnew_setseq_hg$junctions), 0)
    expect_error(ggnew_setseq_hg$G, NA) ## check it works
    expect_equal(length(ggnew_setseq_hg$adj), 0)
    expect_equal(length(ggnew_setseq_hg$A), 0)
    expect_equal(ggnew_setseq_hg$parts, NULL)
    expect_equal(length(ggnew_setseq_hg$seqinfo), 25)   ### checks!
    expect_equal(ggnew_setseq_hg$purity, NULL)
    expect_equal(ggnew_setseq_hg$ploidy, NULL)
    expect_equal(ggnew_setseq_hg$td, NULL)
    expect_equal(length(ggnew_setseq_hg$win), 0)
    expect_equal(ggnew_setseq_hg$ig, NULL)
    ## nullGraph = function(regular=TRUE, genome=NULL)
    ggnew_setseq_nullGraph = ggnew$nullGGraph()
    expect_true(is(ggnew_setseq_nullGraph, 'gGraph'))
    expect_equal(length(ggnew_setseq_nullGraph$segstats), 0)
    expect_equal(dim(ggnew_setseq_nullGraph$edges)[1], 0)
    expect_equal(length(ggnew_setseq_nullGraph$junctions), 0)
    expect_error(ggnew_setseq_nullGraph$G, NA) ## check it works
    expect_equal(length(ggnew_setseq_nullGraph$adj), 0)
    expect_equal(length(ggnew_setseq_nullGraph$A), 0)
    expect_equal(ggnew_setseq_nullGraph$parts, NULL)
    expect_equal(length(ggnew_setseq_nullGraph$seqinfo), 25)   ### checks! "null" means there is no node, you can still have a "space" of possible values when the set is empty
    expect_equal(ggnew_setseq_nullGraph$purity, NULL)
    expect_equal(ggnew_setseq_nullGraph$ploidy, NULL)
    expect_equal(ggnew_setseq_nullGraph$td, NULL)
    expect_equal(length(ggnew_setseq_nullGraph$win), 0)
    expect_equal(ggnew_setseq_nullGraph$ig, NULL)
    ## simpleGraph = function(genome = NULL, chr=FALSE, include.junk=FALSE, ploidy = NULL)
    ggnew_setseq_simpleGraph = ggnew$simpleGraph()
    expect_true(is(ggnew_setseq_simpleGraph, 'gGraph'))
    expect_equal(length(ggnew_setseq_simpleGraph$segstats), 50)
    expect_equal(dim(ggnew_setseq_simpleGraph$edges)[1], 0)
    expect_equal(length(ggnew_setseq_simpleGraph$junctions), 0)
    expect_error(ggnew_setseq_simpleGraph$G, NA) ## check it works
    expect_equal(length(ggnew_setseq_simpleGraph$adj), 2500)
    expect_equal(length(ggnew_setseq_simpleGraph$A), 2500)
    ## expect_equal(ggnew_setseq_simpleGraph$parts, NULL)
    expect_equal(length(ggnew_setseq_simpleGraph$seqinfo), 25)   ### checks! "null" means there is no node, you can still have a "space" of possible values when the set is empty
    expect_equal(ggnew_setseq_simpleGraph$purity, NULL)
    expect_equal(ggnew_setseq_simpleGraph$ploidy, NULL)
    expect_true(is(ggnew_setseq_simpleGraph$td, 'gTrack'))
    expect_equal((ggnew_setseq_simpleGraph$td)$ygap, 2)
    expect_match((ggnew_setseq_simpleGraph$td)$name, 'CN')
    expect_equal(length(ggnew_setseq_simpleGraph$win), 25)
    ## STILL ERROR
    ## >  ggnew_setseq_simpleGraph$ig
    ##Error in log(private$segs$cn, 1.4) : 
    ##  non-numeric argument to mathematical function
    ##In addition: Warning messages:
    ##1: replacing previous import ‘VariantAnnotation::select’ by ‘plotly::select’ when loading ‘skitools’ 
    ##2: replacing previous import ‘ggplot2::last_plot’ by ‘plotly::last_plot’ when loading ‘skitools’ 
    ## dipGraph = function(genome = NULL, chr=FALSE, include.junk=FALSE)
    ggnew_dd = ggnew$dipGraph()
    expect_true(is(ggnew_dd, 'gGraph'))
    expect_equal(length(ggnew_dd$segstats), 50)
    expect_equal(dim(ggnew_dd$edges)[1], 0)
    expect_equal(length(ggnew_dd$junctions), 0)
    expect_error(ggnew_dd$G, NA) ## check it works
    expect_equal(length(ggnew_dd$adj), 2500)
    expect_equal(length(ggnew_dd$A), 2500)
    ## > ggnew_dd$ig
    ## > ggnew_dd$parts
    expect_equal(length(ggnew_dd$seqinfo), 25)   
    expect_equal(ggnew_dd$purity, NULL)
    expect_equal(ggnew_dd$ploidy, 2)  ## checks!
    expect_true(is(ggnew_dd$td, 'gTrack'))
    expect_equal((ggnew_dd$td)$ygap, 2)
    expect_match((ggnew_dd$td)$name, 'CN')
    expect_equal(length(ggnew_dd$win), 25)
    ## dipGraph
    ggnew_dd_junk = ggnew$dipGraph(genome = hg_seqlengths(), chr=TRUE, include.junk=TRUE)
    expect_true(is(ggnew_dd_junk, 'gGraph'))
    expect_equal(length(ggnew_dd_junk$segstats), 50)
    expect_equal(dim(ggnew_dd_junk$edges)[1], 0)
    expect_equal(length(ggnew_dd_junk$junctions), 0)
    expect_error(ggnew_dd_junk$G, NA) ## check it works
    expect_equal(length(ggnew_dd_junk$adj), 2500)
    expect_equal(length(ggnew_dd_junk$A), 2500)
    ## expect_equal(ggnew_dd$parts, NULL)
    expect_equal(length(ggnew_dd_junk$seqinfo), 25)   
    expect_equal(ggnew_dd_junk$purity, NULL)
    expect_equal(ggnew_dd_junk$ploidy, 2)  ## checks!
    expect_true(is(ggnew_dd_junk$td, 'gTrack'))
    expect_equal((ggnew_dd_junk$td)$ygap, 2)
    expect_match((ggnew_dd_junk$td)$name, 'CN')
    expect_equal(length(ggnew_dd_junk$win), 25)
    ## addJuncs = function(junc, cn=TRUE)
    expect_error(ggnew$addJuncs())
    added_juncs = ggnew$addJuncs(junc = junctions)
    expect_true(is(added_juncs, 'gGraph'))
    expect_equal(length(added_juncs$segstats), 1112)
    expect_equal(dim(added_juncs$edges)[1], 1596)
    expect_equal(dim(added_juncs$edges)[2], 15)
    expect_equal(length(added_juncs$junctions), 267)
    expect_error(added_juncs$G, NA) ## check it works
    expect_equal(length(added_juncs$adj), 1236544)
    expect_equal(length(added_juncs$A),  1236544)
    ## expect_equal(added_juncs$parts, NULL)
    expect_equal(length(added_juncs$seqinfo), 25)   
    expect_equal(added_juncs$purity, NULL)
    expect_equal(added_juncs$ploidy, NULL)  ## checks!
    expect_true(is(added_juncs$td, 'gTrack'))
    expect_equal((added_juncs$td)$ygap, 2)
    expect_match((added_juncs$td)$name, 'CN')
    expect_equal(length(added_juncs$win), 25)
    ## addJuncs, cn = FALSE
    ## added_juncs_cnfalse = ggnew$addJuncs(junc = junctions, cn = FALSE) ERROR
    ##
    ## addSegs = function(tile)
    added_segs = ggnew$addSegs(tile = test_segs)
    expect_true(is(added_segs, 'gGraph'))
    expect_equal(length(added_segs$segstats), 1124)
    expect_equal(dim(added_segs$edges)[1], 1074)
    expect_equal(dim(added_segs$edges)[2], 14)
    expect_equal(length(added_segs$junctions), 0)
    expect_error(added_segs$G, NA) ## check it works
    expect_equal(length(added_segs$adj), 1263376)
    expect_equal(length(added_segs$A),  1263376)
    ## expect_equal(added_segs$parts, NULL)
    expect_equal(length(added_segs$seqinfo), 25)   
    expect_equal(added_segs$purity, NULL)
    expect_equal(added_segs$ploidy, NULL)  ## checks!
    expect_true(is(added_segs$td, 'gTrack'))
    expect_equal((added_segs$td)$ygap, 2)
    expect_match((added_segs$td)$name, 'CN')
    expect_equal(length(added_segs$win), 25)
    ## karyograph = function(tile=NULL, juncs=NULL, cn=FALSE, regular=FALSE)
    ## default
    ggnew_karyograph = ggnew$karyograph()
    expect_true(is(ggnew_karyograph, 'gGraph'))
    expect_equal(length(ggnew_karyograph$segstats), 1124)
    expect_equal(dim(ggnew_karyograph$edges)[1], 1074)
    expect_equal(dim(ggnew_karyograph$edges)[2], 23)
    expect_equal(length(ggnew_karyograph$junctions), 0)
    expect_error(ggnew_karyograph$G, NA) ## check it works
    expect_equal(length(ggnew_karyograph$adj), 1263376)
    expect_equal(length(ggnew_karyograph$A),  1263376)
    ## expect_equal(length((ggnew_karyograph$parts)$membership), 50)
    ##expect_equal(length((ggnew_karyograph$parts)$csize), 25)
    ##expect_equal(length(ggnew_karyograph$seqinfo), 25)   
    expect_equal(ggnew_karyograph$purity, NULL)
    expect_equal(ggnew_karyograph$ploidy, NULL)  ## checks!
    expect_true(is(ggnew_karyograph$td, 'gTrack'))
    expect_equal((ggnew_karyograph$td)$ygap, 2)
    expect_match((ggnew_karyograph$td)$name, 'CN')
    expect_equal(length(ggnew_karyograph$win), 25)
    ##
    ## simplify = function(mod=TRUE)
    ggnew_simplify = ggnew$simplify()
    expect_true(is(ggnew_simplify, 'gGraph'))
    expect_equal(length(ggnew_simplify$segstats), 50)
    expect_equal(dim(ggnew_simplify$edges)[1], 0)
    expect_equal(dim(ggnew_simplify$edges)[2], 11)
    expect_equal(length(ggnew_simplify$junctions), 0)
    expect_error(ggnew_simplify$G, NA) ## check it works
    expect_equal(length(ggnew_simplify$adj), 2500)
    expect_equal(length(ggnew_simplify$A),  2500)
    ##expect_equal(length((ggnew_simplify$parts)$membership), 50)
    ##expect_equal(length((ggnew_simplify$parts)$csize), 25)
    expect_equal(length(ggnew_simplify$seqinfo), 25)   
    expect_equal(ggnew_simplify$purity, NULL)
    expect_equal(ggnew_simplify$ploidy, NULL)  ## checks!
    expect_true(is(ggnew_simplify$td, 'gTrack'))
    expect_equal((ggnew_simplify$td)$ygap, 2)
    expect_match((ggnew_simplify$td)$name, 'CN')
    expect_equal(length(ggnew_simplify$win), 25)
    ## simplify = function(mod=FALSE)
    ggnew_simplify = ggnew$simplify(mod=FALSE)
    expect_true(is(ggnew_simplify, 'gGraph'))
    expect_equal(length(ggnew_simplify$segstats), 50)
    expect_equal(dim(ggnew_simplify$edges)[1], 0)
    expect_equal(dim(ggnew_simplify$edges)[2], 11)
    expect_equal(length(ggnew_simplify$junctions), 0)
    expect_error(ggnew_simplify$G, NA) ## check it works
    expect_equal(length(ggnew_simplify$adj), 2500)
    expect_equal(length(ggnew_simplify$A),  2500)
    ##expect_equal(length((ggnew_simplify$parts)$membership), 50)
    ##expect_equal(length((ggnew_simplify$parts)$csize), 25)
    expect_equal(length(ggnew_simplify$seqinfo), 25)   
    expect_equal(ggnew_simplify$purity, NULL)
    expect_equal(ggnew_simplify$ploidy, NULL)  ## checks!
    expect_true(is(ggnew_simplify$td, 'gTrack'))
    expect_equal((ggnew_simplify$td)$ygap, 2)
    expect_match((ggnew_simplify$td)$name, 'CN')
    expect_equal(length(ggnew_simplify$win), 25)
    ## decouple = function(mod=TRUE)
    ggnew_decouple = ggnew$decouple()
    expect_true(is(ggnew_decouple, 'gGraph'))
    expect_equal(length(ggnew_decouple$segstats), 50)
    expect_equal(dim(ggnew_decouple$edges)[1], 0)
    expect_equal(dim(ggnew_decouple$edges)[2], 19)
    expect_equal(length(ggnew_decouple$junctions), 0)
    expect_error(ggnew_decouple$G, NA) ## check it works
    expect_equal(length(ggnew_decouple$adj), 2500)
    expect_equal(length(ggnew_decouple$A),  2500)
    ##expect_equal(length((ggnew_decouple$parts)$membership), 50)
    ##expect_equal(length((ggnew_decouple$parts)$csize), 25)
    expect_equal(length(ggnew_decouple$seqinfo), 25)   
    expect_equal(ggnew_decouple$purity, NULL)
    expect_equal(ggnew_decouple$ploidy, NULL)  ## checks!
    expect_true(is(ggnew_decouple$td, 'gTrack'))
    expect_equal((ggnew_decouple$td)$ygap, 2)
    expect_match((ggnew_decouple$td)$name, 'CN')
    expect_equal(length(ggnew_decouple$win), 25)
    ggnew_decouple = ggnew$decouple(mod=FALSE)
    expect_true(is(ggnew_decouple, 'gGraph'))
    expect_equal(length(ggnew_decouple$segstats), 50)
    expect_equal(dim(ggnew_decouple$edges)[1], 0)
    expect_equal(dim(ggnew_decouple$edges)[2], 19)
    expect_equal(length(ggnew_decouple$junctions), 0)
    expect_error(ggnew_decouple$G, NA) ## check it works
    expect_equal(length(ggnew_decouple$adj), 2500)
    expect_equal(length(ggnew_decouple$A),  2500)
    ##expect_equal(length((ggnew_decouple$parts)$membership), 50)
    ##expect_equal(length((ggnew_decouple$parts)$csize), 25)
    expect_equal(length(ggnew_decouple$seqinfo), 25)   
    expect_equal(ggnew_decouple$purity, NULL)
    expect_equal(ggnew_decouple$ploidy, NULL)  ## checks!
    expect_true(is(ggnew_decouple$td, 'gTrack'))
    expect_equal((ggnew_decouple$td)$ygap, 2)
    expect_match((ggnew_decouple$td)$name, 'CN')
    expect_equal(length(ggnew_decouple$win), 25)
    ## 
    ## add = function(gg, mod=FALSE)
    added = ggnew$add(gg=ggnew_karyograph)
    expect_true(is(added, 'gGraph'))
    expect_equal(length(added$segstats), 100)
    expect_equal(dim(added$edges)[1], 0)
    expect_equal(dim(added$edges)[2], 10)
    expect_equal(length(added$junctions), 0)
    expect_error(added$G, NA) ## check it works
    expect_equal(length(added$adj), 10000)
    expect_equal(length(added$A),  10000)
    ##expect_true(is(added$parts, 'list'))
    ##expect_equal(length(added$parts), 3)
    ##expect_equal(length(added$parts$membership), 100)
    ##expect_equal(length((added$parts)$csize), 25)
    expect_equal(length(added$seqinfo), 25)   
    expect_equal(added$purity, NULL)
    expect_equal(added$ploidy, NULL)  ## checks!
    expect_true(is(added$td, 'gTrack'))
    expect_equal((added$td)$ygap, 2)
    expect_match((added$td)$name, 'CN')
    expect_equal(length(added$win), 25)
    ## add() several times
    five_adds = ggnew$add(gg=ggnew_karyograph)$add(gg=ggnew_karyograph)$add(gg=ggnew_karyograph)$add(gg=ggnew_karyograph)$add(gg=ggnew_karyograph)
    expect_true(is(five_adds, 'gGraph'))
    expect_equal(length(five_adds$segstats), 300)
    expect_equal(dim(five_adds$edges)[1], 0)
    expect_equal(dim(five_adds$edges)[2], 10)
    expect_equal(length(five_adds$junctions), 0)
    expect_error(five_adds$G, NA) ## check it works
    expect_equal(length(five_adds$adj),  90000)
    expect_equal(length(five_adds$A),   90000)
    ##expect_true(is(five_adds$parts, 'list'))
    ##expect_equal(length(five_adds$parts), 3)
    ##expect_equal(length(five_adds$parts$membership), 300)
    ##expect_equal(length((five_adds$parts)$csize), 25)
    expect_equal(length(five_adds$seqinfo), 25)   
    expect_equal(five_adds$purity, NULL)
    expect_equal(five_adds$ploidy, NULL)  ## checks!
    expect_true(is(five_adds$td, 'gTrack'))
    expect_equal((five_adds$td)$ygap, 2)
    expect_match((five_adds$td)$name, 'CN')
    expect_equal(length(five_adds$win), 25) 
    ## print() to STDOUT
    expect_error(ggnew$print(), NA)
    ## plot() nothing is returned here, so let's do this for now:
    expect_error(ggnew$plot(), NA)
    ## window()
    expect_equal(length(ggnew$window()), 25)
    ## 
    ## > ggnew$layout()
    ## Error in log(private$segs$cn, 1.4) : 
    ## summary()
    expect_true(is.character(ggnew$summary()))
    ## length()
    expect_equal(ggnew$length(), NULL)
    ##
    expect_true(is(ggnew$gg2td(), 'gTrack'))
    ## JSON
    ## > ggnew$json()
    ## Error in eval(jsub, SDenv, parent.frame()) : object 'y' not found
    ## In addition: Warning message:
    ## In `[.data.table`(node.dt, , `:=`(y, private$segs$cn[oid])) :
    ##   Adding new column 'y' then assigning NULL (deleting it).
    ## HTML
    ## > ggnew$html()
    ## No gGnome.js repository found on your system.
    ## Error in ggnew$html() : Get from https://github.com/mskilab/gGnome.js
    ## gg2j()
    ## > ggnew$gg2js()
    ## Error in eval(jsub, SDenv, parent.frame()) : object 'y' not found
    ## In addition: Warning message:
    ## In `[.data.table`(node.dt, , `:=`(y, private$segs$cn[oid])) :
    ##   Adding new column 'y' then assigning NULL (deleting it).
    ## components()
    component = ggnew$components(mc.cores=2)
    expect_true(is(component, 'list'))
    expect_equal(length(component), 25)
    expect_equal(length(component$segstats), 0)
    expect_equal(dim(component$edges)[1], NULL)
    expect_equal(dim(component$edges)[2], NULL)
    expect_equal(length(component$junctions), 0)
    expect_equal(component$G, NULL) 
    expect_equal(length(component$adj), 0)
    expect_equal(length(component$A),  0)
    expect_equal(length(component$parts), 0)
    expect_equal(length(component$seqinfo), 0)   
    expect_equal(component$purity, NULL)
    expect_equal(component$ploidy, NULL)  ## checks!
    expect_equal(component$td, NULL)
    expect_equal(length(component$win), 0) 
    ##
    ## subgraph = function(v=numeric(0), na.rm=T, mod=T)
    ## default 
    subgraphed = ggnew$subgraph()
    expect_true(is(subgraphed, 'gGraph'))
    expect_equal(length(subgraphed), 25)
    expect_equal(length(subgraphed$segstats), 50)
    expect_equal(dim(subgraphed$edges)[1], 0)
    expect_equal(dim(subgraphed$edges)[2], 19)
    expect_equal(length(subgraphed$junctions), 0)
    expect_error(subgraphed$G, NA) 
    expect_equal(length(subgraphed$adj), 2500)
    expect_equal(length(subgraphed$A),  2500)
    ##expect_equal(length(subgraphed$parts), 0)
    expect_equal(length(subgraphed$seqinfo), 25)   
    expect_equal(component$purity, NULL)
    expect_equal(component$ploidy, NULL)  ## checks!
    expect_equal(component$td, NULL)
    expect_equal(length(component$win), 0) 
    ##
    ## vertices5K = ggnew$subgraph(v=5000)
    ##
    filled = ggnew$fillin()
    expect_true(is(filled, 'gGraph'))
    expect_equal(length(filled$segstats), 50)
    expect_equal(dim(filled$edges)[1], 0)
    expect_equal(dim(filled$edges)[2], 19)
    expect_equal(length(filled$junctions), 0)
    expect_error(filled$G, NA) 
    expect_equal(length(filled$adj), 2500)
    expect_equal(length(filled$A),  2500)
    ##expect_equal(length(filled$parts), 3)
    expect_equal(length(filled$seqinfo), 25)   
    expect_equal(filled$purity, NULL)
    expect_equal(filled$ploidy, NULL)  ## checks!
    expect_true(is(filled$td, 'gTrack'))
    expect_equal(length(filled$win), 25)
    ## 
    ## trim = function(gr=NULL, mod=FALSE)
    ## default
    trimmed = ggnew$trim()
    expect_true(is(trimmed, 'gGraph'))
    expect_equal(length(trimmed$segstats), 50)
    expect_equal(dim(trimmed$edges)[1], 0)
    expect_equal(dim(trimmed$edges)[2], 19)
    expect_equal(length(trimmed$junctions), 0)
    expect_error(trimmed$G, NA) 
    expect_equal(length(trimmed$adj), 2500)
    expect_equal(length(trimmed$A),  2500)
    expect_equal(length(trimmed$parts), 3)
    expect_equal(length(trimmed$seqinfo), 25)   
    expect_equal(trimmed$purity, NULL)
    expect_equal(trimmed$ploidy, NULL)  ## checks!
    expect_true(is(trimmed$td, 'gTrack'))
    expect_equal(length(trimmed$win), 25)
    ## 
    ## trimmed_mod = ggnew$trim(gr=gr2, mod=TRUE) 
    ## gotg
    ## default
    gotg = ggnew$get.g() 
    expect_true(is(gotg, 'gGraph'))
    expect_equal(length(gotg$segstats), 50)
    expect_equal(dim(gotg$edges)[1], 0)
    expect_equal(dim(gotg$edges)[2], 19)
    expect_equal(length(gotg$junctions), 0)
    expect_error(gotg$G, NA) 
    expect_equal(length(gotg$adj), 2500)
    expect_equal(length(gotg$A),  2500)
    ##expect_equal(length(gotg$parts), 3)
    expect_equal(length(gotg$seqinfo), 25)   
    expect_equal(gotg$purity, NULL)
    expect_equal(gotg$ploidy, NULL)  ## checks!
    expect_true(is(gotg$td, 'gTrack'))
    expect_equal(length(gotg$win), 25)
    ## get.g = function(force=FALSE)
    gotg_forced = ggnew$get.g(force=TRUE) 
    expect_true(is(gotg_forced, 'gGraph'))
    expect_equal(length(gotg_forced$segstats), 50)
    expect_equal(dim(gotg_forced$edges)[1], 0)
    expect_equal(dim(gotg_forced$edges)[2], 19)
    expect_equal(length(gotg_forced$junctions), 0)
    expect_error(gotg_forced$G, NA) 
    expect_equal(length(gotg_forced$adj), 2500)
    expect_equal(length(gotg_forced$A),  2500)
    ##expect_equal(length(gotg_forced$parts), 3)
    expect_equal(length(gotg_forced$seqinfo), 25)   
    expect_equal(gotg_forced$purity, NULL)
    expect_equal(gotg_forced$ploidy, NULL)  ## checks!
    expect_true(is(gotg_forced$td, 'gTrack'))
    expect_equal(length(gotg_forced$win), 25)
    ## with added
    gotgadd = added$get.g(force=TRUE) 
    expect_true(is(gotgadd, 'gGraph'))
    expect_equal(length(gotgadd$segstats), 100)
    expect_equal(dim(gotgadd$edges)[1], 0)
    expect_equal(dim(gotgadd$edges)[2], 19)
    expect_equal(length(gotgadd$junctions), 0)
    expect_error(gotgadd$G, NA) 
    expect_equal(length(gotgadd$adj), 10000)
    expect_equal(length(gotgadd$A),  10000)
    expect_equal(length(gotgadd$parts), 3)  
    expect_equal(length((gotgadd$parts)$membership), 100)
    expect_equal(length((gotgadd$parts)$csize), 25)
    expect_equal(length((gotgadd$parts)$no), 1)  ## 25
    expect_equal(length(gotgadd$seqinfo), 25)   
    expect_equal(gotgadd$purity, NULL)
    expect_equal(gotgadd$ploidy, NULL)  ## checks!
    expect_true(is(gotgadd$td, 'gTrack'))
    expect_equal((gotgadd$td)$ygap, 2)
    expect_match((gotgadd$td)$name, 'CN')
    expect_equal(length(gotgadd$win), 25) 
    ##  hood = function(win, d=NULL, k=NULL, pad=0, bagel=FALSE, ignore.strand=T, verbose=FALSE)  
    ##gr2_win = ggnew$hood(win=gr2)
    expect_error(ggnew$hood(win=grl1)) ## Error in .local(x, y, ...) : setdiff() between a GRanges and a GRangesList object is not supported
    ##grl2_win = ggnew$hood(win=grl.unlist(grl2))
    ##hgseq_win = ggnew$hood(win = si2gr(hg_seqlengths()))
    ##expect_true(is(hgseq_win, 'gGraph'))
    ##expect_equal(length(hgseq_win$segstats), 0)
    ##expect_equal(dim(hgseq_win$edges)[1], 0)
    ##expect_equal(dim(hgseq_win$edges)[2], 19)
    ##expect_equal(length(hgseq_win$junctions), 0)
    ##expect_error(hgseq_win$G, NA) 
    ##expect_equal(length(hgseq_win$adj), 0)
    ##expect_equal(length(hgseq_win$A),  0)
    ##expect_equal(length(hgseq_win$parts), 0)   
    ##expect_equal(hgseq_win$purity, NULL)
    ##expect_equal(hgseq_win$ploidy, NULL)  ## checks!
    ##expect_equal(length(hgseq_win$td), 0)
    ##expect_equal(length(hgseq_win$win), 0) 

})


### try gGraph with inputs above
### gGraph$new(tile = test_segs)

## segstats, edges, grl, td, path, values




### bGraph, default
### should be similar to gGraph

















##test_that('rev.comp works', {
##
##    expect_error(rev.comp())
##    expect_error(rev.comp(data.frame()))
##    gr3 = dt2gr(dt)
##    expect_error(rev.comp(gr3))   ## Error in rev.comp(gr3) : Input must be all strand specific.
##    expect_equal(width(rev.comp(as.vector(gr)[1]), 4)
##    expect_equal(as.character(strand(rev.comp(gr)[1])), "+")
##    expect_equal(width(rev.comp(gr)[2]), 3)
##    expect_equal(as.character(strand(rev.comp(gr)[2])), "+")
##
##})





test_that('capitalize works', {

    str1 = "Foo FOO"
    str2 = "2Foo . $%@"
    str3 = "foobar foo"
    expect_match(capitalize(str1), "Foo FOO")
    expect_match(capitalize(str1, un=TRUE), "foo FOO")
    expect_equal(capitalize(str2), "2Foo . $%@")
    expect_equal(capitalize(str2, un=TRUE), "2Foo . $%@")    ## probably should report this bug to r-lib
    expect_match(capitalize(str3), "Foobar foo")
    expect_match(capitalize(str3, un=TRUE), "foobar foo")

})


  


test_that('ul works', {

    A = matrix(c(2, 4, 3, 1, 5, 7),  nrow=2, ncol=3, byrow = TRUE)    
    expect_equal(ul(A, n=0), NULL)   ### Is this expected behavior? 
    expect_equal(as.integer(ul(A, n=1)), 2)
    expect_equal(dim(ul(A, n=2))[1], 2)
    expect_equal(dim(ul(A, n=2))[2], 2)
    expect_equal(dim(ul(A, n=999))[1], 2)
    expect_equal(dim(ul(A, n=9999))[2], 2)   ### Is this expected behavior? 

})






test_that('e2j works', {

    expect_equal(length(e2j(test_segs, test_es)), 2)

})




test_that('etype works', {

    expect_equal(dim(etype(test_segs, test_es))[1], 12)
    expect_equal(dim(etype(test_segs, test_es))[2], 16) 

})







test_that('write.tab() works', {

    expect_equal(length(write.tab(dt)), 0)  ### throws dt to STDOUT

})







test_that('get.ploidy works', {

    expect_error(get.ploidy(GRangesList()))
    expect_equal(get.ploidy(test_segs), 3)
    ## if (length(cnix <- grep("CN", colnames(mcols(segs)), ignore.case = T)) == 

})





test_that('dedup() works', {

    expect_equal(dedup(c(rep(2, 10.5), rep(3, 20)))[30], "3.20")

})










test_that('chr2num works', {

    expect_equal(as.logical(chr2num("ChrX")), NA)
    expect_equal(chr2num("chrX"), 23)
    expect_equal(chr2num("chrY"), 24)

})



test_that('affine.map works', {

    expect_equal(affine.map(49), 0.5)

})


test_that('gr.flatmap works', {

    expect_equal((gr.flatmap(example_genes, windows=GRanges('1:10000-20000'))$window.segs)$start, 1)
    expect_equal((gr.flatmap(example_genes, windows=GRanges('1:10000-20000'))$window.segs)$end, 10001)
    expect_equal(length((gr.flatmap(example_genes, windows=GRanges('1:10000-20000'))$window.segs)$grl.segs), 0)

})




### XT's tests

##-------------------------------------------------------##
test_that('constructors and essential functions', {
    ## small example, nested tDUP
    ## default
    expect_equal(dim(etype(test_segs, test_es))[1], 12)
    expect_equal(dim(etype(test_segs, test_es))[2], 16)
    expect_equal(unique(as.integer(etype(test_segs, test_es)$toChr)), 5)
    expect_equal(any(etype(test_segs, test_es)$fromLoose), FALSE)
    expect_equal(any(etype(test_segs, test_es)$toLoose), FALSE)
    expect_error(etype(GRangesList(), GRangesList()))  ## Error in etype(GRangesList(), GRangesList()) : Error:segs must be GRanges
    expect_error(etype(GRanges(), GRangesList()))      ## Error in etype(GRanges(), GRangesList()) : Error:es must be data.frame
    expect_error(etype(GRanges(), data.table()))       ## Error: 'from' & 'to' must be in es!
    expect_error(gGraph$new(), NA)  ## test it works
    foo = gGraph$new(segs=test_segs, es=test_es)
    ##expect_equal(dim(foo$edges)[1], 12)
    ##expect_equal(dim(foo$edges)[2], 16)
    expect_equal(max((foo$edges)$cn), 3)
    expect_equal(max((foo$edges)$fromStart), 18593415)
    expect_equal(max((foo$edges)$fromEnd), 18793414)
    foo = gGraph$new(segs=test_segs, es=test_es)
    ##expect_equal(dim(foo$edges)[1], 12)
    ## expect_equal(dim(foo$edges)[2], 16)
    expect_equal(max((foo$edges)$cn), 3)
    expect_equal(max((foo$edges)$fromStart), 18593415)
    expect_equal(max((foo$edges)$fromEnd), 18793414)
    expect_equal(dim((foo$nullGGraph())$edges)[1], 0)
    expect_equal(dim((foo$nullGGraph())$edges)[2], 3)

})




##-------------------------------------------------------##
test_that('gGraph, dipGraph', {

    expect_error(gGraph$new()$dipGraph(), NA)
    expect_equal(nrow(gGraph$new()$dipGraph()$edges), 0)
    expect_equal(length(gGraph$new()$dipGraph()$segstats), length(gUtils::hg_seqlengths())*2)

})



test_that('karyograph', {

    kag.tile = gGraph$new(tile = test_segs)
    expect_true(inherits(kag.tile, "gGraph"))

})




##-------------------------------------------------------##
test_that('gread', {

    jab = system.file('extdata', 'jabba.simple.rds', package="gGnome")
    message("JaBbA result: ", jab)
    prego = system.file('extdata', 'intervalFile.results', package='gGnome')
    message("PREGO results: ", prego)
    weaver = system.file('extdata', 'weaver', package='gGnome')
    message("Weaver results: ", weaver)
    expect_error(gread('no_file_here'))
    jab_bgraph = gread(jab)
    expect_true(is(jab_bgraph, "bGraph"))
    ## preg_bgraph = gread(prego)
    ## expect_true(is(preg_bgraph, "bGraph")) ### 'gGraph'
    ## wv_bgraph = gread(weaver)   
    ## expect_true(is(wv_bgraph, "bGraph"))  ### 'gGraph'
    ## if (is.list(file)){
    list_foo = gread(readRDS(system.file("extdata", "jabba.simple.rds", package="gGnome")))
    expect_true(is(list_foo, 'bGraph'))
    
})





## ##-------------------------------------------------------##
## test_that('gtf2json', {
##     expect_error(gread('no_file_here'))
##     expect_equal(gtf2json(system.file('extdata', 'test.gtf', package='gGnome')), "./gtf.json")
##     system(paste('rm', "./gtf.json"))
## })




## * could not find function "setxor"
##-------------------------------------------------------##

test_that('setxor', {

    A = c(1, 2, 3)
    B = c(1, 4, 5)
    expect_equal(setxor(A, B), c(2, 3, 4, 5))

})





##-------------------------------------------------------##
test_that('special ranges functions for skew-symmetric graph', {


    segments = readRDS(jab)$segstats
    junctions = readRDS(jab)$junctions
    expect_equal(length(seg.fill(GRanges())), 0)
    expect_equal(length(seg.fill(segments)), 2346)
    ## check 'verbose'
    expect_equal(length(seg.fill(segments, verbose=TRUE)), 2346)
    expect_equal(length(seg.fill(segments %Q% (strand=="+"), verbose=TRUE)), 2346)
    expect_equal(dim(hydrogenBonds(segments))[1], length(segments))
    expect_equal(dim(hydrogenBonds(segments))[2], 3)
    expect_equal(unique(hydrogenBonds(segments)$type), 'hydrogen')

})


## Error: Test failed: 'gWalks'
## * length(gw <<- as(grl, "gWalks")) not equal to sum(values(grl)$cn > 0).
## 1/1 mismatches
## [1] 32 - 630 == -598
## * `bg <<- as(gw, "bGraph")` threw an error.
## Message: Error: Given edge data is not skew-symmetric!!!
## Class:   simpleError/error/condition
## * object 'bg' not found
## 1: expect_equal(length(bg$junctions), sum(values(junctions)$cn > 0)) at :12
## 2: quasi_label(enquo(object), label)
## 3: eval_bare(get_expr(quo), get_env(quo))
 



##-------------------------------------------------------##
##test_that('gWalks', {
##
##    jab = system.file('extdata', 'jabba.simple.rds', package="gGnome")
##    message("JaBbA result: ", jab)
##    segments = readRDS(jab)$segstats
##    junctions = readRDS(jab)$junctions
##    grl = system.file("extdata", "gw.grl.rds", package="gGnome")
##    message("Walks for testing:", grl)
##    grl = readRDS(grl)
##    expect_equal(length(gw <<- as(grl, "gWalks")), sum(values(grl)$cn>0))
##    expect_error(bg <<- as(gw, "bGraph"), NA)
##    expect_equal(length(bg$junctions), sum(values(junctions)$cn>0))
##    expect_true(inherits(gw.simp <<- gw$simplify(mod=FALSE), "gWalks"))
##    expect_error(bg.simp <<- as(gw.simp, "bGraph"), NA)
##    expect_error(bg.dc <<- bg.simp$decouple(mod=FALSE), NA)
##    ## expect_equal(length(bg.dc$junctions), length(bg$junctions)) 291>269
##    ## why does simplifying gwalks then decouple create more junctions????
##
##})


## I think downloading data without warnings is evil

## trying URL 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/GRCh37_mapping/gencode.v27lift37.basic.annotation.gff3.gz'
## Content type 'unknown' length 39550148 bytes (37.7 MB)
## ==================================================
## trying URL 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/GRCh37_mapping/gencode.v27lift37.basic.annotation.gff3.gz'
## Content type 'unknown' length 39550148 bytes (37.7 MB)
## ==================================================



##-------------------------------------------------------##
##test_that('fusions', {
##    juncs = system.file('extdata', 'testing_junctions.rds', package="gGnome")
##    message("Junctions for testing: ", juncs)
##    juncs = readRDS(juncs)

    ## make sure the gene annotation can be loaded
##    expect_error(cds <<- read_gencode(type = "cds"), NA)
##    expect_error(fusions())
##    expect_error(fusions(junc = juncs, cds = cds), NA) ## no problem
##})

## ##-------------------------------------------------------##
## test_that('graph distance and proximity', {
##     query = readRDS()
##     expect_error()
## })


