
test_that("load_start adds the input parameters and variables to gc_env", {
    load_start("SF")
    expect_true(all(
    c("n.i", "aging.rate.m", "theta", "diag.rate", "births", "out_all", 
    "n.total", "p.low", "s.dist.dat.f", "diag.index", "y.f.23", "omega.t", 
    "index", "s.index", "s.dist.dat.m", "diag.subpop.rate", "o.f", 
    "cal.end", "y.index", "p.sa.f.y", "prev.nhanes", "diag.rate.sd", 
    "p.dist", "o.m", "y.f", "rr.diag.subpop.sd", "dat.s.dist", "n.s.pop.sa", 
    "site", "p.symp.dat", "o.f.1", "n.s.dist.sa", "p.low.msm", "y.m", 
    "o.f.2", "var.p.msm.ssun", "o.f.3", "n.s.pop", "n.dist", "prior.param1", 
    "n.sa", "diag.age.sex.rate", "prior.param2", "var.symp.ssun", 
    "z.index", "cal.start", "start.year", "aging", "p.sa.f.o", "y.m.1", 
    "y.m.2", "diag.rr", "msw", "y.m.3", "p.k", "inc.index", "y.m.4", 
    "n.s.a", "m1", "m2", "m3", "prev.extra.cat.dat", "m4", "f1", 
    "births.sa", "f2", "f3", "age.dist.dat", "rr.diag.subpop", "f4", 
    "priors", "y.m.msw", "p.s.dist", "model.end", "var.s.dist.dat", 
    "end.year", "tstep", "pop1", "p.sa.f", "pop2", "males", "pop3", 
    "s.dist.sd", "pop4", "sd.prior", "females", "yinit", "p.sa.m", 
    "p.symp.ssun", "o.f.23", "births.nsa", "nsa.index", "p.msm.ssun", 
    "ind", "p.sa.m.y", "aging.nsa", "p.msm", "o.m.1", "o.m.2", "o.m.3", 
    "cal.period", "p.msm.dat", "o.m.4", "n.dist.sa", "o.m.msw", "age.cat", 
    "n.s.dist", "n.ns.a", "symp.index", "i", "j", "var.prev.nhanes", 
    "k", "l", "n.nsa", "omega", "aging.rate.f", "p.sa.m.o", "y.f.1", 
    "p.sa", "p.sa.array", "p.s.1", "init.Y", "y.f.2", "p.s.2", "y.f.3", 
    "p.s.3") %in% names(gc_env)))


})


