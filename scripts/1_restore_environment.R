a <- commandArgs(trailingOnly=TRUE)
resume <- length(a)==4 && a[4]=='--resume'
if(!(length(a)==3 || resume)) stop('Usage: 1_restore_environment.R LOCKFILE WORKDIR LIBRARY [--resume]')
stopifnot(requireNamespace('renv',quietly=TRUE),file.exists(a[1]))
if(!resume) stopifnot(!dir.exists(a[2]),!dir.exists(a[3]))
if(resume) stopifnot(dir.exists(a[2]),dir.exists(a[3]))
lock <- normalizePath(a[1])
if(!resume) {dir.create(a[2],recursive=TRUE); dir.create(a[3],recursive=TRUE)}
work <- normalizePath(a[2]); lib <- normalizePath(a[3])
.libPaths(lib,include.site=FALSE)
Sys.setenv(R_LIBS=lib,R_LIBS_USER=lib,R_LIBS_SITE=lib)
Sys.setenv(RENV_PATHS_ROOT=file.path(work,'root'),RENV_PATHS_CACHE=file.path(work,'cache'),
           RENV_PATHS_SOURCE=file.path(work,'source'),RENV_PATHS_SANDBOX=file.path(work,'base-sandbox'),RENV_CONFIG_CACHE_ENABLED='FALSE')
options(repos=c(CRAN='https://cran.r-project.org'),timeout=300)
stopifnot(as.character(utils::packageVersion('renv'))=='1.3.1')
sandbox_result<-get('renv_sandbox_activate',asNamespace('renv'))(project=work)
if(inherits(sandbox_result,'error'))stop(sandbox_result)
.libPaths(lib,include.site=FALSE)
writeLines(.libPaths(),file.path(work,'restore_library_paths.txt'))
renv::restore(project=work,library=lib,lockfile=lock,prompt=FALSE,transactional=FALSE)
expected <- renv::lockfile_read(lock)$Packages
actual <- installed.packages(lib.loc=c(lib,.Library))
check <- data.frame(package=names(expected),expected=vapply(expected,function(x)x$Version,''))
check$actual <- actual[match(check$package,rownames(actual)),'Version']
check$library <- actual[match(check$package,rownames(actual)),'LibPath']
check$match <- !is.na(check$actual) & check$expected==check$actual
write.csv(check,file.path(work,'version_verification.csv'),row.names=FALSE)
stopifnot(all(check$match))
writeLines(capture.output(sessionInfo()),file.path(work,'restore_session.txt'))
cat('RESTORE_VERIFIED:',nrow(check),'packages\n')
